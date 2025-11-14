#pragma once

#include <iostream>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <string>
#include <vector>
#include <utility>

#include "atom_netlist.h"
#include "physical_types.h"
#include "globals.h"
#include "cluster_legalizer.h"

void try_open_logfile();

extern std::ofstream g_logfile;

// pb_graph_node, port_name, and pin_port_bit uniquely identify a connection inside a cluster
// TODO: Need to determine if primitive has input port equivalence and handle accordingly
struct PackSignatureConnection {
    const t_pb_graph_node* pb_graph_node;
    std::string source_port_name; // XXX it would be best to use StringID here to save space
    BitIndex source_bit_index;
    std::string sink_port_name;
    BitIndex sink_bit_index;

    bool operator==(PackSignatureConnection const& rhs) const {
        if (this->pb_graph_node != rhs.pb_graph_node) return false;
        if (this->source_bit_index != rhs.source_bit_index) return false;
        if (this->sink_bit_index != rhs.sink_bit_index) return false;
        if (this->source_port_name != rhs.source_port_name) return false;
        return this->sink_port_name == rhs.sink_port_name;
    }

    bool operator<(PackSignatureConnection const& rhs) const {
        if (this->pb_graph_node != rhs.pb_graph_node) return this->pb_graph_node < rhs.pb_graph_node;
        if (this->source_bit_index != rhs.source_bit_index) return this->source_bit_index < rhs.source_bit_index;
        if (this->sink_bit_index != rhs.sink_bit_index) return this->sink_bit_index < rhs.sink_bit_index;
        if (this->source_port_name != rhs.source_port_name) return this->source_port_name < rhs.source_port_name;
        return this->sink_port_name < rhs.sink_port_name;
    }

    // XXX For characterization -- removable later
    std::string to_string() {
        return std::format("<\"{}\": {:#08x}, <\"{}\", {}>, <\"{}\", {}>>",
                           pb_graph_node->pb_type->name,
                           reinterpret_cast<uintptr_t>(pb_graph_node),
                           source_port_name,
                           source_bit_index,
                           sink_port_name,
                           sink_bit_index);
    }
};

struct PackSignatureExternalIO {
    std::vector<std::vector<PackSignatureConnection>> cluster_inputs;
    std::vector<PackSignatureConnection> cluster_outputs;

    // identify the legalization cluster IDs that terminate at this node
    // XXX is this required for anything besides stats?
    std::vector<LegalizationClusterId> successful_legalization_cluster_ids;
    std::vector<bool> successful_legalization_cluster_detailedness;

    std::vector<LegalizationClusterId> failed_legalization_cluster_ids;
    std::vector<bool> failed_legalization_cluster_detailedness;

    bool operator==(PackSignatureExternalIO const& rhs) const {
        if (this->cluster_inputs.size() != rhs.cluster_inputs.size()) return false;
        if (this->cluster_outputs.size() != rhs.cluster_outputs.size()) return false;
        for (size_t i = 0; i < this->cluster_inputs.size(); i++) {
            if (this->cluster_inputs[i].size() != rhs.cluster_inputs[i].size()) return false;
            for (size_t j = 0; j < this->cluster_inputs[i].size(); j++) {
                if (this->cluster_inputs[i][j] != rhs.cluster_inputs[i][j]) return false;
            }
        }
        for (size_t i = 0; i < this->cluster_outputs.size(); i++) {
            if (this->cluster_outputs[i] != rhs.cluster_outputs[i]) return false;
        }
        return true;
    }
};

// TODO: PackSignaturePrimitive and PackSignatureTreeNode should be able to be refactored into one struct if the tree roots are handled different.
//       Is that desirable?
struct PackSignaturePrimitive {
    const t_pb_graph_node* pb_graph_node;
    std::vector<PackSignatureConnection> intracluster_sources_to_primitive_inputs;
    std::vector<PackSignatureConnection> intracluster_sinks_of_primitive_outputs;
    AtomBlockId atom_block_id; // not part of signature, but referenced when backtracking

    PackSignaturePrimitive() {}

    PackSignaturePrimitive(PackSignaturePrimitive& other) {
        this->pb_graph_node = other.pb_graph_node;
        this->intracluster_sources_to_primitive_inputs = other.intracluster_sources_to_primitive_inputs;
        this->intracluster_sinks_of_primitive_outputs = other.intracluster_sinks_of_primitive_outputs;
    }

    bool operator==(PackSignaturePrimitive const& rhs) const {
        if (this->pb_graph_node != rhs.pb_graph_node) return false;
        if (this->intracluster_sources_to_primitive_inputs.size() != rhs.intracluster_sources_to_primitive_inputs.size()) return false;
        if (this->intracluster_sinks_of_primitive_outputs.size() != rhs.intracluster_sinks_of_primitive_outputs.size()) return false;
        for (size_t i = 0; i < this->intracluster_sources_to_primitive_inputs.size(); i++) {
            if (this->intracluster_sources_to_primitive_inputs[i] != rhs.intracluster_sources_to_primitive_inputs[i]) return false;
        }
        for (size_t i = 0; i < this->intracluster_sinks_of_primitive_outputs.size(); i++) {
            if (this->intracluster_sinks_of_primitive_outputs[i] != rhs.intracluster_sinks_of_primitive_outputs[i]) return false;
        }
        return true;
    }
};

struct PackSignatureTreeNode {
    size_t visits = 0; // XXX Only useful for characterization. Can be removed for final implementation.

    PackSignaturePrimitive* primitive = nullptr;
    PackSignatureTreeNode* parent = nullptr;

    std::vector<PackSignaturePrimitive*> child_primitives;
    std::vector<PackSignatureTreeNode*> child_nodes;

    std::vector<PackSignatureExternalIO*> leaf_nodes;

    PackSignatureTreeNode() {}
    ~PackSignatureTreeNode() {
        for (PackSignatureTreeNode* child_node : child_nodes) delete child_node;
        for (PackSignaturePrimitive* child_primitive : child_primitives) delete child_primitive;
        for (PackSignatureExternalIO* leaf_node : leaf_nodes) delete leaf_node;
    }
};


class PackSignatureTree {
public:
    bool detailed_legalization = false;
    std::chrono::duration<double> speculative_legalization_success_duration;
    std::chrono::duration<double> speculative_legalization_failure_duration;
    std::chrono::duration<double> detailed_legalization_success_duration;
    std::chrono::duration<double> detailed_legalization_failure_duration;

    size_t memory_usage_scratch = 0;
    size_t total_memory_used = 0;

    PackSignatureTree() {}

    ~PackSignatureTree() {
        for (PackSignatureTreeNode* signature_root_node : signatures_) delete signature_root_node;
    }

    void start_pack_signature(const t_logical_block_type* cluster_logical_block_type);
    void add_primitive(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id);

    void finalize_path(LegalizationClusterId legalization_cluster_id);
    void fail_path(LegalizationClusterId legalization_cluster_id);

    void log_equivalent(); // XXX characterization only

private:
    PackSignaturePrimitive* generate_primitive(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id);

    // signatures_[i] corresponds to the logical_block_type pointed to by cluster_logical_block_types_[i].
    std::vector<const t_logical_block_type*> cluster_logical_block_types_;
    std::vector<PackSignatureTreeNode*> signatures_;

    // Current node for signature being constructed
    PackSignatureTreeNode* at_node_;

    // External IO bookkeeping
    struct ExternalOutputRecord {
        PackSignatureConnection connection;
        size_t external_sinks_count;
    };

    std::map<AtomNetId, std::vector<PackSignatureConnection>> input_nets_;
    std::map<AtomNetId, ExternalOutputRecord> output_nets_;
};

extern PackSignatureTree g_pack_signatures; // XXX this should not be a global in the end. Part of ClusterLegalizer?
