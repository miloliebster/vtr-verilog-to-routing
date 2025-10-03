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
    std::string port_name; // TODO Just using StringId would be preferable, but access to this would need to be exposed from netlist class. Friend class?
    BitIndex pin_port_bit;

    bool operator==(PackSignatureConnection const& rhs) const {
        if (this->pb_graph_node != rhs.pb_graph_node) return false;
        if (this->pin_port_bit != rhs.pin_port_bit) return false;
        if (this->port_name != rhs.port_name) return false;
        return true;
    }

    bool operator<(PackSignatureConnection const& rhs) const {
        if (this->pb_graph_node >= rhs.pb_graph_node) return false;
        if (this->pin_port_bit >= rhs.pin_port_bit) return false;
        if (this->port_name.compare(rhs.port_name) >= 0) return false;
        return true;
    }

    // For characterization -- removable later
    std::string to_string() {
        return std::format("<<\"{}\": {:#08x}>, <\"{}\", {}>>", pb_graph_node->pb_type->name, reinterpret_cast<uintptr_t>(pb_graph_node), port_name, pin_port_bit);
    }
};

struct PackSignaturePrimitive {
    const t_pb_graph_node* pb_graph_node;
    size_t num_incoming_sources; // TODO this needs to be on a pin-by-pin basis
    size_t num_sinks; // TODO this needs to be on a pin-by-pin basis
    std::vector<PackSignatureConnection> intracluster_sources_to_primitive_inputs;
    std::vector<PackSignatureConnection> intracluster_sinks_of_primitive_outputs;
    AtomBlockId atom_block_id; // not part of signature per se, but referenced when backtracking

    PackSignaturePrimitive() {}

    PackSignaturePrimitive(PackSignaturePrimitive& other) {
        this->pb_graph_node = other.pb_graph_node;
        this->num_incoming_sources = other.num_incoming_sources;
        this->num_sinks = other.num_sinks;
        this->intracluster_sources_to_primitive_inputs = other.intracluster_sources_to_primitive_inputs;
        this->intracluster_sinks_of_primitive_outputs = other.intracluster_sinks_of_primitive_outputs;
    }

    bool operator==(PackSignaturePrimitive const& rhs) const {
        if (this->pb_graph_node != rhs.pb_graph_node) return false;
        if (this->num_incoming_sources != rhs.num_incoming_sources) return false;
        if (this->num_sinks != rhs.num_sinks) return false;
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
    std::vector<LegalizationClusterId> legalization_cluster_ids;
    PackSignaturePrimitive* primitive = nullptr;
    PackSignatureTreeNode* parent = nullptr;

    std::vector<PackSignaturePrimitive*> child_primitives;
    std::vector<PackSignatureTreeNode*> child_nodes;

    PackSignatureTreeNode() {}
    ~PackSignatureTreeNode() {
        for (PackSignatureTreeNode* child_node : child_nodes) delete child_node;
        for (PackSignaturePrimitive* child_primitive : child_primitives) delete child_primitive;
    }
};

class PackSignatureTree {
public:
    PackSignatureTree() {}

    ~PackSignatureTree() {
        for (PackSignatureTreeNode* signature_root_node : signatures_) delete signature_root_node;
    }

    void start_pack_signature(const t_logical_block_type* cluster_logical_block_type);
    void add_primitive(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id);

    inline void finalize_path(LegalizationClusterId legalization_cluster_id) {
        at_node_->legalization_cluster_ids.push_back(legalization_cluster_id);
    }

    void log_equivalent(); // XXX characterization only

private:
    PackSignaturePrimitive* generate_primitive(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id);

    // signatures_[i] corresponds to the logical_block_type pointed to by cluster_logical_block_types_[i].
    std::vector<const t_logical_block_type*> cluster_logical_block_types_;
    std::vector<PackSignatureTreeNode*> signatures_;

    // Current node for signature being constructed
    PackSignatureTreeNode* at_node_;
};

extern PackSignatureTree g_pack_signatures; // XXX this should not be a global in the end. Part of ClusterLegalizer?
