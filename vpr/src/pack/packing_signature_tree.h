#pragma once

#include <iostream>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "vtr_strong_id.h"
#include "atom_netlist.h"
#include "physical_types.h"
#include "globals.h"
#include "cluster_legalizer.h"

// Define to resolve circular dependancy with cluster_legalizer.h
typedef vtr::StrongId<struct legalization_cluster_id_tag, size_t> LegalizationClusterId;

typedef std::tuple<int, std::string, BitIndex> PstPin;

struct PstPinHash {
    size_t operator()(PstPin pin) const noexcept {
        return std::hash<int>{}(std::get<0>(pin)) ^ std::hash<std::string>{}(std::get<1>(pin)) ^ std::hash<int>{}((int)std::get<2>(pin));
    }
};

typedef vtr::StrongId<struct pst_pin_id_tag> PstPinId;

struct PstConnection {
    PstPinId source;
    PstPinId sink;

    bool operator==(PstConnection const& rhs) const {
        return (this->source == rhs.source && this->sink == rhs.sink);
    }

    bool operator<(PstConnection const& rhs) const {
        return (this->source != rhs.source) ? this->source < rhs.source : this->sink < rhs.sink;
    }

    std::string to_string() const {
        return "<pin" + std::to_string((int)source) + ", pin" + std::to_string((int)sink) + ">";
    }
};

struct ExternalConnectivityNode {
    std::vector<std::vector<PstPinId>> cluster_inputs;
    std::vector<PstPinId> cluster_outputs;

    bool legal;

    int successful_clusters;// XXX
    int successful_detailed_clusters;// XXX

    int failed_clusters;// XXX
    int failed_detailed_clusters;// XXX

    ExternalConnectivityNode()
        : legal(false), successful_clusters(0), successful_detailed_clusters(0), failed_clusters(0), failed_detailed_clusters(0) {}

    bool operator==(ExternalConnectivityNode const& rhs) const {
        if (this->cluster_inputs != rhs.cluster_inputs) return false;
        return (this->cluster_outputs == rhs.cluster_outputs);
    }
};

struct PrimitiveSignatureNode {
    int primitive_num;
    std::vector<PstConnection> intracluster_sources_to_primitive_inputs;
    std::vector<PstConnection> intracluster_sinks_of_primitive_outputs;

    std::vector<PrimitiveSignatureNode*> child_psn;
    std::vector<ExternalConnectivityNode*> child_ecn;

    PrimitiveSignatureNode() : primitive_num(-1) {}
    ~PrimitiveSignatureNode() {
        for (auto psn : child_psn) delete psn;
        for (auto ecn : child_ecn) delete ecn;
    }

    bool operator==(PrimitiveSignatureNode const& rhs) const {
        if (this->primitive_num != rhs.primitive_num) return false;
        if (this->intracluster_sources_to_primitive_inputs != rhs.intracluster_sources_to_primitive_inputs) return false;
        if (this->intracluster_sinks_of_primitive_outputs != rhs.intracluster_sinks_of_primitive_outputs) return false;
        return true;
    }
};

enum e_packing_signature_legality {
    UNKNOWN,
    LEGAL,
    ILLEGAL
};

class PackingSignatureTree {
public:
    bool detailed_legalization = false; // XXX
    std::chrono::duration<double> speculative_legalization_success_duration; // XXX
    std::chrono::duration<double> speculative_legalization_failure_duration; // XXX
    std::chrono::duration<double> detailed_legalization_success_duration; // XXX
    std::chrono::duration<double> detailed_legalization_failure_duration; // XXX
    std::chrono::duration<double> signature_processing_duration; // XXX

    bool routed = false;

    PackingSignatureTree() {}
    ~PackingSignatureTree() {
        for (auto psn : packing_signatures_) delete psn;
    }

    void start_packing_signature(const t_logical_block_type* cluster_logical_block_type);
    void add_psn(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id);
    void set_checkpoint();
    void rollback_to_checkpoint();

    e_packing_signature_legality check_legality();
    void mark_signature_as_legal(LegalizationClusterId legalization_cluster_id = LegalizationClusterId::INVALID());
    void mark_signature_as_illegal(LegalizationClusterId legalization_cluster_id = LegalizationClusterId::INVALID());
    void finalize_path(LegalizationClusterId legalization_cluster_id);
    void fail_path(LegalizationClusterId legalization_cluster_id);

    void log_equivalent(); // XXX characterization only

private:
    PrimitiveSignatureNode* create_psn(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id);
    ExternalConnectivityNode* create_ecn();

    // packing_signatures_[i] corresponds to the logical_block_type pointed to by cluster_logical_block_types_[i].
    std::vector<const t_logical_block_type*> cluster_logical_block_types_;
    std::vector<PrimitiveSignatureNode*> packing_signatures_;

    // Current node for signature being constructed
    PrimitiveSignatureNode* at_psn_;

    std::unordered_map<AtomBlockId, int> packed_atoms_;

    // Lookup of 4-byte PstPinId from unique PstPin tuple.
    // This is to reduce the memory footprint of encoding connections in the PST.
    std::unordered_map<PstPin, PstPinId, PstPinHash> pin_mappings_;

    inline PstPinId get_pin_mapping(PstPin pin) {
        auto got = pin_mappings_.find(pin);
        if (got == pin_mappings_.end()) {
            PstPinId pin_id = PstPinId(pin_mappings_.size());
            pin_mappings_[pin] = pin_id;
            return pin_id;
        }
        return got->second;
    }

    // External IO bookkeeping
    struct EcnRecord {
        PstPinId connection;
        size_t external_sinks_count;
    };

    std::map<AtomNetId, std::vector<PstPinId>> input_nets_;
    std::map<AtomNetId, EcnRecord> output_nets_;

    // Checkpoint bookkeeping for rolling back to a legal packing signature
    PrimitiveSignatureNode* checkpoint_psn_;
    std::unordered_map<AtomNetId, size_t> checkpoint_input_nets_;
    std::unordered_map<AtomNetId, size_t> checkpoint_decremented_output_nets_;
    std::vector<AtomNetId> checkpoint_new_output_nets_;
    std::vector<AtomBlockId> checkpoint_new_atoms_;

};
