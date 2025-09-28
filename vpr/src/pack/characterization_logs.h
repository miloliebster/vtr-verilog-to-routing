#pragma once

#include "atom_netlist.h"
#include "physical_types.h"
#include "globals.h"
#include "cluster_legalizer.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <string>
#include <vector>
#include <utility>

void try_open_logfile();

extern std::ofstream g_logfile;

struct PrimitiveSignature {
    const t_pb_graph_node* primitive_location;
    size_t num_drivers;
    size_t num_driven;
    std::vector<const t_pb_graph_node*> intracluster_drivers;
    std::vector<const t_pb_graph_node*> intracluster_driven;
    AtomBlockId blk_id; // not part of signature per se, but referenced when backtracking

    bool operator==(PrimitiveSignature const& rhs) const {
        if (this->primitive_location != rhs.primitive_location) return false;
        if (this->num_drivers != rhs.num_drivers) return false;
        if (this->num_driven != rhs.num_driven) return false;
        if (this->intracluster_drivers.size() != rhs.intracluster_drivers.size()) return false;
        if (this->intracluster_driven.size() != rhs.intracluster_driven.size()) return false;
        for (size_t i = 0; i < this->intracluster_drivers.size(); i++) {
            if (this->intracluster_drivers[i] != rhs.intracluster_drivers[i]) return false;
        }
        for (size_t i = 0; i < this->intracluster_driven.size(); i++) {
            if (this->intracluster_driven[i] != rhs.intracluster_driven[i]) return false;
        }
        return true;
    }
};

struct PackSignatureNode {
    size_t visits = 0;
    std::vector<LegalizationClusterId> legalization_cluster_ids;
    PrimitiveSignature* signature = nullptr;
    PackSignatureNode* parent = nullptr;

    std::vector<PrimitiveSignature*> child_signatures;
    std::vector<PackSignatureNode*> child_nodes;
};

class PackSignatures {
public:
    void start_signature(const t_logical_block_type* cluster_type);
    void add_primitive(const t_pb_graph_node* primitive_location, const AtomBlockId blk_id);

    // // TODO: can maybe remove redundant checks and searched once I understand what is happening
    // inline void walk_back(AtomBlockId blk_id) {
    //     PackSignatureNode* probe = at_node_;
    //
    //     while (probe->parent != nullptr && probe->signature->blk_id != blk_id) probe = probe->parent;
    //     if (probe->parent == nullptr) return;
    //     VTR_ASSERT(probe->signature->blk_id == blk_id);
    //     VTR_ASSERT(probe == at_node_);
    //     VTR_ASSERT(at_node_->visits == 1);
    //
    //     at_node_ = at_node_->parent;
    //     VTR_ASSERT(at_node_->child_nodes.back() == probe);
    //     VTR_ASSERT(at_node_->child_signatures.back() == probe->signature);
    //
    //     delete at_node_->child_signatures.back();
    //     at_node_->child_signatures.pop_back();
    //     delete at_node_->child_nodes.back();
    //     at_node_->child_nodes.pop_back();
    // }

    inline void finalize_path(LegalizationClusterId legalization_cluster_id) {
        at_node_->legalization_cluster_ids.push_back(legalization_cluster_id);
        finalization_calls_++;
    }

    void log_equivalent_placement_dependent();

private:
    size_t finalization_calls_ = 0;
    PrimitiveSignature* generate_primitive_signature(const t_pb_graph_node* primitive_location, const AtomBlockId atom_id);

    // cluster_signatures_[i] corresponds to the logical_block_type pointed to by cluster_types_[i].
    std::vector<const t_logical_block_type*> cluster_types_;
    std::vector<PackSignatureNode*> cluster_signatures_;

    // Current node for signature being constructed
    PackSignatureNode* at_node_;

    // pb_graph_nodes are enumerated to provide an order-dependent, but placement-independent, representation for signatures.
    // The index that a pb_graph_node pointer is stored is its ID.
    std::vector<t_pb_graph_node*> pb_graph_nodes_ordered_;

    // Maps an atom to the pb_graph_node where it was placed, since this relation is not preserved by the packing code.
    std::vector<std::pair<AtomBlockId, PrimitiveSignature*>> atom_to_primitive_signature_;
};

extern PackSignatures g_pack_signatures;
