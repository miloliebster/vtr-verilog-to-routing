#include <format>
#include <algorithm>

#include "atom_netlist.h"
#include "physical_types.h"
#include "globals.h"
#include "cluster_legalizer.h"

#include "characterization_logs.h"

std::ofstream g_logfile;
bool logfile_open = false;

PackSignatures g_pack_signatures;

void try_open_logfile() {
    if (logfile_open) return;

    char date_string[64] = {};
    time_t t = time(NULL);
    strftime(date_string, sizeof(date_string), "%F_%H-%M-%S", localtime(&t));
    std::string logfile_path = "char_";
    logfile_path += date_string;
    logfile_path += ".txt";
    g_logfile.open(logfile_path);
    std::atexit([](){ g_logfile.close(); });
    logfile_open = true;
}

size_t clb_starts = 0;

void PackSignatures::start_signature(const t_logical_block_type* cluster_type) {
    // g_logfile << "START SIGNATURE" << std::endl;
    pb_graph_nodes_ordered_.clear();
    atom_to_primitive_signature_.clear();

    clb_starts++;

    auto it = std::find(cluster_types_.begin(), cluster_types_.end(), cluster_type);
    if (it != cluster_types_.end()) {
        // existing cluster type
        at_node_ = cluster_signatures_[std::distance(cluster_types_.begin(), it)];
        at_node_->visits++;
        return;
    }
    // new cluster type
    at_node_ = new PackSignatureNode;
    at_node_->visits = 1;
    cluster_types_.push_back(cluster_type);
    cluster_signatures_.push_back(at_node_);
}

void PackSignatures::add_primitive(const t_pb_graph_node* primitive_location, const AtomBlockId blk_id) {
    VTR_ASSERT(at_node_ != nullptr);

    PrimitiveSignature* signature = this->generate_primitive_signature(primitive_location, blk_id);

    // Determine whether path with this primitive signature already exists
    for(size_t i = 0; i < at_node_->child_signatures.size(); i++) {
        if (*at_node_->child_signatures[i] == *signature) {
            VTR_ASSERT(at_node_->child_nodes.size() > i);
            delete signature;
            at_node_->child_signatures[i]->blk_id = blk_id;
            at_node_ = at_node_->child_nodes[i];
            at_node_->visits++;
            return;
        }
    }

    // This is a new path, so add signature and create a new child node
    at_node_->child_signatures.push_back(signature);
    atom_to_primitive_signature_.push_back(std::make_pair(blk_id, signature));

    PackSignatureNode* new_node = new PackSignatureNode;
    new_node->parent = at_node_;
    new_node->signature = signature;
    at_node_->child_nodes.push_back(new_node);

    at_node_ = new_node;
    at_node_->visits = 1;
}

PrimitiveSignature* PackSignatures::generate_primitive_signature(const t_pb_graph_node* primitive_location, const AtomBlockId blk_id) {
    PrimitiveSignature* signature = new PrimitiveSignature;
    signature->primitive_location = primitive_location;
    signature->blk_id = blk_id;

    const AtomNetlist& atom_netlist = g_vpr_ctx.atom().netlist();

    AtomNetlist::pin_range input_pins = atom_netlist.block_input_pins(blk_id);
    signature->num_drivers = input_pins.size();

    // Identify any source primitives which have already been mapped to this cluster
    for (AtomPinId pin_id : input_pins) {
        AtomNetId net_id = atom_netlist.pin_net(pin_id);
        AtomBlockId driver_blk_id = atom_netlist.net_driver_block(net_id);
        VTR_ASSERT(driver_blk_id != AtomBlockId::INVALID());

        // Follow tree back to root, looking for if any placed primitives drive this primitive
        for (PackSignatureNode* probe = at_node_; probe->parent != nullptr; probe = probe->parent) {
            if (driver_blk_id == probe->signature->blk_id) {
                signature->intracluster_drivers.push_back(probe->signature->primitive_location);
            }
        }
    }

    // pin_range order could differ between equivalent clusters, so sort the input primitive pointers
    // to ensure that signatures are comparable.
    std::sort(signature->intracluster_drivers.begin(), signature->intracluster_drivers.end());

    // TODO: signature needs to encode distinctions between different output pins (good enough for now for characterization)
    signature->num_driven = 0;
    for (AtomPinId pin_id : atom_netlist.block_output_pins(blk_id)) {
        AtomNetId net_id = atom_netlist.pin_net(pin_id);
        signature->num_driven += atom_netlist.net_sinks(net_id).size();
    }

    // Since this block could have arbitrarily many sinks, rather than check all the sinks of this block
    // to see if they are already in the cluster, the search space is reduced if we instead check to see
    // if the source of any of the blocks already packed in this cluster is this block.
    for (PackSignatureNode* probe = at_node_; probe->parent != nullptr; probe = probe->parent) {
        input_pins = atom_netlist.block_input_pins(probe->signature->blk_id);
        for (AtomPinId pin_id : input_pins) {
            AtomNetId net_id = atom_netlist.pin_net(pin_id);
            AtomBlockId driver_blk_id = atom_netlist.net_driver_block(net_id);
            VTR_ASSERT(driver_blk_id != AtomBlockId::INVALID());

            if (driver_blk_id == blk_id) {
                const t_pb_graph_node* driven_primitive_location = probe->signature->primitive_location;
                signature->intracluster_driven.push_back(driven_primitive_location);
                break;
            }
        }
    }

    return signature;
}

size_t total_finalized_clusters = 0;

static void recurse_placement_dependent(
    const PackSignatureNode* node,
    size_t depth
) {
    if (node->parent != nullptr) {
        for (size_t i = 0; i < depth; i++) g_logfile << "| ";

        std::string intracluster_drivers_string = "{ ";
        for (size_t i = 0; i < node->signature->intracluster_drivers.size(); i++) {
            intracluster_drivers_string += std::format("<{:#08x}, {}>",
                                                       reinterpret_cast<uintptr_t>(node->signature->intracluster_drivers[i]),
                                                       node->signature->intracluster_drivers[i]->pb_type->name);
            if (i < node->signature->intracluster_drivers.size() - 1) {
                intracluster_drivers_string += ", ";
            }
        }
        intracluster_drivers_string += " }";

        std::string intracluster_driven_string = "{ ";
        for (size_t i = 0; i < node->signature->intracluster_driven.size(); i++) {
            intracluster_driven_string += std::format("<{:#08x}, {}>",
                                                       reinterpret_cast<uintptr_t>(node->signature->intracluster_driven[i]),
                                                       node->signature->intracluster_driven[i]->pb_type->name);
            if (i < node->signature->intracluster_driven.size() - 1) {
                intracluster_driven_string += ", ";
            }
        }
        intracluster_driven_string += " }";

        g_logfile << std::format("primitive_location: <{:#08x}, {}>, num_drivers: {}, drivers: {}, num_driven: {}, driven: {}, visits: {}",
                                 reinterpret_cast<uintptr_t>(node->signature->primitive_location),
                                 node->signature->primitive_location->pb_type->name,
                                 node->signature->num_drivers,
                                 intracluster_drivers_string,
                                 node->signature->num_driven,
                                 intracluster_driven_string,
                                 node->visits);
    }

    if (node->legalization_cluster_ids.size() > 0) {
        g_logfile << ", [" << node->legalization_cluster_ids.size() << " CLUSTERS]: { ";
        for (LegalizationClusterId id : node->legalization_cluster_ids) {
            g_logfile << id << " ";
        }
        g_logfile << "}";
        total_finalized_clusters += node->legalization_cluster_ids.size();
    }

    g_logfile << std::endl;

    // End of path reached
    if (node->child_nodes.empty()) return;

    for (PackSignatureNode* child_node : node->child_nodes) {
        // if (child_node->invalid) continue;
        recurse_placement_dependent(child_node, depth + 1);
    }
}

// static void recurse_placement_agnostic(const PackSignatureNode* node, std::string& packing_signature_string, size_t& visits) {}

void PackSignatures::log_equivalent_placement_dependent() {
    g_logfile << "CLB_STARTS: " << clb_starts << std::endl << std::endl;
    for (size_t i = 0; i < cluster_types_.size(); i++) {
        if (cluster_types_[i]->name != "clb") continue;
        g_logfile << std::format("cluster_pb_type: <{:#08x}, {}>", reinterpret_cast<uintptr_t>(cluster_types_[i]), cluster_types_[i]->name);
        // drop_invalid(cluster_signatures_[i]);
        recurse_placement_dependent(cluster_signatures_[i], 0);
        g_logfile << "TOTAL FINALIZED CLUSTERS: " << total_finalized_clusters << std::endl << std::endl;
    }
    g_logfile << "FINALIZATION CALLS: " << finalization_calls_ << std::endl;
}













