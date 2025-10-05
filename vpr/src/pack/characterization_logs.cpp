#include <format>
#include <algorithm>

#include "atom_netlist.h"
#include "physical_types.h"
#include "globals.h"
#include "cluster_legalizer.h"

#include "characterization_logs.h"

PackSignatureTree g_pack_signatures; // TODO this should not be a global in the end. Part of ClusterLegalizer?

void PackSignatureTree::start_pack_signature(const t_logical_block_type* cluster_logical_block_type) {
    for (size_t i = 0; i < cluster_logical_block_types_.size(); i++) {
        if (cluster_logical_block_types_[i] != cluster_logical_block_type) continue;

        // existing cluster type
        at_node_ = signatures_[i];
        at_node_->visits++;
        return;
    }

    // new cluster type
    at_node_ = new PackSignatureTreeNode;
    at_node_->visits = 1; // XXX characterization
    cluster_logical_block_types_.push_back(cluster_logical_block_type);
    signatures_.push_back(at_node_);
}

void PackSignatureTree::add_primitive(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id) {
    VTR_ASSERT(at_node_ != nullptr);

    PackSignaturePrimitive* primitive = this->generate_primitive(primitive_pb_graph_node, atom_block_id);

    // Determine whether path with this primitive primitive already exists.
    // Similar packing primitives are likely to appear close to eachother due to greedy candidate selection,
    // so iterate over the list in reverse to take better advantage of this locality.
    for (ssize_t i = at_node_->child_primitives.size() - 1; i >= 0 ; i--) {
        if (*at_node_->child_primitives[i] == *primitive) {
            delete primitive;
            at_node_->child_primitives[i]->atom_block_id = atom_block_id;
            at_node_ = at_node_->child_nodes[i];
            at_node_->visits++; // XXX characterization
            return;
        }
    }

    // This is a new diverging path, so add the primitive to the tree and create a new child node
    at_node_->child_primitives.push_back(primitive);

    PackSignatureTreeNode* new_node = new PackSignatureTreeNode;
    new_node->parent = at_node_;
    new_node->primitive = primitive;
    at_node_->child_nodes.push_back(new_node);

    at_node_ = new_node;
    at_node_->visits = 1; // XXX characterization
}

PackSignaturePrimitive* PackSignatureTree::generate_primitive(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id) {
    PackSignaturePrimitive* primitive = new PackSignaturePrimitive;
    primitive->pb_graph_node = primitive_pb_graph_node;
    primitive->atom_block_id = atom_block_id;

    const AtomNetlist& atom_netlist = g_vpr_ctx.atom().netlist();
    AtomNetlist::pin_range source_input_pins = atom_netlist.block_input_pins(atom_block_id);
    primitive->num_incoming_sources = source_input_pins.size();

    // Identify if any primitives that have already been placed drive this new primitive
    for (AtomPinId primitive_input_pin_id : source_input_pins) {
        AtomNetId primitive_input_pin_net_id = atom_netlist.pin_net(primitive_input_pin_id);
        AtomBlockId source_atom_block_id = atom_netlist.net_driver_block(primitive_input_pin_net_id);
        VTR_ASSERT(source_atom_block_id != AtomBlockId::INVALID());

        // Walk up the tree back to the seed primitive, checking if any already placed primitive is the source of this pin
        for (PackSignatureTreeNode* probe = at_node_; probe->parent != nullptr; probe = probe->parent) {
            if (source_atom_block_id == probe->primitive->atom_block_id) {
                AtomPinId source_pin_id = atom_netlist.net_driver(primitive_input_pin_net_id);
                AtomPortId source_port_id = atom_netlist.pin_port(source_pin_id);

                PackSignatureConnection source_connection = {
                    probe->primitive->pb_graph_node,
                    atom_netlist.port_name(source_port_id),
                    atom_netlist.pin_port_bit(source_pin_id)
                };

                primitive->intracluster_sources_to_primitive_inputs.push_back(source_connection);
            }
        }
    }

    // pin_range order could differ between equivalent clusters, so sort the sources list to ensure that primitives are comparable.
    std::sort(primitive->intracluster_sources_to_primitive_inputs.begin(), primitive->intracluster_sources_to_primitive_inputs.end());

    primitive->num_sinks = 0;
    for (AtomPinId primitive_output_pin_id : atom_netlist.block_output_pins(atom_block_id)) {
        AtomNetId primitive_output_net_id = atom_netlist.pin_net(primitive_output_pin_id);
        primitive->num_sinks += atom_netlist.net_sinks(primitive_output_net_id).size();
    }

    // Since this block could have arbitrarily many sinks, rather than check all the sinks of this block
    // to see if they are already in the cluster, the search space is reduced if we instead check to see
    // if the source of any of the blocks already packed in this cluster is this block.
    for (PackSignatureTreeNode* probe = at_node_; probe->parent != nullptr; probe = probe->parent) {
        AtomNetlist::pin_range potential_sink_input_pins = atom_netlist.block_input_pins(probe->primitive->atom_block_id);

        for (AtomPinId potential_sink_pin_id : potential_sink_input_pins) {
            AtomNetId potential_sink_net_id = atom_netlist.pin_net(potential_sink_pin_id);
            AtomBlockId source_atom_block_id = atom_netlist.net_driver_block(potential_sink_net_id);
            VTR_ASSERT(source_atom_block_id != AtomBlockId::INVALID());

            if (source_atom_block_id == atom_block_id) {
                AtomPinId primitive_output_pin_id = atom_netlist.net_driver(potential_sink_net_id);
                AtomPortId primitive_output_port_id = atom_netlist.pin_port(primitive_output_pin_id);
                PackSignatureConnection sink_connection = {
                    probe->primitive->pb_graph_node,
                    atom_netlist.port_name(primitive_output_port_id),
                    atom_netlist.pin_port_bit(primitive_output_pin_id)
                };
                primitive->intracluster_sinks_of_primitive_outputs.push_back(sink_connection);
                break;
            }
        }
    }

    return primitive;
}

// ================================================================
// START CHARACTERIZATION ONLY CODE
// ================================================================

std::ofstream g_logfile;
bool logfile_open = false;

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

size_t total_finalized_clusters = 0;

static void recurse_placement_dependent(
    const PackSignatureTreeNode* node,
    size_t depth
) {
    if (node->parent != nullptr) {
        for (size_t i = 0; i < depth; i++) g_logfile << "| ";

        std::string intracluster_sources_to_primitive_inputs_string = "{ ";
        for (size_t i = 0; i < node->primitive->intracluster_sources_to_primitive_inputs.size(); i++) {
            intracluster_sources_to_primitive_inputs_string += node->primitive->intracluster_sources_to_primitive_inputs[i].to_string();
            if (i < node->primitive->intracluster_sources_to_primitive_inputs.size() - 1) {
                intracluster_sources_to_primitive_inputs_string += ", ";
            }
        }
        intracluster_sources_to_primitive_inputs_string += " }";

        std::string intracluster_sinks_of_primitive_outputs_string = "{ ";
        for (size_t i = 0; i < node->primitive->intracluster_sinks_of_primitive_outputs.size(); i++) {
            intracluster_sinks_of_primitive_outputs_string += node->primitive->intracluster_sinks_of_primitive_outputs[i].to_string();
            if (i < node->primitive->intracluster_sinks_of_primitive_outputs.size() - 1) {
                intracluster_sinks_of_primitive_outputs_string += ", ";
            }
        }
        intracluster_sinks_of_primitive_outputs_string += " }";

        g_logfile << std::format("<\"{}\", {:#08x}>: {{ num_incoming_sources: {}, drivers: {}, num_sinks: {}, driven: {}, visits: {} }}",
                                 node->primitive->pb_graph_node->pb_type->name,
                                 reinterpret_cast<uintptr_t>(node->primitive->pb_graph_node),
                                 node->primitive->num_incoming_sources,
                                 intracluster_sources_to_primitive_inputs_string,
                                 node->primitive->num_sinks,
                                 intracluster_sinks_of_primitive_outputs_string,
                                 node->visits);
    }

    if (node->successful_legalization_cluster_ids.size() > 0) {
        g_logfile << " [" << node->successful_legalization_cluster_ids.size() << " CLUSTERS]: { ";
        for (size_t i = 0; i < node->successful_legalization_cluster_ids.size(); i++) {
            g_logfile << node->successful_legalization_cluster_ids[i] << ((node->successful_legalization_cluster_detailedness[i]) ? "!" : "") << " ";
        }
        g_logfile << "}";
        total_finalized_clusters += node->successful_legalization_cluster_ids.size();
    }

    if (node->failed_legalization_cluster_ids.size() > 0) {
        g_logfile << " [" << node->failed_legalization_cluster_ids.size() << " FAILED]: { ";
        for (size_t i = 0; i < node->failed_legalization_cluster_ids.size(); i++) {
            g_logfile << node->failed_legalization_cluster_ids[i] << ((node->failed_legalization_cluster_detailedness[i]) ? "!" : "") << " ";
        }
        g_logfile << "}";
    }

    g_logfile << std::endl;

    // End of path reached
    if (node->child_nodes.empty()) return;

    for (PackSignatureTreeNode* child_node : node->child_nodes) {
        recurse_placement_dependent(child_node, depth + 1);
    }
}

void PackSignatureTree::log_equivalent() {
    for (size_t i = 0; i < cluster_logical_block_types_.size(); i++) {
        if (cluster_logical_block_types_[i]->name != "clb") continue;
        g_logfile << std::format("cluster_pb_type: <{:#08x}, {}>", reinterpret_cast<uintptr_t>(cluster_logical_block_types_[i]), cluster_logical_block_types_[i]->name);
        recurse_placement_dependent(signatures_[i], 0);
        g_logfile << "TOTAL FINALIZED CLUSTERS: " << total_finalized_clusters << std::endl << std::endl;
    }
    g_logfile << "SPECULATIVE LEGALIZATION SUCCESS TIME: " << speculative_legalization_success_duration << std::endl;
    g_logfile << "SPECULATIVE LEGALIZATION FAILURE TIME: " << speculative_legalization_failure_duration << std::endl;
    g_logfile << "DETAILED LEGALIZATION SUCCESS TIME: " << detailed_legalization_success_duration << std::endl;
    g_logfile << "DETAILED LEGALIZATION FAILURE TIME: " << detailed_legalization_failure_duration << std::endl;
}
