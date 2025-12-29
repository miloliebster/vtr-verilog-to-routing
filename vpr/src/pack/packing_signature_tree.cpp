#include <format>
#include <algorithm>

#include "atom_netlist.h"
#include "physical_types.h"
#include "globals.h"
#include "cluster_legalizer.h"

#include <chrono>
#include "packing_signature_tree.h"

void PackingSignatureTree::start_packing_signature(const t_logical_block_type* cluster_logical_block_type) {
    // Reset external IO bookkeeping
    input_nets_.clear();
    output_nets_.clear();
    packed_atoms_.clear();
    checkpoint_input_nets_.clear();
    checkpoint_new_output_nets_.clear();
    checkpoint_decremented_output_nets_.clear();
    checkpoint_new_atoms_.clear();
    legalization_checks = 0;

    routed = false;

    for (size_t i = 0; i < cluster_logical_block_types_.size(); i++) {
        if (cluster_logical_block_types_[i] != cluster_logical_block_type) continue;
        // Existing cluster type
        at_psn_ = packing_signatures_[i];
        return;
    }

    // New cluster type
    at_psn_ = new PrimitiveSignatureNode;
    checkpoint_psn_ = at_psn_;
    cluster_logical_block_types_.push_back(cluster_logical_block_type);
    packing_signatures_.push_back(at_psn_);
}

void PackingSignatureTree::add_psn(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id) {
    auto start_time = std::chrono::high_resolution_clock::now(); // XXX
    VTR_ASSERT(at_psn_ != nullptr);

    PrimitiveSignatureNode* new_psn = this->create_psn(primitive_pb_graph_node, atom_block_id);

    // Determine whether path with this primitive already exists.
    // Similar packing primitives are likely to appear close to each other due to greedy candidate selection,
    // so iterate over the list in reverse to take better advantage of this locality.
    for (ssize_t i = at_psn_->child_psn.size() - 1; i >= 0; i--) {
        PrimitiveSignatureNode* child_psn = at_psn_->child_psn[i];
        if (*child_psn == *new_psn) {
            delete new_psn;
            at_psn_ = child_psn;
            auto end_time = std::chrono::high_resolution_clock::now(); // XXX
            signature_processing_duration += end_time - start_time; // XXX
            return;
        }
    }

    // This is a new diverging path, so add the primitive to the tree and create a new child node
    at_psn_->child_psn.push_back(new_psn);
    at_psn_ = new_psn;

    auto end_time = std::chrono::high_resolution_clock::now(); // XXX
    signature_processing_duration += end_time - start_time;
}

PrimitiveSignatureNode* PackingSignatureTree::create_psn(const t_pb_graph_node* primitive_pb_graph_node, const AtomBlockId atom_block_id) {
    PrimitiveSignatureNode* psn = new PrimitiveSignatureNode;
    psn->primitive_num = primitive_pb_graph_node->primitive_num;
    packed_atoms_[atom_block_id] = psn->primitive_num;
    checkpoint_new_atoms_.push_back(atom_block_id);

    const AtomNetlist& atom_netlist = g_vpr_ctx.atom().netlist();
    AtomNetlist::pin_range primitive_input_pins = atom_netlist.block_input_pins(atom_block_id);

    // Determine how many sinks each net that the output pins for this block drive.
    // If we find that these sinks are inside the cluster later, we decrement the external sinks count for the net.
    // Once packing is complete, if the value for the net is >0, then we know there are sinks outside of the cluster
    for (AtomPinId primitive_output_pin_id : atom_netlist.block_output_pins(atom_block_id)) {
        AtomPortId primitive_output_port_id = atom_netlist.pin_port(primitive_output_pin_id);
        AtomNetId primitive_output_net_id = atom_netlist.pin_net(primitive_output_pin_id);

        PstPin primitive_output_psn_pin = std::make_tuple(primitive_pb_graph_node->primitive_num,
                                                          atom_netlist.port_name(primitive_output_port_id),
                                                          atom_netlist.pin_port_bit(primitive_output_pin_id));

        output_nets_[primitive_output_net_id] = EcnRecord(
            this->get_pin_mapping(primitive_output_psn_pin),
            atom_netlist.net_sinks(primitive_output_net_id).size()
        );

        checkpoint_new_output_nets_.push_back(primitive_output_net_id);
    }

    for (AtomPinId primitive_input_pin_id : primitive_input_pins) {
        AtomNetId primitive_input_pin_net_id = atom_netlist.pin_net(primitive_input_pin_id);
        AtomPortId primitive_input_pin_port_id = atom_netlist.pin_port(primitive_input_pin_id);
        AtomBlockId source_atom_block_id = atom_netlist.net_driver_block(primitive_input_pin_net_id);
        VTR_ASSERT(source_atom_block_id != AtomBlockId::INVALID());

        // Create record of groupings of primitive pins driven by the same net.
        // This is used by the final path node to identify pins which are driven by an external net.
        PstPin input_psn_pin = std::make_tuple(primitive_pb_graph_node->primitive_num,
                                               atom_netlist.port_name(primitive_input_pin_port_id),
                                               atom_netlist.pin_port_bit(primitive_input_pin_id));

        input_nets_[primitive_input_pin_net_id].push_back(this->get_pin_mapping(input_psn_pin));
        checkpoint_input_nets_[primitive_input_pin_net_id]++;

        // Identify if any primitives that have already been placed drive this new primitive
        // Walk up the tree back to the seed primitive, checking if any already placed primitive is the source of this pin
        auto got = packed_atoms_.find(source_atom_block_id);
        if (got != packed_atoms_.end()) {
            AtomPinId source_pin_id = atom_netlist.net_driver(primitive_input_pin_net_id);
            AtomPortId source_port_id = atom_netlist.pin_port(source_pin_id);
            int source_primitive_num = got->second;

            PstPin source_psn_pin = std::make_tuple(source_primitive_num,
                                                    atom_netlist.port_name(source_port_id),
                                                    atom_netlist.pin_port_bit(source_pin_id));
            PstPin sink_psn_pin = std::make_tuple(psn->primitive_num,
                                                  atom_netlist.port_name(primitive_input_pin_port_id),
                                                  atom_netlist.pin_port_bit(primitive_input_pin_id));
            PstConnection source_connection{ this->get_pin_mapping(source_psn_pin), this->get_pin_mapping(sink_psn_pin) };
            psn->intracluster_sources_to_primitive_inputs.push_back(source_connection);

            EcnRecord& record = output_nets_[primitive_input_pin_net_id];
            if (record.external_sinks_count <= 0) {
                log_equivalent();
            }
            VTR_ASSERT(record.external_sinks_count > 0);
            record.external_sinks_count--;
            checkpoint_decremented_output_nets_[primitive_input_pin_net_id]++;
        }
    }

    // pin_range order could differ between equivalent clusters, so sort the sources list to ensure that primitives are comparable.
    std::sort(psn->intracluster_sources_to_primitive_inputs.begin(), psn->intracluster_sources_to_primitive_inputs.end());

    // Since this block could have arbitrarily many sinks, rather than check all the sinks of this block
    // to see if they are already in the cluster, the search space is reduced if we instead check to see
    // if the source of any of the blocks already packed in this cluster is this block.
    for (auto maybe_sink_atom : packed_atoms_) {
        AtomBlockId maybe_sink_atom_block_id = maybe_sink_atom.first;
        AtomNetlist::pin_range potential_sink_input_pins = atom_netlist.block_input_pins(maybe_sink_atom_block_id);
        int maybe_sink_primitive_num = maybe_sink_atom.second;

        for (AtomPinId potential_sink_pin_id : potential_sink_input_pins) {
            AtomNetId potential_sink_net_id = atom_netlist.pin_net(potential_sink_pin_id);
            AtomBlockId source_atom_block_id = atom_netlist.net_driver_block(potential_sink_net_id);
            VTR_ASSERT(source_atom_block_id != AtomBlockId::INVALID());

            if (source_atom_block_id == atom_block_id) {
                AtomPinId primitive_output_pin_id = atom_netlist.net_driver(potential_sink_net_id);
                AtomPortId primitive_output_port_id = atom_netlist.pin_port(primitive_output_pin_id);
                AtomPinId sink_pin_id = potential_sink_pin_id;
                AtomPortId sink_port_id = atom_netlist.pin_port(sink_pin_id);

                PstPin source_psn_pin = std::make_tuple(psn->primitive_num,
                                                        atom_netlist.port_name(primitive_output_port_id),
                                                        atom_netlist.pin_port_bit(primitive_output_pin_id));
                PstPin sink_psn_pin = std::make_tuple(maybe_sink_primitive_num,
                                                      atom_netlist.port_name(sink_port_id),
                                                      atom_netlist.pin_port_bit(sink_pin_id));
                PstConnection sink_connection{ this->get_pin_mapping(source_psn_pin), this->get_pin_mapping(sink_psn_pin) };
                psn->intracluster_sinks_of_primitive_outputs.push_back(sink_connection);

                EcnRecord& record = output_nets_[potential_sink_net_id];
                VTR_ASSERT(record.external_sinks_count > 0);
                record.external_sinks_count--;

                break;
            }
        }
    }

    return psn;
}

void PackingSignatureTree::set_checkpoint() {
    auto start_time = std::chrono::high_resolution_clock::now(); // XXX
    checkpoint_psn_ = at_psn_;
    checkpoint_input_nets_.clear();
    checkpoint_new_output_nets_.clear();
    checkpoint_decremented_output_nets_.clear();
    checkpoint_new_atoms_.clear();
    auto end_time = std::chrono::high_resolution_clock::now(); // XXX
    signature_processing_duration += end_time - start_time; // XXX
}

void PackingSignatureTree::rollback_to_checkpoint() {
    auto start_time = std::chrono::high_resolution_clock::now(); // XXX
    at_psn_ = checkpoint_psn_;
    VTR_ASSERT(at_psn_ != nullptr);
    for (auto it = checkpoint_input_nets_.begin(); it != checkpoint_input_nets_.end(); it++) {
        if (input_nets_[it->first].size() == it->second) {
            input_nets_.erase(it->first);
        } else {
            input_nets_[it->first].resize(input_nets_[it->first].size() - it->second);
        }
    }
    for (auto it = checkpoint_decremented_output_nets_.begin(); it != checkpoint_decremented_output_nets_.end(); it++) {
        output_nets_[it->first].external_sinks_count += it->second;
    }
    for (auto net : checkpoint_new_output_nets_) {
        output_nets_.erase(net);
    }
    for  (auto atom : checkpoint_new_atoms_) {
        packed_atoms_.erase(atom);
    }
    auto end_time = std::chrono::high_resolution_clock::now(); // XXX
    signature_processing_duration += end_time - start_time; // XXX
}

ExternalConnectivityNode* PackingSignatureTree::create_ecn() {
    ExternalConnectivityNode* new_ecn = new ExternalConnectivityNode;

    for (auto input_net : input_nets_) {
        if (output_nets_.count(input_net.first) != 0) continue; // net is driven from inside cluster; not an external source
        std::sort(input_net.second.begin(), input_net.second.end());
        new_ecn->cluster_inputs.push_back(input_net.second);
    }
    std::sort(new_ecn->cluster_inputs.begin(), new_ecn->cluster_inputs.end(), [](auto a, auto b) { return a < b; });

    for (auto output_net : output_nets_) {
        if (output_net.second.external_sinks_count == 0) continue; // net only drives pins inside cluster
        new_ecn->cluster_outputs.push_back(output_net.second.connection);
    }
    std::sort(new_ecn->cluster_outputs.begin(), new_ecn->cluster_outputs.end());

    return new_ecn;
}

e_packing_signature_legality PackingSignatureTree::check_legality() {
    auto start_time = std::chrono::high_resolution_clock::now(); // XXX
    legalization_checks++;

    ExternalConnectivityNode* ecn = this->create_ecn();
    for (ssize_t i = at_psn_->child_ecn.size() - 1; i >= 0 ; i--) {
        ExternalConnectivityNode* child_ecn = at_psn_->child_ecn[i];
        if (*child_ecn == *ecn) {
            delete ecn;
            return (at_psn_->child_ecn[i]->legal) ? e_packing_signature_legality::LEGAL : e_packing_signature_legality::ILLEGAL;
        }
    }
    delete ecn;
    return e_packing_signature_legality::UNKNOWN;

    auto end_time = std::chrono::high_resolution_clock::now(); // XXX
    signature_processing_duration += end_time - start_time; // XXX
}

void PackingSignatureTree::mark_signature_as_legal(LegalizationClusterId legalization_cluster_id) {
    auto start_time = std::chrono::high_resolution_clock::now(); // XXX
                                                                 //
    ExternalConnectivityNode* ecn = this->create_ecn();
    ecn->legal = true;
    ecn->detailed_legalization_checks = legalization_checks; // XXX

    for (ssize_t i = at_psn_->child_ecn.size() - 1; i >= 0 ; i--) {
        ExternalConnectivityNode* child_ecn = at_psn_->child_ecn[i];
        if (*child_ecn == *ecn) {
            if (legalization_cluster_id != LegalizationClusterId::INVALID()) { // XXX
                child_ecn->successful_clusters++; // XXX
                if (this->detailed_legalization) child_ecn->successful_detailed_clusters++; // XXX
            } // XXX
            delete ecn;
            return;
        }
    }
    if (legalization_cluster_id != LegalizationClusterId::INVALID()) { // XXX
        ecn->successful_clusters++; // XXX
        if (this->detailed_legalization) ecn->successful_detailed_clusters++; // XXX
    } // XXX
    at_psn_->child_ecn.push_back(ecn);

    auto end_time = std::chrono::high_resolution_clock::now(); // XXX
    signature_processing_duration += end_time - start_time; // XXX
}

void PackingSignatureTree::mark_signature_as_illegal(LegalizationClusterId legalization_cluster_id) {
    auto start_time = std::chrono::high_resolution_clock::now(); // XXX
                                                                 //
    ExternalConnectivityNode* ecn = this->create_ecn();
    ecn->legal = false;
    ecn->detailed_legalization_checks = legalization_checks; // XXX

    for (ssize_t i = at_psn_->child_ecn.size() - 1; i >= 0 ; i--) {
        ExternalConnectivityNode* child_ecn = at_psn_->child_ecn[i];
        if (*child_ecn == *ecn) {
            if (legalization_cluster_id != LegalizationClusterId::INVALID()) { // XXX
                child_ecn->failed_clusters++; // XXX
            } // XXX
            delete ecn;
            return;
        }
    }
    if (legalization_cluster_id != LegalizationClusterId::INVALID()) { // XXX
        ecn->failed_clusters++; // XXX
    } // XXX
    at_psn_->child_ecn.push_back(ecn);

    auto end_time = std::chrono::high_resolution_clock::now(); // XXX
    signature_processing_duration += end_time - start_time; // XXX
}

void PackingSignatureTree::finalize_path(LegalizationClusterId legalization_cluster_id) {
    this->mark_signature_as_legal(legalization_cluster_id);
}

void PackingSignatureTree::fail_path(LegalizationClusterId legalization_cluster_id) {
    this->mark_signature_as_illegal(legalization_cluster_id);
}

// ================================================================
// START CHARACTERIZATION ONLY CODE
// ================================================================

size_t total_finalized_clusters = 0;

static void recurse_placement_dependent(const PrimitiveSignatureNode* psn, size_t depth, std::ofstream& logfile) {
    if (psn->primitive_num != -1) {
        for (size_t i = 0; i < depth; i++) logfile << "| ";

        std::string intracluster_sources_to_primitive_inputs_string = "{ ";
        for (size_t i = 0; i < psn->intracluster_sources_to_primitive_inputs.size(); i++) {
            intracluster_sources_to_primitive_inputs_string += psn->intracluster_sources_to_primitive_inputs[i].to_string();
            if (i < psn->intracluster_sources_to_primitive_inputs.size() - 1) {
                intracluster_sources_to_primitive_inputs_string += ", ";
            }
        }
        intracluster_sources_to_primitive_inputs_string += " }";

        std::string intracluster_sinks_of_primitive_outputs_string = "{ ";
        for (size_t i = 0; i < psn->intracluster_sinks_of_primitive_outputs.size(); i++) {
            intracluster_sinks_of_primitive_outputs_string += psn->intracluster_sinks_of_primitive_outputs[i].to_string();
            if (i < psn->intracluster_sinks_of_primitive_outputs.size() - 1) {
                intracluster_sinks_of_primitive_outputs_string += ", ";
            }
        }
        intracluster_sinks_of_primitive_outputs_string += " }";

        logfile << std::format("loc{}: {{ drivers: {}, driven: {} }}",
                                 psn->primitive_num,
                                 intracluster_sources_to_primitive_inputs_string,
                                 intracluster_sinks_of_primitive_outputs_string);
        logfile << std::endl;
    }

    for (auto ecn : psn->child_ecn) {
        for (size_t i = 0; i < depth + 1; i++) logfile << "| ";

        std::string cluster_inputs_string = "{ ";
        for (size_t i = 0; i < ecn->cluster_inputs.size(); i++) {
            cluster_inputs_string += "{ ";
            for (size_t j = 0; j < ecn->cluster_inputs[i].size(); j++) {
                cluster_inputs_string += "pin";
                cluster_inputs_string += std::to_string((int)ecn->cluster_inputs[i][j]);
                if (j < ecn->cluster_inputs[i].size() - 1) {
                    cluster_inputs_string += ", ";
                }
            }
            cluster_inputs_string += " }";
            if (i < ecn->cluster_inputs.size() - 1) {
                cluster_inputs_string += ", ";
            }
        }
        cluster_inputs_string += " }";

        std::string cluster_outputs_string = "{ ";
        for (size_t i = 0; i < ecn->cluster_outputs.size(); i++) {
            cluster_outputs_string += "pin";
            cluster_outputs_string += std::to_string((int)ecn->cluster_outputs[i]);
            if (i < ecn->cluster_outputs.size() - 1) {
                cluster_outputs_string += ", ";
            }
        }
        cluster_outputs_string += " }";

        logfile << std::format("IO: {{ legal: {}, cluster_inputs: {}, cluster_outputs: {} }}",
                                 ecn->legal,
                                 cluster_inputs_string,
                                 cluster_outputs_string);

        if (ecn->successful_clusters > 0) {
            logfile << " [" << ecn->successful_clusters << " CLUSTERS ; " << ecn->successful_detailed_clusters << " DETAILED" << " ~ " << ecn->detailed_legalization_checks;
        }

        if (ecn->failed_clusters > 0) {
            logfile << " [" << ecn->failed_clusters << " FAILED";
        }

        logfile << "]";
        total_finalized_clusters += ecn->successful_clusters;
        logfile << std::endl;

        VTR_ASSERT(
            (ecn->failed_clusters    == 0  && ecn->successful_clusters == 0) ||
            (ecn->failed_clusters     > 0  && ecn->successful_clusters == 0) ||
            (ecn->successful_clusters > 0  && ecn->failed_clusters     == 0)
        );
    }

    for (PrimitiveSignatureNode* child_psn : psn->child_psn) {
        recurse_placement_dependent(child_psn, depth + 1, logfile);
    }
}

void PackingSignatureTree::log_equivalent() {
    std::ofstream logfile;
    logfile.open("pst.txt");
    for (size_t i = 0; i < cluster_logical_block_types_.size(); i++) {
        logfile << std::format("cluster_pb_type: <{:#08x}, {}>", reinterpret_cast<uintptr_t>(cluster_logical_block_types_[i]), cluster_logical_block_types_[i]->name) << std::endl;
        recurse_placement_dependent(packing_signatures_[i], 0, logfile);
        logfile << "TOTAL FINALIZED CLUSTERS: " << total_finalized_clusters << std::endl << std::endl;
    }
    logfile << "SPECULATIVE LEGALIZATION SUCCESS TIME: " << speculative_legalization_success_duration << std::endl;
    logfile << "SPECULATIVE LEGALIZATION FAILURE TIME: " << speculative_legalization_failure_duration << std::endl;
    logfile << "DETAILED LEGALIZATION SUCCESS TIME: " << detailed_legalization_success_duration << std::endl;
    logfile << "DETAILED LEGALIZATION FAILURE TIME: " << detailed_legalization_failure_duration << std::endl << std::endl;

    logfile << "SIGNATURE PROCESSING TIME: " << signature_processing_duration << std::endl << std::endl;

    for (auto pin : pin_mappings_) {
        logfile << "PIN " << pin.second << " (" << std::get<0>(pin.first) << "," << std::get<1>(pin.first) << "," << (int)std::get<2>(pin.first) << ")" << std::endl;
    }
    logfile.close();
}
