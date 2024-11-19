
#include "ap_clusterer.h"
#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>
#include "ShowSetup.h"
#include "ap_netlist.h"
#include "atom_netlist.h"
#include "check_netlist.h"
#include "cluster.h"
#include "cluster_legalizer.h"
#include "cluster_util.h"
#include "clustered_netlist.h"
#include "globals.h"
#include "pack.h"
#include "pack_types.h"
#include "partial_placement.h"
#include "place_constraints.h"
#include "prepack.h"
#include "vpr_api.h"
#include "vpr_types.h"
#include "vtr_assert.h"

APGainCluster::APGainCluster(const t_pack_molecule* mol,
                             const AtomNetlist& atom_netlist) : molecules({mol}) {
    // Compute the internal and external pins.

    // Internal pins are pins connected to a block in the molecule.
    for (AtomBlockId blk_id : mol->atom_block_ids) {
        if (!blk_id.is_valid())
            continue;
        internal_pins.insert(atom_netlist.block_pins(blk_id).begin(), atom_netlist.block_pins(blk_id).end());
    }

    // External pins are pins connected to blocks connected to a block in the
    // molecule, but are not in the molecule itself.
    // Get all of the nets connected to this molecule from the internal pins.
    std::unordered_set<AtomNetId> connected_nets;
    for (AtomPinId internal_pin : internal_pins) {
        connected_nets.insert(atom_netlist.pin_net(internal_pin));
    }
    // Get all of the pins in all nets which are not internal pins.
    for (AtomNetId connected_net : connected_nets) {
        for (AtomPinId net_pin : atom_netlist.net_pins(connected_net)) {
            // If this is an internal pin, it cannot be an external pin.
            if (internal_pins.count(net_pin) != 0)
                continue;
            external_pins.insert(net_pin);
        }
    }

    // Collect the nets
    for (AtomPinId internal_pin : internal_pins) {
        AtomNetId net = atom_netlist.pin_net(internal_pin);
        bool all_pins_internal = true;
        for (AtomPinId net_pin : atom_netlist.net_pins(net)) {
            if (internal_pins.count(net_pin) == 0) {
                all_pins_internal = false;
                break;
            }
        }
        if (all_pins_internal) {
            internal_nets.insert(net);
        } else {
            external_nets.insert(net);
        }
    }
}

APClusterGainCalculator::APClusterGainCalculator(const AtomNetlist& atom_netlist,
                                                 const Prepacker& prepacker,
                                                 const PartialPlacement& p_placement) : atom_netlist_(atom_netlist),
                                                                                        prepacker_(prepacker),
                                                                                        p_placement_(p_placement) {
    // Create a gain cluster for each molecule now. This allows us to precompute
    // data for each molecule at the start to make future computations easier.
    std::vector<t_pack_molecule*> molecules = prepacker.get_molecules_vector();
    for (const t_pack_molecule* mol : molecules) {
        APGainClusterId new_gain_cluster_id = APGainClusterId(gain_clusters_.size());
        molecule_gain_cluster_[mol] = new_gain_cluster_id;
        gain_clusters_.emplace_back(mol, atom_netlist);
        VTR_ASSERT(gain_clusters_[new_gain_cluster_id].molecules.size() == 1);
        VTR_ASSERT(gain_clusters_[new_gain_cluster_id].molecules[0] == mol);
    }
}

APGainClusterId APClusterGainCalculator::create_gain_cluster(const t_pack_molecule* mol) {
    VTR_ASSERT(mol != nullptr);
    VTR_ASSERT(molecule_gain_cluster_.count(mol) != 0);

    // Get the precomputed molecule gain cluster and ensure it has not been
    // clustered somewhere else.
    APGainClusterId molecule_gain_cluster_id = molecule_gain_cluster_[mol];
    VTR_ASSERT(molecule_gain_cluster_id.is_valid());
    APGainCluster& molecule_gain_cluster = gain_clusters_[molecule_gain_cluster_id];
    // FIXME: Change this to "is_molecule_clustered" maybe...
    VTR_ASSERT(molecule_gain_cluster.valid);

    VTR_ASSERT(molecule_gain_cluster.molecules.size() == 1);
    VTR_ASSERT(molecule_gain_cluster.molecules[0] == mol);

    // Create a new gain cluster which is identical to the molecule's gain
    // cluster.
    APGainClusterId new_gain_cluster_id = APGainClusterId(gain_clusters_.size());
    gain_clusters_.push_back(molecule_gain_cluster);

    VTR_ASSERT(gain_clusters_[new_gain_cluster_id].molecules.size() == 1);
    VTR_ASSERT(gain_clusters_[new_gain_cluster_id].molecules[0] == mol);

    // Mark the molecule gain cluster as clustered.
    // NOTE: The previous insertion invalidates all references. Need to get a new
    //       reference.
    gain_clusters_[molecule_gain_cluster_id].valid = false;

    return new_gain_cluster_id;
}

void APClusterGainCalculator::add_mol_to_gain_cluster(const t_pack_molecule* mol,
                                                      APGainClusterId gain_cluster_id) {
    VTR_ASSERT(mol != nullptr);
    VTR_ASSERT(gain_cluster_id.is_valid() && (size_t)gain_cluster_id < gain_clusters_.size());

    // Get the current gain cluster to insert the molecule into.
    APGainCluster& gain_cluster = gain_clusters_[gain_cluster_id];

    // Get the cluster containing only the molecule.
    VTR_ASSERT(molecule_gain_cluster_.count(mol) != 0);
    APGainClusterId molecule_gain_cluster_id = molecule_gain_cluster_[mol];
    VTR_ASSERT(molecule_gain_cluster_id.is_valid());
    APGainCluster& molecule_gain_cluster = gain_clusters_[molecule_gain_cluster_id];

    // Insert the molecule into the gain cluster.
    gain_cluster.molecules.push_back(mol);
    // Mark this molecule as clustered.
    VTR_ASSERT(molecule_gain_cluster.valid);
    molecule_gain_cluster.valid = false;

    // Compute the new external pins.
    // Get the intersection of the two external pins. These are the shared
    // external pins between the two clusters.
    std::unordered_set<AtomPinId> shared_external_pins;
    std::set_intersection(gain_cluster.external_pins.begin(),
                          gain_cluster.external_pins.end(),
                          molecule_gain_cluster.external_pins.begin(),
                          molecule_gain_cluster.external_pins.end(),
                          std::inserter(shared_external_pins, shared_external_pins.begin()));
    // The new external pins is the union of the two external pins, without the
    // shared external pins.
    gain_cluster.external_pins.insert(molecule_gain_cluster.external_pins.begin(),
                                      molecule_gain_cluster.external_pins.end());
    for (AtomPinId shared_external_pin : shared_external_pins) {
        gain_cluster.external_pins.erase(shared_external_pin);
    }

    // Compute the new internal pins.
    // The new internal pins are simply the union of the internal pins of both
    // clusters and the shared external pins.
    gain_cluster.internal_pins.insert(molecule_gain_cluster.internal_pins.begin(),
                                      molecule_gain_cluster.internal_pins.end());
    gain_cluster.internal_pins.insert(shared_external_pins.begin(),
                                      shared_external_pins.end());

    // Compute the nets.
    // The only external nets which can become internal are the nets that are
    // shared between the two clusters.
    std::unordered_set<AtomNetId> shared_external_nets;
    std::set_intersection(gain_cluster.external_nets.begin(),
                          gain_cluster.external_nets.end(),
                          molecule_gain_cluster.external_nets.begin(),
                          molecule_gain_cluster.external_nets.end(),
                          std::inserter(shared_external_nets, shared_external_nets.begin()));
    // Combine the internal and external nets.
    gain_cluster.internal_nets.insert(molecule_gain_cluster.internal_nets.begin(),
                                      molecule_gain_cluster.internal_nets.end());
    gain_cluster.external_nets.insert(molecule_gain_cluster.external_nets.begin(),
                                      molecule_gain_cluster.external_nets.end());
    // Check if the shared external nets are now internal. If so, move them to
    // the correct set.
    for (AtomNetId net : shared_external_nets) {
        // FIXME: This is duplicate code.
        bool all_pins_internal = true;
        for (AtomPinId net_pin : atom_netlist_.net_pins(net)) {
            if (gain_cluster.internal_pins.count(net_pin) == 0) {
                all_pins_internal = false;
                break;
            }
        }
        if (all_pins_internal) {
            gain_cluster.external_nets.erase(net);
            gain_cluster.internal_nets.insert(net);
        }
    }


}

float APClusterGainCalculator::get_gain(APGainClusterId gain_cluster_id,
                                        const t_pack_molecule* mol) {
    VTR_ASSERT(mol != nullptr);
    VTR_ASSERT(gain_cluster_id.is_valid() && (size_t)gain_cluster_id < gain_clusters_.size());

    // Get the current gain cluster to insert the molecule into.
    APGainCluster& gain_cluster = gain_clusters_[gain_cluster_id];

    // Get the cluster containing only the molecule.
    VTR_ASSERT(molecule_gain_cluster_.count(mol) != 0);
    APGainClusterId molecule_gain_cluster_id = molecule_gain_cluster_[mol];
    VTR_ASSERT(molecule_gain_cluster_id.is_valid());
    APGainCluster& molecule_gain_cluster = gain_clusters_[molecule_gain_cluster_id];
    VTR_ASSERT(molecule_gain_cluster.valid);

    // Compute the Pin Sharing Gain as described by Singhal et al.
    float pin_gain = 0.f;
    size_t total_nets = 0;
    // Internal nets are simple. The gain increases by 1 for each internal net.
    pin_gain += gain_cluster.internal_nets.size();
    pin_gain += molecule_gain_cluster.internal_nets.size();
    total_nets += gain_cluster.internal_nets.size();
    total_nets += molecule_gain_cluster.internal_nets.size();
    // Collect the shared external nets. These will be handled special.
    std::unordered_set<AtomNetId> shared_external_nets;
    // External nets which are not shared.
    for (AtomNetId net : gain_cluster.external_nets) {
        if (molecule_gain_cluster.internal_nets.count(net) != 0) {
            continue;
            shared_external_nets.insert(net);
        }
        // Count the number of internal pins in this net.
        size_t num_internal_pins = 0;
        for (AtomPinId pin : atom_netlist_.net_pins(net)) {
            if (gain_cluster.internal_pins.count(pin) != 0) {
                num_internal_pins++;
            }
        }
        VTR_ASSERT(num_internal_pins != 0);
        pin_gain += static_cast<float>(num_internal_pins - 1) / static_cast<float>(atom_netlist_.net_pins(net).size() - 1);
        total_nets++;
    }
    // FIXME: Duplicate code!
    for (AtomNetId net : molecule_gain_cluster.external_nets) {
        if (gain_cluster.internal_nets.count(net) != 0) {
            continue;
            shared_external_nets.insert(net);
        }
        // Count the number of internal pins in this net.
        size_t num_internal_pins = 0;
        for (AtomPinId pin : atom_netlist_.net_pins(net)) {
            if (molecule_gain_cluster.internal_pins.count(pin) != 0) {
                num_internal_pins++;
            }
        }
        VTR_ASSERT(num_internal_pins != 0);
        pin_gain += static_cast<float>(num_internal_pins - 1) / static_cast<float>(atom_netlist_.net_pins(net).size() - 1);
        total_nets++;
    }
    // Shared external nets are handled specal.
    for (AtomNetId net : shared_external_nets) {
        // Count the number of internal pins in this net.
        size_t num_internal_pins = 0;
        for (AtomPinId pin : atom_netlist_.net_pins(net)) {
            if (molecule_gain_cluster.internal_pins.count(pin) != 0 ||
                gain_cluster.internal_pins.count(pin) != 0) {
                num_internal_pins++;
            }
        }
        VTR_ASSERT(num_internal_pins != 0);
        pin_gain += static_cast<float>(num_internal_pins - 1) / static_cast<float>(atom_netlist_.net_pins(net).size() - 1);
        total_nets++;
    }
    // Divide the gain by the total number of nets.
    // TODO: The LSC paper recommended dividing by 1 + c * total_nets where
    //       c is an empirically-determined constant. To investigate.
    pin_gain /= static_cast<float>(total_nets);

    // TODO: Implement the WL gain.
    float wl_gain = 0.f;

    // Return a weighted mean of the gain terms.
    float total_gain = (pin_gain * pin_gain_weight_) + (wl_gain * wl_gain_weight_);
    total_gain /= (pin_gain_weight_ + wl_gain_weight_);
    return total_gain;
}

void APClusterGainCalculator::destroy_gain_cluster(APGainClusterId gain_cluster_id) {
    VTR_ASSERT(gain_cluster_id.is_valid() && (size_t)gain_cluster_id < gain_clusters_.size());

    // Mark all of the molecules in this gain cluster as unclustered.
    APGainCluster& gain_cluster = gain_clusters_[gain_cluster_id];
    for (const t_pack_molecule* mol : gain_cluster.molecules) {
        VTR_ASSERT(molecule_gain_cluster_.count(mol) != 0);
        APGainClusterId molecule_gain_cluster_id = molecule_gain_cluster_[mol];
        VTR_ASSERT(molecule_gain_cluster_id.is_valid());
        APGainCluster& molecule_gain_cluster = gain_clusters_[molecule_gain_cluster_id];
        VTR_ASSERT(!molecule_gain_cluster.valid);
        molecule_gain_cluster.valid = true;
    }

    // Clean up all the data in the object. We will not remove it from the
    // vector for now. There is no reason to iterate over all of the gain
    // clusters, so we can tolerate some invalid members.
    gain_cluster.external_pins.clear();
    gain_cluster.internal_pins.clear();
    gain_cluster.external_nets.clear();
    gain_cluster.internal_nets.clear();
    gain_cluster.molecules.clear();
}

void APClusterGainCalculator::clean_gain_cluster(APGainClusterId gain_cluster_id) {
    VTR_ASSERT(gain_cluster_id.is_valid() && (size_t)gain_cluster_id < gain_clusters_.size());

    // Clean up all of the data for each molecule.
    APGainCluster& gain_cluster = gain_clusters_[gain_cluster_id];
    for (const t_pack_molecule* mol : gain_cluster.molecules) {
        VTR_ASSERT(molecule_gain_cluster_.count(mol) != 0);
        APGainClusterId molecule_gain_cluster_id = molecule_gain_cluster_[mol];
        VTR_ASSERT(molecule_gain_cluster_id.is_valid());
        APGainCluster& molecule_gain_cluster = gain_clusters_[molecule_gain_cluster_id];
        VTR_ASSERT(!molecule_gain_cluster.valid);
        molecule_gain_cluster.external_pins.clear();
        molecule_gain_cluster.internal_pins.clear();
        molecule_gain_cluster.external_nets.clear();
        molecule_gain_cluster.internal_nets.clear();
        molecule_gain_cluster.molecules.clear();
    }

    // Clean up all the data in the object. We will not remove it from the
    // vector for now. There is no reason to iterate over all of the gain
    // clusters, so we can tolerate some invalid members.
    gain_cluster.external_pins.clear();
    gain_cluster.internal_pins.clear();
    gain_cluster.external_nets.clear();
    gain_cluster.internal_nets.clear();
    gain_cluster.molecules.clear();
}

GreedyAPClusterer::GreedyAPClusterer(const APNetlist& netlist,
                                     t_vpr_setup& vpr_setup,
                                     const t_arch* arch,
                                     const AtomNetlist& atom_netlist,
                                     const Prepacker& prepacker,
                                     const std::vector<t_logical_block_type>& logical_block_types,
                                     std::vector<t_lb_type_rr_node>* lb_type_rr_graphs,
                                     const t_model* user_models,
                                     const t_model* library_models,
                                     const t_packer_opts& packer_opts)
   : APClusterer(netlist),
     high_fanout_thresholds_(packer_opts.high_fanout_threshold),
     vpr_setup_(vpr_setup),
     arch_(arch),
     atom_netlist_(atom_netlist),
     prepacker_(prepacker),
     packer_opts_(packer_opts),
     cluster_legalizer_(atom_netlist,
                        prepacker,
                        logical_block_types,
                        lb_type_rr_graphs,
                        user_models,
                        library_models,
                        packer_opts.target_external_pin_util,
                        high_fanout_thresholds_,
                        ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE,
                        packer_opts.enable_pin_feasibility_filter,
                        packer_opts.feasible_block_array_size,
                        packer_opts.pack_verbosity),
     primitive_candidate_block_types_(identify_primitive_candidate_block_types()) {

    // Initialize the molecule to block lookup.
    for (APBlockId block : netlist.blocks()) {
        mol_block_[netlist.block_molecule(block)] = block;
    }
}

APBlockId GreedyAPClusterer::select_seed_block() {
    // For now we are using a simple scheme that just finds the block with the
    // largest number of pins.
    // TODO: Papers recommend using the block with the largest number of nets.
    //       This is not the same.
    // TODO: This does not take into account the current clustering. It really
    //       should.
    // TODO: I wonder if we can do something with gain here. Just getting the
    //       highest gain APBlock. The gain calculator may have this for free!

    // FIXME: This should really not be an int, it should be a size_t; but need
    //        to handle the case where the number of pins is zero.
    int best_num_pins = -1;
    APBlockId best_unclustered_block;
    for (APBlockId block : netlist_.blocks()) {
        int num_block_pins = netlist_.block_pins(block).size();
        // FIXME: The netlist stores a const pointer to mol; but the cluster
        //        legalizer does not accept this. Need to fix one or the other.
        t_pack_molecule* block_molecule = const_cast<t_pack_molecule*>(netlist_.block_molecule(block));
        if (!cluster_legalizer_.is_mol_clustered(block_molecule) && num_block_pins > best_num_pins) {
            best_unclustered_block = block;
            best_num_pins = num_block_pins;
        }
    }

    // Ensure that a best was found, or the best unclustered block is invalid.
    VTR_ASSERT(best_num_pins != -1 || !best_unclustered_block.is_valid());

    // If no unclustered block can be found, then return invalid block.
    return best_unclustered_block;
}

APBlockId GreedyAPClusterer::get_highest_gain_compatible_neighbor(APGainClusterId gain_cluster_id,
                                                                  LegalizationClusterId leg_cluster_id) {
    // TODO: Should maintain a netlist where we can quickly lookup the neighbors.
    // This should be abstracted into the gain calculator class. The greedy
    // clusterer would just ask which has the lowest gain. This can be pre-computed
    // to save lot of time.

    // FIXME: We should ignore high-fanout nets!

    // For now just get the closest neighbor. This is for testing code.
    // FIXME: Should precompute the neighbors. This can probably be done in the
    //        gain calculator.
    APBlockId highest_gain_neighbor;
    float highest_gain = std::numeric_limits<float>::lowest();
    const std::vector<t_pack_molecule*>& cluster_molecules = cluster_legalizer_.get_cluster_molecules(leg_cluster_id);
    for (const t_pack_molecule* molecule : cluster_molecules) {
        APBlockId block = mol_block_[molecule];

        for (APPinId pin : netlist_.block_pins(block)) {
            APNetId net = netlist_.pin_net(pin);
            
            // FIXME: This should really be on the atom netlist not AP netlist.
            if (netlist_.net_pins(net).size() > 50)
                continue;

            for (APPinId neighbor_pin : netlist_.net_pins(net)) {
                APBlockId neighbor_block = netlist_.pin_block(neighbor_pin);
                if (neighbor_block == block)
                    continue;
                t_pack_molecule* neighbor_mol = const_cast<t_pack_molecule*>(netlist_.block_molecule(neighbor_block));
                if (cluster_legalizer_.is_mol_clustered(neighbor_mol))
                    continue;
                if (!cluster_legalizer_.is_molecule_compatible(neighbor_mol, leg_cluster_id))
                    continue;
                float gain = gain_calculator_->get_gain(gain_cluster_id, neighbor_mol);
                if (gain > highest_gain) {
                    highest_gain_neighbor = neighbor_block;
                    highest_gain = gain;
                }
            }
        }
    }

    VTR_ASSERT(highest_gain != std::numeric_limits<float>::lowest() ||
               !highest_gain_neighbor.is_valid());

    return highest_gain_neighbor;
}

bool GreedyAPClusterer::add_block_to_cluster(APBlockId block,
                                             LegalizationClusterId cluster_id) {
    VTR_ASSERT_DEBUG(block.is_valid());
    VTR_ASSERT_DEBUG(cluster_id.is_valid());

    // FIXME: The netlist stores a const pointer to mol; but the cluster
    //        legalizer does not accept this. Need to fix one or the other.
    t_pack_molecule* block_molecule = const_cast<t_pack_molecule*>(netlist_.block_molecule(block));

    e_block_pack_status pack_status = cluster_legalizer_.add_mol_to_cluster(block_molecule, cluster_id);
    
    return pack_status == e_block_pack_status::BLK_PASSED;
}

LegalizationClusterId GreedyAPClusterer::start_new_cluster(APBlockId seed) {
    VTR_ASSERT_DEBUG(seed.is_valid());

    // FIXME: The netlist stores a const pointer to mol; but the cluster
    //        legalizer does not accept this. Need to fix one or the other.
    t_pack_molecule* seed_molecule = const_cast<t_pack_molecule*>(netlist_.block_molecule(seed));
    AtomBlockId root_atom = seed_molecule->atom_block_ids[seed_molecule->root];
    const t_model* root_model = atom_netlist_.block_model(root_atom);

    auto itr = primitive_candidate_block_types_.find(root_model);
    VTR_ASSERT(itr != primitive_candidate_block_types_.end());
    const std::vector<t_logical_block_type_ptr>& candidate_types = itr->second;
    // FIXME: Sort based on balance!

    for (t_logical_block_type_ptr type : candidate_types) {
        int num_modes = type->pb_graph_head->pb_type->num_modes;
        for (int mode = 0; mode < num_modes; mode++) {
            e_block_pack_status pack_status = e_block_pack_status::BLK_STATUS_UNDEFINED;
            LegalizationClusterId new_cluster_id;
            std::tie(pack_status, new_cluster_id) = cluster_legalizer_.start_new_cluster(seed_molecule, type, mode);
            if (pack_status == e_block_pack_status::BLK_PASSED) {
                VTR_ASSERT_DEBUG(new_cluster_id.is_valid());
                return new_cluster_id;
            }
        }
    }

    // This should never happen.
    VPR_FATAL_ERROR(VPR_ERROR_AP,
                    "Unable to create a cluster for the given AP block");
    return LegalizationClusterId();
}

bool GreedyAPClusterer::grow_cluster(APBlockId seed,
                                     ClusterLegalizationStrategy strategy) {
    VTR_ASSERT(seed.is_valid());

    // Set the cluster legalizer to the provided strategy.
    cluster_legalizer_.set_legalization_strategy(strategy);

    // Start a new cluster containing just the seed block.
    LegalizationClusterId new_cluster_id = start_new_cluster(seed);
    VTR_ASSERT(new_cluster_id.is_valid());

    // Start the new cluster in the gain calculator.
    APGainClusterId gain_cluster_id = gain_calculator_->create_gain_cluster(netlist_.block_molecule(seed));

    // Add blocks to the cluster.
    while (true) {
        // Get the highest gain compatible neighbor as a candidate to add to the
        // cluster. If one cannot be found, break.
        APBlockId candidate_block = get_highest_gain_compatible_neighbor(gain_cluster_id, new_cluster_id);
        if (!candidate_block.is_valid())
            break;

        // Try to add the candidate block to the cluster. If it cannot be
        // legally added, break.
        bool add_success = add_block_to_cluster(candidate_block, new_cluster_id);
        if (!add_success)
            break;

        gain_calculator_->add_mol_to_gain_cluster(netlist_.block_molecule(candidate_block), gain_cluster_id);
    }

    // If the cluster legalization strategy skipped intra LB routing, need to
    // check intra LB routing at the end. If it fails, destroy the cluster.
    if (strategy == ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE) {
        bool cluster_is_legal = cluster_legalizer_.check_cluster_legality(new_cluster_id);
        if (!cluster_is_legal) {
            gain_calculator_->destroy_gain_cluster(gain_cluster_id);
            cluster_legalizer_.destroy_cluster(new_cluster_id);
            cluster_legalizer_.compress();
            return false;
        }
    }

    gain_calculator_->clean_gain_cluster(gain_cluster_id);

    // Clean the cluster. Since we know we will never modify the cluster passed
    // this point, we can remove some of the data from the legalizer to save
    // memory.
    cluster_legalizer_.clean_cluster(new_cluster_id);

    // Cluster has grown successfully.
    return true;
}

void GreedyAPClusterer::create_clusters(const PartialPlacement& p_placement) {
    // Ensure that the cluster legalizer is empty. Cannot run the clusterer
    // more than once.
    VTR_ASSERT(cluster_legalizer_.clusters().empty());

    // TODO: The partial placement should be fed into the gain calculator once
    //       it is created. It may need to be passed into grow_cluster or
    //       select seed_block.
    gain_calculator_ = std::make_unique<APClusterGainCalculator>(atom_netlist_,
                                                                 prepacker_,
                                                                 p_placement);

    // Select the starting seed block.
    APBlockId seed_block = select_seed_block();
    // Greedily form the clusters.
    while (seed_block.is_valid()) {
        // Grow the cluster by only performing intra LB routing at end.
        bool success = grow_cluster(seed_block,
                                    ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE);
        if (!success) {
            // If growing the cluster failed (meaning that the final intra LB
            // routing failed), try to grow the cluster again but run intra LB
            // routing for every molecule that is inserted.
            success = grow_cluster(seed_block,
                                   ClusterLegalizationStrategy::FULL);
            // Some cluster must be created. Even a cluster of just 1 molecule
            // should always exist.
            VTR_ASSERT(success);
        }

        // Select the next seed block.
        seed_block = select_seed_block();
    }

    // Destory the gain calculator. No longer needed.
    gain_calculator_.reset();

    // Check and output the clustering
    std::unordered_set<AtomNetId> is_clock = alloc_and_load_is_clock();
    check_and_output_clustering(cluster_legalizer_, packer_opts_, is_clock, arch_);

    // Reset the cluster legalizer. This is required to load the packing.
    cluster_legalizer_.reset();
    // Regenerate the clustered netlist from the file generated previously.
    // FIXME: This writing and loading from a file is wasteful. Should generate
    //        the clusters directly from the cluster legalizer.
    vpr_load_packing(vpr_setup_, *arch_);
    load_cluster_constraints();
    const ClusteredNetlist& clb_nlist = g_vpr_ctx.clustering().clb_nlist;

    // Verify the packing and print some info.
    check_netlist(packer_opts_.pack_verbosity);
    writeClusteredNetlistStats(vpr_setup_.FileNameOpts.write_block_usage);
    print_pb_type_count(clb_nlist);
}

