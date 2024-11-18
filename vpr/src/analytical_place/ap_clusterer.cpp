
#include "ap_clusterer.h"
#include <limits>
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

APBlockId GreedyAPClusterer::get_highest_gain_compatible_neighbor(LegalizationClusterId cluster_id,
                                                                  const PartialPlacement& p_placement) {
    // TODO: Should maintain a netlist where we can quickly lookup the neighbors.
    // This should be abstracted into the gain calculator class. The greedy
    // clusterer would just ask which has the lowest gain. This can be pre-computed
    // to save lot of time.

    // For now just get the closest neighbor. This is for testing code.
    APBlockId closest_neighbor;
    double closest_neighbor_distance = std::numeric_limits<double>::max();
    const std::vector<t_pack_molecule*>& cluster_molecules = cluster_legalizer_.get_cluster_molecules(cluster_id);
    for (const t_pack_molecule* molecule : cluster_molecules) {
        APBlockId block = mol_block_[molecule];

        for (APPinId pin : netlist_.block_pins(block)) {
            APNetId net = netlist_.pin_net(pin);
            for (APPinId neighbor_pin : netlist_.net_pins(net)) {
                APBlockId neighbor_block = netlist_.pin_block(neighbor_pin);
                if (neighbor_block == block)
                    continue;
                t_pack_molecule* neighbor_mol = const_cast<t_pack_molecule*>(netlist_.block_molecule(neighbor_block));
                if (cluster_legalizer_.is_mol_clustered(neighbor_mol))
                    continue;
                if (!cluster_legalizer_.is_molecule_compatible(neighbor_mol, cluster_id))
                    continue;
                // FIXME: This should really be the distance to the centroid.
                double dx = p_placement.block_x_locs[neighbor_block] - p_placement.block_x_locs[block];
                double dy = p_placement.block_y_locs[neighbor_block] - p_placement.block_y_locs[block];
                double dist = (dx * dx) + (dy * dy);
                if (dist < closest_neighbor_distance) {
                    closest_neighbor = neighbor_block;
                    closest_neighbor_distance = dist;
                }
            }
        }
    }

    VTR_ASSERT(closest_neighbor_distance != std::numeric_limits<double>::max() ||
               !closest_neighbor.is_valid());

    return closest_neighbor;
}

void GreedyAPClusterer::destroy_cluster(LegalizationClusterId cluster_id) {
    VTR_ASSERT_DEBUG(cluster_id.is_valid());
    cluster_legalizer_.destroy_cluster(cluster_id);
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
                                     const PartialPlacement& p_placement,
                                     ClusterLegalizationStrategy strategy) {
    VTR_ASSERT(seed.is_valid());

    // Set the cluster legalizer to the provided strategy.
    cluster_legalizer_.set_legalization_strategy(strategy);

    // Start a new cluster containing just the seed block.
    LegalizationClusterId new_cluster_id = start_new_cluster(seed);
    VTR_ASSERT(new_cluster_id.is_valid());

    // Add blocks to the cluster.
    while (true) {
        // Get the highest gain compatible neighbor as a candidate to add to the
        // cluster. If one cannot be found, break.
        APBlockId candidate_block = get_highest_gain_compatible_neighbor(new_cluster_id,
                                                                         p_placement);
        if (!candidate_block.is_valid())
            break;

        // Try to add the candidate block to the cluster. If it cannot be
        // legally added, break.
        bool add_success = add_block_to_cluster(candidate_block, new_cluster_id);
        if (!add_success)
            break;
    }

    // If the cluster legalization strategy skipped intra LB routing, need to
    // check intra LB routing at the end. If it fails, destroy the cluster.
    if (strategy == ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE) {
        bool cluster_is_legal = cluster_legalizer_.check_cluster_legality(new_cluster_id);
        if (!cluster_is_legal) {
            destroy_cluster(new_cluster_id);
            cluster_legalizer_.compress();
            return false;
        }
    }

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

    // Select the starting seed block.
    APBlockId seed_block = select_seed_block();
    // Greedily form the clusters.
    while (seed_block.is_valid()) {
        // Grow the cluster by only performing intra LB routing at end.
        bool success = grow_cluster(seed_block,
                                    p_placement,
                                    ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE);
        if (!success) {
            // If growing the cluster failed (meaning that the final intra LB
            // routing failed), try to grow the cluster again but run intra LB
            // routing for every molecule that is inserted.
            success = grow_cluster(seed_block,
                                   p_placement,
                                   ClusterLegalizationStrategy::FULL);
            // Some cluster must be created. Even a cluster of just 1 molecule
            // should always exist.
            VTR_ASSERT(success);
        }

        // Select the next seed block.
        seed_block = select_seed_block();
    }

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

