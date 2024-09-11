
#include "full_legalizer.h"
#include <cmath>
#include <list>
#include <unordered_set>
#include <vector>
#include "PartialPlacement.h"
#include "ShowSetup.h"
#include "ap_netlist_fwd.h"
#include "check_netlist.h"
#include "cluster.h"
#include "cluster_legalizer.h"
#include "cluster_util.h"
#include "clustered_netlist.h"
#include "globals.h"
#include "initial_placement.h"
#include "logic_types.h"
#include "net_cost_handler.h"
#include "pack.h"
#include "physical_types.h"
#include "place_constraints.h"
#include "vpr_api.h"
#include "vpr_context.h"
#include "vpr_types.h"
#include "vtr_assert.h"
#include "vtr_ndmatrix.h"
#include "vtr_strong_id.h"
#include "vtr_vector.h"

namespace {

// FIXME: Bring this into the device_grid.
//  - Maybe this can be called DeviceRootTileId or something.
struct device_tile_id_tag {};
typedef vtr::StrongId<device_tile_id_tag, size_t> DeviceTileId;

class APClusterPlacer {
private:
    // Get the macro for the given cluster block.
    t_pl_macro get_macro(ClusterBlockId clb_blk_id) {
        // Basically stolen from initial_placement.cpp:place_one_block
        // TODO: Make this a cleaner interface and share the code.
        int imacro;
        get_imacro_from_iblk(&imacro, clb_blk_id, g_vpr_ctx.placement().pl_macros);
        // If this block is part of a macro, return it.
        if (imacro != -1)
            return g_vpr_ctx.placement().pl_macros[imacro];
        // If not, create a "fake" macro with a single element.
        t_pl_macro_member macro_member;
        t_pl_offset block_offset(0, 0, 0, 0);
        macro_member.blk_index = clb_blk_id;
        macro_member.offset = block_offset;
        
        t_pl_macro pl_macro;
        pl_macro.members.push_back(macro_member);
        return pl_macro;
    }

public:
    APClusterPlacer() {
        // FIXME: This was stolen from place/place.cpp
        //        it used a static method, just taking what I think I will need.
        // FIXME: WILL LIKELY NEED MORE! DID NOT ALLOCATE PLACEMENT MACROS!

        auto& block_locs = g_vpr_ctx.mutable_placement().mutable_block_locs();
        auto& grid_blocks = g_vpr_ctx.mutable_placement().mutable_grid_blocks();
        auto& blk_loc_registry = g_vpr_ctx.mutable_placement().mutable_blk_loc_registry();
        init_placement_context(block_locs, grid_blocks);

        // stolen from place/place.cpp:alloc_and_load_try_swap_structs
        size_t num_nets = g_vpr_ctx.clustering().clb_nlist.nets().size();
        // FIXME: set cube_bb to false by hand, should be passed in.
        g_vpr_ctx.mutable_placement().cube_bb = false;
        init_try_swap_net_cost_structs(num_nets, false /*cube_bb*/);
        g_vpr_ctx.mutable_placement().compressed_block_grids = create_compressed_block_grids();

        // Initialize the macros
        const t_arch* arch = g_vpr_ctx.device().arch;
        g_vpr_ctx.mutable_placement().pl_macros = alloc_and_load_placement_macros(arch->Directs, arch->num_directs);

        // TODO: The next few steps will be basically a direct copy of the initial
        //       placement code since it does everything we need! It would be nice
        //       to share the code.

        // Clear the grid locations (stolen from initial_placement)
        // FIXME: Should I have stole this?
        VTR_LOG("CLEARING GRID LOCS\n");
        clear_all_grid_locs(blk_loc_registry);

        // Deal with the placement constraints.
        VTR_LOG("Propogating constraints\n");
        propagate_place_constraints();

        VTR_LOG("Marking fixed blocks\n");
        mark_fixed_blocks(blk_loc_registry);

        VTR_LOG("Allocating and loading compressed cluster constraints\n");
        alloc_and_load_compressed_cluster_constraints();    
    }

    // Given a cluster and tile it wants to go into, try to place the cluster
    // at this tile's postion.
    bool place_cluster(ClusterBlockId clb_blk_id, const t_physical_tile_loc& tile_loc, int sub_tile) {
        const DeviceContext& device_ctx = g_vpr_ctx.device();
        // FIXME: THIS MUST TAKE INTO ACCOUNT THE CONSTRAINTS AS WELL!!!
        //  - Right now it is just implied.
        const auto& block_locs = g_vpr_ctx.placement().block_locs();
        auto& blk_loc_registry = g_vpr_ctx.mutable_placement().mutable_blk_loc_registry();
        VTR_ASSERT(!is_block_placed(clb_blk_id, block_locs) && "Block already placed. Is this intentional?");
        t_pl_macro pl_macro = get_macro(clb_blk_id);
        t_pl_loc to_loc;
        to_loc.x = tile_loc.x;
        to_loc.y = tile_loc.y;
        to_loc.layer = tile_loc.layer_num;
        // Special case where the tile has no sub-tiles. It just cannot be placed.
        if (device_ctx.grid.get_physical_type(tile_loc)->sub_tiles.size() == 0)
            return false;
        VTR_ASSERT(sub_tile >= 0 && sub_tile < device_ctx.grid.get_physical_type(tile_loc)->capacity);
        // FIXME: GET THIS FROM THE PARTIAL PLACEMENT!!!! Or do this better.
        //  - May need to try all the sub-tiles in a location.
        //  - https://github.com/AlexandreSinger/vtr-verilog-to-routing/blob/feature-analytical-placer/vpr/src/place/initial_placement.cpp#L755
        to_loc.sub_tile = sub_tile;
        // FIXME: NEED TO VERIFY THAT THIS FOLLOWS THE PLACEMENT CONSTRAINTS!
        //          - CANNOT PLACE A CLUSTER IN A PLACE IT CANNOT EXIST.
        return try_place_macro(pl_macro, to_loc, blk_loc_registry);
    }

    // This is not the best way of doing things, but its the simplest. Given a
    // cluster, just find somewhere for it to go.
    // TODO: Make this like the initial placement code where we first try
    //       centroid, then random, then exhaustive.
    bool exhaustively_place_cluster(ClusterBlockId clb_blk_id) {
        // TODO: Again, try_place_macro_randomly/exhaustively may be very very
        //       helpful here.
        const auto& block_locs = g_vpr_ctx.placement().block_locs();
        auto& blk_loc_registry = g_vpr_ctx.mutable_placement().mutable_blk_loc_registry();
        VTR_ASSERT(!is_block_placed(clb_blk_id, block_locs) && "Block already placed. Is this intentional?");
        t_pl_macro pl_macro = get_macro(clb_blk_id);
        const PartitionRegion& pr = is_cluster_constrained(clb_blk_id) ? g_vpr_ctx.floorplanning().cluster_constraints[clb_blk_id] : get_device_partition_region();
        t_logical_block_type_ptr block_type = g_vpr_ctx.clustering().clb_nlist.block_type(clb_blk_id);
        // FIXME: We really should get this from the place context, not the device context.
        //      - Stealing it for now to get this to work.
        enum e_pad_loc_type pad_loc_type = g_vpr_ctx.device().pad_loc_type;
        return try_place_macro_exhaustively(pl_macro, pr, block_type, pad_loc_type, blk_loc_registry);
    }
};
} // namespace

// FIXME: This was stolen from pack.cpp; should probably be absorbed into the
//        cluster legalizer.
static size_t count_models(const t_model* models) {
    if (models == nullptr)
        return 0;

    size_t num_models = 0;
    const t_model* cur_model = models;
    while (cur_model != nullptr) {
        num_models++;
        cur_model = cur_model->next;
    }
    return num_models;
}

static LegalizationClusterId create_new_cluster(t_pack_molecule* seed_molecule,
                                                ClusterLegalizer& cluster_legalizer,
                                                const std::map<const t_model*, std::vector<t_logical_block_type_ptr>>& primitive_candidate_block_types) {
    const AtomContext& atom_ctx = g_vpr_ctx.atom();
    // This was stolen from pack/cluster_util.cpp:start_new_cluster
    // It tries to find a block type and mode for the given molecule.
    // TODO: This should take into account the tile this molecule wants to be
    //       placed into.
    // TODO: The original implementation sorted based on balance. Perhaps this
    //       should do the same.
    AtomBlockId root_atom = seed_molecule->atom_block_ids[seed_molecule->root];
    const t_model* root_model = atom_ctx.nlist.block_model(root_atom);

    auto itr = primitive_candidate_block_types.find(root_model);
    VTR_ASSERT(itr != primitive_candidate_block_types.end());
    const std::vector<t_logical_block_type_ptr>& candidate_types = itr->second;

    for (t_logical_block_type_ptr type : candidate_types) {
        int num_modes = type->pb_graph_head->pb_type->num_modes;
        for (int mode = 0; mode < num_modes; mode++) {
            e_block_pack_status pack_status = e_block_pack_status::BLK_STATUS_UNDEFINED;
            LegalizationClusterId new_cluster_id;
            std::tie(pack_status, new_cluster_id) = cluster_legalizer.start_new_cluster(seed_molecule, type, mode);
            if (pack_status == e_block_pack_status::BLK_PASSED)
                return new_cluster_id;
        }
    }
    // This should never happen.
    // FIXME: Turn this into an error.
    VTR_ASSERT(false);
    return LegalizationClusterId();
}

// TODO: Organize this better by splitting into two sub-modules: Pack and Place?
//  - Perhaps we may want to combine the two.
void FullLegalizer::legalize(PartialPlacement& p_placement) {
    // PACKING:
    // Initialize the cluster legalizer (Packing)
    
    // FIXME: The legalization strategy is currently set to full. Should handle
    //        this better to make it faster.
    size_t num_models = count_models(user_models_) + count_models(library_models_);
    t_pack_high_fanout_thresholds high_fanout_thresholds(packer_opts_.high_fanout_threshold);
    ClusterLegalizer cluster_legalizer(atom_netlist_,
                                       prepacker_,
                                       logical_block_types_,
                                       lb_type_rr_graphs_,
                                       num_models,
                                       packer_opts_.target_external_pin_util,
                                       high_fanout_thresholds,
                                       ClusterLegalizationStrategy::FULL,
                                       packer_opts_.enable_pin_feasibility_filter,
                                       packer_opts_.feasible_block_array_size,
                                       packer_opts_.pack_verbosity);
    // Create clusters for each tile.
    //  Start by giving each root tile a unique ID.
    size_t grid_width = device_grid_.width();
    size_t grid_height = device_grid_.height();
    vtr::NdMatrix<DeviceTileId, 2> tile_grid({grid_width, grid_height});
    size_t num_device_tiles = 0;
    for (size_t x = 0; x < grid_width; x++) {
        for (size_t y = 0; y < grid_height; y++) {
            // Ignoring 3D placement for now.
            t_physical_tile_loc tile_loc(x, y, 0);
            // Ignore non-root locations
            size_t width_offset = device_grid_.get_width_offset(tile_loc);
            size_t height_offset = device_grid_.get_height_offset(tile_loc);
            if (width_offset != 0 || height_offset != 0) {
                tile_grid[x][y] = tile_grid[x - width_offset][y - height_offset];
                continue;
            }
            tile_grid[x][y] = DeviceTileId(num_device_tiles);
            num_device_tiles++;
        }
    }
    //  Next, collect the AP blocks which will go into each root tile
    VTR_ASSERT_SAFE(p_placement.verify_locs(ap_netlist_, g_vpr_ctx.device()));
    vtr::vector<DeviceTileId, std::vector<APBlockId>> blocks_in_tiles(num_device_tiles);
    for (APBlockId ap_blk_id : ap_netlist_.blocks()) {
        // FIXME: Add these conversions to the PartialPlacement class.
        size_t blk_x_loc = std::floor(p_placement.block_x_locs[ap_blk_id]);
        size_t blk_y_loc = std::floor(p_placement.block_y_locs[ap_blk_id]);
        VTR_ASSERT(p_placement.block_layer_nums[ap_blk_id] == 0);
        DeviceTileId tile_id = tile_grid[blk_x_loc][blk_y_loc];
        blocks_in_tiles[tile_id].push_back(ap_blk_id);
    }
    //  Create the legalized clusters per tile.
    std::map<const t_model*, std::vector<t_logical_block_type_ptr>> primitive_candidate_block_types = identify_primitive_candidate_block_types();

    for (size_t tile_id_idx = 0; tile_id_idx < num_device_tiles; tile_id_idx++) {
        DeviceTileId tile_id = DeviceTileId(tile_id_idx);
        // Create the molecule list
        std::list<t_pack_molecule*> mol_list;
        for (APBlockId ap_blk_id : blocks_in_tiles[tile_id]) {
            // FIXME: The netlist stores a const pointer to mol; but the cluster
            //        legalizer does not accept this. Need to fix one or the other.
            // For now, using const_cast cause lazy.
            t_pack_molecule* mol = const_cast<t_pack_molecule*>(ap_netlist_.block_molecule(ap_blk_id));
            mol_list.push_back(mol);
        }
        // Clustering algorithm: Create clusters one at a time.
        while (!mol_list.empty()) {
            // Arbitrarily choose the first molecule as a seed molecule.
            t_pack_molecule* seed_mol = mol_list.front();
            mol_list.pop_front();
            // Use the seed molecule to create a cluster for this tile.
            LegalizationClusterId new_cluster_id = create_new_cluster(seed_mol, cluster_legalizer, primitive_candidate_block_types);
            // Insert all molecules that you can into the cluster.
            // NOTE: If the mol_list was somehow sorted, we can just stop at
            //       first failure!
            auto it = mol_list.begin();
            while (it != mol_list.end()) {
                t_pack_molecule* mol = *it;
                // FIXME: IS THIS NEEDED?
                if (!check_free_primitives_for_molecule_atoms(mol, cluster_legalizer.get_cluster_placement_stats(new_cluster_id), cluster_legalizer)) {
                    ++it;
                    continue;
                }
                // Try to insert it. If successful, remove from list.
                e_block_pack_status pack_status = cluster_legalizer.add_mol_to_cluster(mol, new_cluster_id);
                if (pack_status == e_block_pack_status::BLK_PASSED) {
                    it = mol_list.erase(it);
                } else {
                    ++it;
                }
            }
            // Once all molecules have been inserted, clean the cluster.
            // NOTE: Due to current limitations on the cluster legalizer, cannot
            //       build two clusters at the same time. This will be fixed later.
            cluster_legalizer.clean_cluster(new_cluster_id);
        }
    }

    // Check and output the clustering.
    // FIXME: This should be done in the constructor or somewhere safer.
    std::unordered_set<AtomNetId> is_clock = alloc_and_load_is_clock();
    check_and_output_clustering(cluster_legalizer, packer_opts_, is_clock, arch_);
    // Reset the cluster legalizer. This is required to load the packing.
    cluster_legalizer.reset();
    // Regenerate the clustered netlist from the file generated previously.
    vpr_load_packing(vpr_setup_, *arch_);
    load_cluster_constraints();
    const ClusteredNetlist& clb_nlist = g_vpr_ctx.clustering().clb_nlist;

    // Verify the packing and print some info
    check_netlist(packer_opts_.pack_verbosity);
    writeClusteredNetlistStats(vpr_setup_.FileNameOpts.write_block_usage);
    print_pb_type_count(clb_nlist);

    // PLACING:
    // Create a lookup from the AtomBlockId to the APBlockId
    vtr::vector<AtomBlockId, APBlockId> atom_to_ap_block(atom_netlist_.blocks().size());
    for (APBlockId ap_blk_id : ap_netlist_.blocks()) {
        const t_pack_molecule* blk_mol = ap_netlist_.block_molecule(ap_blk_id);
        for (AtomBlockId atom_blk_id : blk_mol->atom_block_ids) {
            // Ensure that this block is not in any other AP block. That would
            // be weird.
            VTR_ASSERT(!atom_to_ap_block[atom_blk_id].is_valid());
            atom_to_ap_block[atom_blk_id] = ap_blk_id;
        }
    }
    // Move the clusters to legal positions
    APClusterPlacer ap_cluster_placer;
    VTR_LOG("Placing clusters...\n");
    std::vector<ClusterBlockId> unplaced_clusters;
    for (ClusterBlockId cluster_blk_id : clb_nlist.blocks()) {
        // Assume that the cluster will always want to be placed wherever the
        // first atom in the cluster wants to be placed.
        // FIXME: This assumption does not always hold! Will need to unify the
        //        cluster legalizer and the clustered netlist!
        const std::unordered_set<AtomBlockId>& atoms_in_cluster = g_vpr_ctx.clustering().atoms_lookup[cluster_blk_id];
        VTR_ASSERT(atoms_in_cluster.size() > 0);
        AtomBlockId first_atom_blk = *atoms_in_cluster.begin();
        APBlockId first_ap_blk = atom_to_ap_block[first_atom_blk];
        size_t blk_x_loc = p_placement.block_x_locs[first_ap_blk];
        size_t blk_y_loc = p_placement.block_y_locs[first_ap_blk];
        size_t blk_layer_num = p_placement.block_layer_nums[first_ap_blk];
        size_t blk_sub_tile = p_placement.block_sub_tiles[first_ap_blk];
        t_physical_tile_loc tile_loc(blk_x_loc, blk_y_loc,blk_layer_num);
        bool placed = ap_cluster_placer.place_cluster(cluster_blk_id, tile_loc, blk_sub_tile);
        if (placed)
            continue;
        // FIXME: Should now try all sub-tiles at this tile location.
        //  - May need to try all the sub-tiles in a location.
        //  - https://github.com/AlexandreSinger/vtr-verilog-to-routing/blob/feature-analytical-placer/vpr/src/place/initial_placement.cpp#L755
        
        // Add to list of unplaced clusters.
        unplaced_clusters.push_back(cluster_blk_id);
    }

    VTR_LOG("Placing unplaced clusters...\n");
    for (ClusterBlockId clb_blk_id : unplaced_clusters) {
        bool success = ap_cluster_placer.exhaustively_place_cluster(clb_blk_id);
        if (!success) {
            throw vtr::VtrError("Unable to find valid place for cluster in AP placement!");
        }
    }

    // Print some statistics about what happened here. This will be useful to
    // improve other algorithms.
    VTR_LOG("Number of clusters which needed to be moved: %zu\n", unplaced_clusters.size());

    // TODO: Print a breakdown per block type. We may find that specific block
    //       types are always conflicting.
    
    // FIXME: Allocate and load moveable blocks?
    //      - This may be needed to perform SA. Not needed right now.

    // FIXME: Check initial placement legality?
}

