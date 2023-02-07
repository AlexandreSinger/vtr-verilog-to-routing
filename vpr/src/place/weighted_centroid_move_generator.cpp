#include "weighted_centroid_move_generator.h"
#include "globals.h"
#include "directed_moves_util.h"
#include "place_constraints.h"
#include "move_utils.h"

e_create_move WeightedCentroidMoveGenerator::propose_move(t_pl_blocks_to_be_moved& blocks_affected, e_move_type& /*move_type*/, t_logical_block_type& blk_type, float rlim, const t_placer_opts& placer_opts, const PlacerCriticalities* criticalities) {
    ClusterBlockId b_from;
    auto& cluster_ctx = g_vpr_ctx.clustering();
    if (blk_type.index == -1) { //If the block type is unspecified, choose any random block to be swapped with another random block
        b_from = pick_from_block();
        if (!b_from) {
            return e_create_move::ABORT; //No movable block found
        }
        blk_type.index = convert_logical_to_agent_block_type(cluster_ctx.clb_nlist.block_type(b_from)->index);
    } else { //If the block type is specified, choose a random block with blk_type to be swapped with another random block
        b_from = pick_from_block(blk_type);
    }

    auto& place_ctx = g_vpr_ctx.placement();
    auto& device_ctx = g_vpr_ctx.device();

    auto& place_move_ctx = g_placer_ctx.mutable_move();

    t_pl_loc from = place_ctx.block_locs[b_from].loc;
    auto cluster_from_type = cluster_ctx.clb_nlist.block_type(b_from);
    auto grid_from_type = device_ctx.grid[from.x][from.y].type;
    VTR_ASSERT(is_tile_compatible(grid_from_type, cluster_from_type));

    t_range_limiters range_limiters;
    range_limiters.original_rlim = rlim;
    range_limiters.first_rlim = place_move_ctx.first_rlim;
    range_limiters.dm_rlim = placer_opts.place_dm_rlim;

    t_pl_loc to, centroid;

    /* Calculate the weighted centroid */
    calculate_centroid_loc(b_from, true, centroid, criticalities);

    /* Find a  */
    if (!find_to_loc_centroid(cluster_from_type, from, centroid, range_limiters, to, b_from)) {
        return e_create_move::ABORT;
    }

    e_create_move create_move = ::create_move(blocks_affected, b_from, to);

    //Check that all of the blocks affected by the move would still be in a legal floorplan region after the swap
    if (!floorplan_legal(blocks_affected)) {
        return e_create_move::ABORT;
    }

    return create_move;
}
