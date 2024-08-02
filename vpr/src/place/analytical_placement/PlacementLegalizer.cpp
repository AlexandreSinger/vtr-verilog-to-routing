
#include "PlacementLegalizer.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "PartialPlacement.h"
#include "atom_netlist.h"
#include "cluster_util.h"
#include "clustered_netlist_utils.h"
#include "device_grid.h"
#include "globals.h"
#include "physical_types.h"
#include "place_util.h"
#include "re_cluster.h"
#include "vpr_context.h"
#include "vpr_types.h"
#include "vtr_assert.h"
#include "vtr_list.h"
#include "vtr_log.h"

namespace {

// A volume of primitives. A P-dimensional quantity (where P is the total
// number of unieque primitives) which is used to specify the capacity of a
// tile and the "size" of a molecule in primitives.
struct PrimitiveVolume {
    // Internally the primitive volume is a vector.
    // FIXME: This vector is likely very sparse; may be more efficient to
    //        store as a list.
    std::vector<float> primitive_lengths;
    PrimitiveVolume(size_t num_unique_primitives) : primitive_lengths(num_unique_primitives, 0.f) {}

    void add_length(size_t primitive_id, float amount) {
        primitive_lengths[primitive_id] += amount;
    }

    float get_length(size_t primitive_id) const {
        return primitive_lengths[primitive_id];
    }

    PrimitiveVolume& operator+=(const PrimitiveVolume& other) {
        size_t num_primitives = primitive_lengths.size();
        for (size_t i = 0; i < num_primitives; i++) {
            primitive_lengths[i] += other.get_length(i);
        }
        return *this;
    }

    void print() {
        for (float length : primitive_lengths) {
            VTR_LOG("%.1f\t", length);
        }
        VTR_LOG("\n");
    }
};

class ArchModel {
public:


    ArchModel(const AtomNetlist &netlist) : atom_netlist(netlist) {
        // FIXME: The models are not sufficient. We need to distinguish 6-LUT from 5-LUT!
        //          - Same problem for RAMs and DSPs.
        // FIXME: LUTs have an implicit mode called "wire". What do we do in that case?
        t_model* lib_models = g_vpr_ctx.device().arch->model_library;
        t_model* curr_model = lib_models;
        while (curr_model != nullptr) {
            vtr::t_linked_vptr* pb_types = curr_model->pb_types;
            vtr::t_linked_vptr* curr_pb_type = pb_types;
            while (curr_pb_type != nullptr) {
                t_pb_type* pb_type = (t_pb_type *)curr_pb_type->data_vptr;
                VTR_LOG("%s\n", pb_type->name);
                curr_pb_type = curr_pb_type->next;
            }
            blk_model_to_id[curr_model] = id_to_blk_model.size();
            id_to_blk_model.push_back(curr_model);
            curr_model = curr_model->next;
        }
        // Add the models defined in the arch file
        t_model* models = g_vpr_ctx.device().arch->models;
        curr_model = models;
        while (curr_model != nullptr) {
            blk_model_to_id[curr_model] = id_to_blk_model.size();
            id_to_blk_model.push_back(curr_model);
            curr_model = curr_model->next;
        }
        num_unique_atoms = id_to_blk_model.size();
    }

    // FIXME: Need a struct to hold the capacity / size information in a clean way such that we can perform math operations on them.

    void get_mode_permutations_recurse(std::vector<float> permutation, std::vector<std::vector<float>> &permutations, const t_mode &mode, int pb_child_idx, int cummulative_num_pb,
                                std::vector<std::vector<std::vector<float>>> &child_pb_permutations_memoization) {
        VTR_LOG("\t\tRecurse id: %d\n", pb_child_idx);
        // If this is the last child pb, this is a full permutation. Add it.
        if (pb_child_idx == mode.num_pb_type_children) {
            VTR_LOG("\t\tAdding permutation!\n");
            permutations.push_back(std::move(permutation));
            return;
        }
        // Get the permutations for the child_pb of the mode.
        // Note: We have to memoize here since there will be so many permutations,
        //       a recomputing the model permutations will get expensive fast.
        if (child_pb_permutations_memoization.size() <= (size_t)pb_child_idx) {
            child_pb_permutations_memoization.resize(pb_child_idx + 1);
            get_model_permutations(&mode.pb_type_children[pb_child_idx], child_pb_permutations_memoization[pb_child_idx], cummulative_num_pb);
        }
        // For each of the child_pb permutations, accumulate into the partial permutation for this mode.
        for (const auto &child_pb_permutation : child_pb_permutations_memoization[pb_child_idx]) {
            std::vector<float> new_permutation = permutation;
            for (size_t i = 0; i < num_unique_atoms; i++) {
                new_permutation[i] += child_pb_permutation[i];
            }
            // Recurse to construct full permutation from partial.
            get_mode_permutations_recurse(std::move(new_permutation), permutations, mode, pb_child_idx + 1, cummulative_num_pb, child_pb_permutations_memoization);
        }
    }

    // Entry point for getting the mode permutations.
    void get_mode_permutations(const t_mode &mode, std::vector<std::vector<float>> &permutations, int cummulative_num_pb) {
        VTR_LOG("Getting permutations of mode: %s\n", mode.name);
        if (mode.num_pb_type_children == 0) {
            VTR_LOG("\tMode has no children. No permutations for this mode.\n");
            return;
        }
        std::vector<float> mode_permutation_seed(num_unique_atoms, 0.f);
        std::vector<std::vector<std::vector<float>>> child_pb_permutations_memoization;
        get_mode_permutations_recurse(mode_permutation_seed, permutations, mode, 0, cummulative_num_pb, child_pb_permutations_memoization);
    }

    // NOTE: This is very expensive to calculate. Should minimize how often it
    //       is called. Also the returned vector will be very large. Should limit
    //       how much it is stored.
    void get_model_permutations(t_pb_type *pb_type_ptr, std::vector<std::vector<float>> &permutations, int cummulative_num_pb) {
        VTR_LOG("Getting Model Permutations of %s\n", pb_type_ptr->name);
        int num_modes = pb_type_ptr->num_modes;
        int num_instances = cummulative_num_pb * pb_type_ptr->num_pb;
        // If this is a leaf/primitive, only one permutation is possible.
        if (num_modes == 0 && pb_type_ptr->model != nullptr) {
            VTR_LOG("\tModel is a leaf.\n");
            std::vector<float> leaf_capacity(num_unique_atoms, 0.f);
            leaf_capacity[blk_model_to_id[pb_type_ptr->model]] += num_instances;
            permutations.push_back(std::move(leaf_capacity));
            return;
        }
        // If this is not a leaf/primitive, need to explore its modes and get their permutations.
        // Since a pb_type cannot be multiple modes at once, we take the union
        // of their permutations.
        VTR_LOG("\tModel is not a leaf.\n");
        for (int mode_idx = 0; mode_idx < num_modes; mode_idx++) {
            const t_mode& mode = pb_type_ptr->modes[mode_idx];
            get_mode_permutations(mode, permutations, num_instances);
        }
    }

    // Note: There is an implicit simplification here regarding sub-tiles.
    //       A tile may have multiple sub-tiles in the same location, each has their 
    //       own sub-tile index. Here, we are using sub-tile ID, where each id is unique.
    std::vector<std::vector<float>> get_capacity(int tile_x, int tile_y, int tile_layer, int sub_tile_id) {
        // FIXME: How do we take into account the height/width of a tile?
        const auto &device_ctx = g_vpr_ctx.device();
        t_physical_tile_type_ptr tile_type_ptr = device_ctx.grid.get_physical_type({tile_x, tile_y, tile_layer});
        if (tile_type_ptr->is_empty()) {
            std::vector<float> empty_capacity(num_unique_atoms, 0.f);
            return {std::move(empty_capacity)};
        }
        const std::vector<t_sub_tile>& sub_tiles = tile_type_ptr->sub_tiles;
        const t_sub_tile& sub_tile = sub_tiles[sub_tile_id];
        // FIXME: Need to handle sub-tile capacity. The capacity means we can have
        //        N copies of whatever combination of permutations we want.
        // int sub_tile_cap = sub_tile.capacity.total();
        std::vector<std::vector<float>> model_permutations;
        for (t_logical_block_type_ptr log_blk_ty_ptr : sub_tile.equivalent_sites) {
            t_pb_type* pb_type = log_blk_ty_ptr->pb_type;
            get_model_permutations(pb_type, model_permutations, 1);
        }
        
        // Essentially, given an atom block model, how many models can fit in this tile location?
        //
        // device_ctx.grid.get_physical_type
        //
        // t_physical_tile_type: libs/libarchfpga/src/physical_types.h
        //      -> sub_tiles
        // t_sub_tile: libs/libarchfpga/src/physical_types.h
        //      -> equivalent sites
        //      -> capacity?
        // t_logical_block_type: libs/libarchfpga/src/physical_types.h
        //      -> pb_type
        // t_pb_type: libs/libarchfpga/src/physical_types.h
        //      -> MODEL!!!!! (but this may not be the actual model we want...
        //      -> This is where the mode is found.
        
        // Something interesting is how to handle sub_tiles and equivalent sites.
        // It sounds like sub_tiles are just different blocks on a single tile.
        // It sounds like "equivalent sites" are the mode. (to confirm with VB)
        // Not sure how to hanlde the equivalent sites. It could make this very complicated.
        // For now the architecture we are using to test does not have equivalent sites.
        
        // Something to confirm: For the atom netlist, what models can be found?
        //  - Do they already have a mode?
        //  - Are they guarenteed to be at the base of the heirarchy?
        //  - How do equivalent sites fit into this?

        // Do need to handle modes somehow... Maybe leave this up to the legalizer?
        // - we could pass back a vector of possible vectors
        // - We could just take the max and call it a day, no need to be that accurate for a partial legalizer!

        // TODO: I wonder if equivalent sites can be baked into the abstraction somehow.
        return model_permutations;
    }

    std::vector<float> get_size(t_pack_molecule* mol) {
        std::vector<float> size(num_unique_atoms, 0.f);
        for (const AtomBlockId &blk_id : mol->atom_block_ids) {
            const t_model* blk_model = atom_netlist.block_model(blk_id);
            size_t id = blk_model_to_id[blk_model];
            size[id] += 1;
        }
        return size;
    }

    std::vector<const t_model*> id_to_blk_model;
private:
    const AtomNetlist& atom_netlist;
    size_t num_unique_atoms;
    std::unordered_map<const t_model*, size_t> blk_model_to_id;
};

// TODO: This should probably be moved into its own class for ArchModel.
struct Bin {
    // The IDs of the nodes contained within this bin.
    std::unordered_set<size_t> contained_node_ids;
    // Is this bin blocked for any reason.
    // TODO: From the perspective of the atom being moved, if it cannot be
    //       moved into this block it should be considered as blocked. So this
    //       value may need to be computed at runtime.
    bool is_blocked = false;

    // For now store the bin locations in the struct.
    // FIXME: This can be optimized.
    double x_pos;
    double y_pos;
};

} // namespace

static inline size_t getBinId(size_t x, size_t y) {
    size_t grid_width = g_vpr_ctx.device().grid.width();
    return x + (y * grid_width);
}

static inline void getNeighborsOfBin(const Bin& bin, std::vector<size_t>& neighbors) {
    size_t grid_width = g_vpr_ctx.device().grid.width();
    size_t grid_height = g_vpr_ctx.device().grid.width();
    if (bin.x_pos >= 1)
        neighbors.push_back(getBinId(bin.x_pos - 1, bin.y_pos));
    if (bin.x_pos <= grid_width - 2)
        neighbors.push_back(getBinId(bin.x_pos + 1, bin.y_pos));
    if (bin.y_pos >= 1)
        neighbors.push_back(getBinId(bin.x_pos, bin.y_pos - 1));
    if (bin.y_pos <= grid_height - 2)
        neighbors.push_back(getBinId(bin.x_pos, bin.y_pos + 1));
}

// Helper method to compute the phi term in the durav algorithm.
static inline double computeMaxMovement(size_t iter) {
    return 100 * (iter + 1) * (iter + 1);
}

// FIXME: REMOVE THIS!
#define BIN_CAPACITY 10

// Helper method to get the supply of the bin.
//  returns 0 if the bin is empty
//  else returns node-size - 1
// FIXME: This implies that each bin can only hold one thing.
static inline size_t getSupply(const Bin& b) {
    if (b.contained_node_ids.empty())
        return 0;
    size_t bin_capacity = BIN_CAPACITY;
    if (b.contained_node_ids.size() <= bin_capacity)
        return 0;
    return b.contained_node_ids.size() - bin_capacity;
}

// Helper method to get the closest node in the src bin to the sink bin.
// NOTE: This is from the original solved position of the node in the src bin.
static inline size_t getClosestNode(const PartialPlacement& p_placement,
                                    std::vector<Bin>& bins,
                                    size_t src_bin_id,
                                    size_t sink_bin_id,
                                    double &smallest_block_dist) {
    double sink_bin_x = bins[sink_bin_id].x_pos;
    double sink_bin_y = bins[sink_bin_id].y_pos;

    const Bin& src_bin = bins[src_bin_id];
    VTR_ASSERT(!src_bin.contained_node_ids.empty()
                && "Cannot get closest block of empty bin.");
    size_t closest_node_id = 0;
    smallest_block_dist = std::numeric_limits<double>::infinity();
    for (size_t node_id : src_bin.contained_node_ids) {
        // NOTE: Slight change from original implementation. Not getting bin pos.
        double orig_node_bin_x = p_placement.node_loc_x[node_id];
        double orig_node_bin_y = p_placement.node_loc_y[node_id];
        double dx = orig_node_bin_x - sink_bin_x;
        double dy = orig_node_bin_y - sink_bin_y;
        double block_dist = (dx * dx) + (dy * dy);
        if (block_dist < smallest_block_dist) {
            closest_node_id = node_id;
            smallest_block_dist = block_dist;
        }
    }
    return closest_node_id;
}

// Helper method to compute the cost of moving objects from the src bin to the
// sink bin.
static inline double computeCost(const PartialPlacement& p_placement,
                                 std::vector<Bin>& bins,
                                 double psi,
                                 size_t src_bin_id,
                                 size_t sink_bin_id) {
    // If the src bin is empty, then there is nothing to move.
    if (bins[src_bin_id].contained_node_ids.empty())
        return std::numeric_limits<double>::infinity();

    // Get the weight, which is proportional to the size of the set that
    // contains the src.
    // FIXME: What happens when this is zero?
    size_t src_supply = getSupply(bins[src_bin_id]);
    double wt = static_cast<double>(src_supply);

    // Get the minimum quadratic movement to move a cell from the src bin to
    // the sink bin.
    // NOTE: This assumes no diagonal movements.
    double quad_mvmt;
    getClosestNode(p_placement, bins, src_bin_id, sink_bin_id, quad_mvmt);

    // If the movement is larger than psi, return infinity
    if (quad_mvmt >= psi)
        return std::numeric_limits<double>::infinity();

    // Following the equation from the Darav paper.
    return wt * quad_mvmt;
}

// Helper method to get the paths to flow nodes between bins.
static inline void getPaths(const PartialPlacement& p_placement,
                            std::vector<Bin>& bins,
                            double psi,
                            size_t starting_bin_id,
                            std::vector<std::vector<size_t>>& paths) {
    // Create a visited vector.
    std::vector<bool> bin_visited(bins.size());
    std::fill(bin_visited.begin(), bin_visited.end(), false);
    bin_visited[starting_bin_id] = true;
    // Create a cost array. A path can be uniquely identified by its tail bin
    // cost.
    std::vector<double> bin_cost(bins.size());
    std::fill(bin_cost.begin(), bin_cost.end(), 0.0);
    // Create a starting path.
    std::vector<size_t> starting_path;
    starting_path.push_back(starting_bin_id);
    // Create a FIFO queue.
    std::queue<std::vector<size_t>> queue;
    queue.push(std::move(starting_path));

    size_t demand = 0;
    size_t starting_bin_supply = getSupply(bins[starting_bin_id]);
    while (!queue.empty() && demand < starting_bin_supply) {
        std::vector<size_t> &p = queue.front();
        size_t tail_bin_id = p.back();
        // Get the neighbors of the tail bin.
        std::vector<size_t> neighbor_bin_ids;
        getNeighborsOfBin(bins[tail_bin_id], neighbor_bin_ids);
        for (size_t neighbor_bin_id : neighbor_bin_ids) {
            // TODO: This is where the bin being blocked can be calculated.
            // FIXME: I do not think this should be here...
            //        This may be encapsulated in the neighbor connections.
            if (bins[neighbor_bin_id].is_blocked)
                continue;
            if (bin_visited[neighbor_bin_id])
                continue;
            double cost = computeCost(p_placement, bins, psi, tail_bin_id, neighbor_bin_id);
            if (cost < std::numeric_limits<double>::infinity()) {
                std::vector<size_t> p_copy(p);
                bin_cost[neighbor_bin_id] = bin_cost[tail_bin_id] + cost;
                p_copy.push_back(neighbor_bin_id);
                bin_visited[neighbor_bin_id] = true;
                // NOTE: Needed to change this from the algorithm since it was
                //       skipping partially filled bins. Maybe indicative of a
                //       bug somewhere.
                // if (bins[neighbor_bin_id].contained_node_ids.empty()) {
                if (bins[neighbor_bin_id].contained_node_ids.size() < BIN_CAPACITY) {
                    paths.push_back(std::move(p_copy));
                    demand += BIN_CAPACITY - bins[neighbor_bin_id].contained_node_ids.size();
                } else {
                    queue.push(std::move(p_copy));
                }
            }
        }
        queue.pop();
    }

    // Sort the paths in increasing order of cost.
    std::sort(paths.begin(), paths.end(), [&](std::vector<size_t>& a, std::vector<size_t>& b) {
        return bin_cost[a.back()] < bin_cost[b.back()];
    });
}

// Helper method to move cells along a given path.
static inline void moveCells(const PartialPlacement& p_placement,
                             std::vector<Bin>& bins,
                             double psi,
                             std::vector<size_t>& path) {
    VTR_ASSERT(!path.empty());
    size_t src_bin_id = path[0];
    std::stack<size_t> s;
    s.push(src_bin_id);
    size_t path_size = path.size();
    for (size_t path_index = 1; path_index < path_size; path_index++) {
        size_t sink_bin_id = path[path_index];
        double cost = computeCost(p_placement, bins, psi, src_bin_id, sink_bin_id);
        if (cost == std::numeric_limits<double>::infinity())
            return;
        //  continue;
        src_bin_id = sink_bin_id;
        s.push(sink_bin_id);
    }
    size_t sink_bin_id = s.top();
    s.pop();
    while (!s.empty()) {
        src_bin_id = s.top();
        s.pop();
        // Minor change to the algorithm proposed by Darav et al., find the
        // closest point in src to sink and move it to sink (instead of sorting
        // the whole list which is wasteful).
        // TODO: Verify this.
        Bin& src_bin = bins[src_bin_id];
        VTR_ASSERT(!src_bin.contained_node_ids.empty() && "Cannot move a node from an empty bin.");
        double smallest_block_dist;
        size_t closest_node_id = getClosestNode(p_placement, bins, src_bin_id, sink_bin_id, smallest_block_dist);
        src_bin.contained_node_ids.erase(closest_node_id);
        Bin& sink_bin = bins[sink_bin_id];
        sink_bin.contained_node_ids.insert(closest_node_id);

        sink_bin_id = src_bin_id;
    }
}

// Flow-Based Legalizer based off work by Darav et al.
// https://dl.acm.org/doi/10.1145/3289602.3293896
// FIXME: This currently is not working since it treats LUTs and FFs as the same
//        block.
void FlowBasedLegalizer::legalize(PartialPlacement &p_placement) {
    VTR_LOG("Running Flow-Based Legalizer\n");

    ArchModel arch_model(p_placement.atom_netlist);

    const std::vector<std::vector<float>> sub_tile_capacity = arch_model.get_capacity(2,1,0,0);
    for (const t_model* model : arch_model.id_to_blk_model) {
        VTR_LOG("%s\t", model->name);
    }
    VTR_LOG("\n");
    for (const auto &capacity : sub_tile_capacity) {
        for (float cap : capacity) {
            VTR_LOG("%f\t", cap);
        }
        VTR_LOG("\n");
    }

    VTR_ASSERT(false);

    // Initialize the bins with the node positions.
    // TODO: Encapsulate this into a class.
    size_t grid_width = g_vpr_ctx.device().grid.width();
    size_t grid_height = g_vpr_ctx.device().grid.height();
    VTR_LOG("\tGrid size = (%zu, %zu)\n", grid_width, grid_height);
    std::vector<Bin> bins(grid_width * grid_height);
    for (size_t i = 0; i < grid_width; i++) {
        for (size_t j = 0; j < grid_height; j++) {
            size_t bin_id = getBinId(i, j);
            bins[bin_id].x_pos = i;
            bins[bin_id].y_pos = j;
        }
    }
    // FIXME: What about the fixed blocks? Shouldnt they be put into the bins?
    for (size_t node_id = 0; node_id < p_placement.num_moveable_nodes; node_id++) {
        double node_x = p_placement.node_loc_x[node_id];
        double node_y = p_placement.node_loc_y[node_id];
        // FIXME: This conversion is not necessarily correct. The bin position
        //        should be the center of the tile, not the corner.
        size_t bin_id = getBinId(node_x, node_y);
        bins[bin_id].contained_node_ids.insert(node_id);
    }
    // Label the IO tiles as obstacles.
    // FIXME: This should be dynamically calculated, as explained above.
    // Do not know how to do this properly. For now assume the IO is on the
    // perimeter without the corners.
    for (size_t i = 0; i < grid_width; i++) {
        bins[getBinId(i, 0)].is_blocked = true;
        bins[getBinId(i, grid_height - 1)].is_blocked = true;
    }
    for (size_t j = 0; j < grid_height; j++) {
        bins[getBinId(0, j)].is_blocked = true;
        bins[getBinId(grid_width - 1, j)].is_blocked = true;
    }

    // Run the flow-based spreader.
    size_t flowBasedIter = 0;
    while (true) {
        if (flowBasedIter > 1000) {
            VTR_LOG("HIT MAX ITERATION!!!\n");
            break;
        }
        // Compute the max movement.
        double psi = computeMaxMovement(flowBasedIter);

        // Get the overfilled bins and sort them in increasing order of supply.
        std::vector<size_t> overfilled_bin_ids;
        size_t num_bins = bins.size();
        for (size_t bin_id = 0; bin_id < num_bins; bin_id++) {
            if (getSupply(bins[bin_id]) > 0)
                overfilled_bin_ids.push_back(bin_id);
        }
        if (overfilled_bin_ids.empty()) {
            // VTR_LOG("No overfilled bins! No spreading needed!\n");
            break;
        }
        std::sort(overfilled_bin_ids.begin(), overfilled_bin_ids.end(), [&](size_t a, size_t b) {
            return getSupply(bins[a]) < getSupply(bins[b]);
        });
        // VTR_LOG("Num overfilled bins: %zu\n", overfilled_bin_ids.size());
        // VTR_LOG("\tLargest overfilled bin supply: %zu\n", getSupply(bins[overfilled_bin_ids.back()]));
        // VTR_LOG("\tpsi = %f\n", psi);

        for (size_t bin_id : overfilled_bin_ids) {
            // Get the list of candidate paths based on psi.
            //  NOTE: The paths are sorted by increasing cost within the
            //        getPaths method.
            std::vector<std::vector<size_t>> paths;
            getPaths(p_placement, bins, psi, bin_id, paths);

            // VTR_LOG("\tNum paths: %zu\n", paths.size());
            for (std::vector<size_t>& path : paths) {
                if (getSupply(bins[bin_id]) == 0)
                    continue;
                // Move cells over the paths.
                //  NOTE: This will modify the bins. (actual block positions
                //        will not change (yet)).
                moveCells(p_placement, bins, psi, path);
            }
        }

        // TODO: Get the total cell displacement for debugging.

        flowBasedIter++;
    }
    VTR_LOG("Flow-Based Legalizer finished in %zu iterations.\n", flowBasedIter + 1);

    // Update the partial placement with the spread blocks.
    // TODO: Look into only updating the nodes that actually moved.
    for (const Bin& bin : bins) {
        for (size_t node_id : bin.contained_node_ids) {
            // TODO: This movement is a bit strange. It should probably move to
            //       the center of the block or the edge closest to the solved
            //       position.
            double node_x = p_placement.node_loc_x[node_id];
            double node_y = p_placement.node_loc_y[node_id];
            double offset_x = node_x - std::floor(node_x);
            double offset_y = node_y - std::floor(node_y);
            p_placement.node_loc_x[node_id] = bin.x_pos + offset_x;
            p_placement.node_loc_y[node_id] = bin.y_pos + offset_y;
        }
    }
}

/*
// Helper function that removes the placement information of a cluster located
// at a given sub_tile.
// This was taken from place/cut_spreader.cpp
// NOTE: This may be useless since we can use place/place_util.cpp:load_grid_blocks_from_block_locs
static inline void unbind_subtile(t_pl_loc sub_tile) {
    PlacementContext& place_ctx = g_vpr_ctx.mutable_placement();
    ClusterBlockId blk = place_ctx.grid_blocks.block_at_location(sub_tile);
    VTR_ASSERT(blk != EMPTY_BLOCK_ID);
    VTR_ASSERT(place_ctx.block_locs[blk].is_fixed == false);
    // Clear the block locs.
    place_ctx.block_locs[blk].loc = t_pl_loc{};
    // Clear the grid blocks.
    place_ctx.grid_blocks.set_block_at_location(sub_tile, EMPTY_BLOCK_ID);
    place_ctx.grid_blocks.set_usage({sub_tile.x, sub_tile.y, sub_tile.layer}, place_ctx.grid_blocks.get_usage({sub_tile.x, sub_tile.y, sub_tile.layer}) - 1);
}

// Helper function that binds the given cluster to a subtile.
// This was taken from place/cut_spreader.cpp
// NOTE: This is very similar to place/place_util.cpp:set_block_location...
static inline void bind_subtile(ClusterBlockId blk, t_pl_loc sub_tile) {
    PlacementContext& place_ctx = g_vpr_ctx.mutable_placement();
    VTR_ASSERT(place_ctx.grid_blocks.block_at_location(sub_tile) == EMPTY_BLOCK_ID);
    VTR_ASSERT(place_ctx.block_locs[blk].is_fixed == false);
    // Set the block locs.
    place_ctx.block_locs[blk].loc = sub_tile;
    // Set the grid blocks.
    place_ctx.grid_blocks.set_block_at_location(sub_tile, blk);
    place_ctx.grid_blocks.set_usage({sub_tile.x, sub_tile.y, sub_tile.layer}, place_ctx.grid_blocks.get_usage({sub_tile.x, sub_tile.y, sub_tile.layer}) + 1);
}
*/

// TODO: Move this into its own file.
void FullLegalizer::legalize(PartialPlacement& p_placement) {
    (void)p_placement;
    VTR_LOG("Running Full Legalizer\n");
    // Since we are currently assuming that clustering has already occured,
    // this legalizer should use the reclustering API to handle creating legal
    // clusters.
    // NOTE: Pre-packing must have been performed before this.

    // NOTE: The reclustering API changes who goes in which ClusterBlockId.
    //       The placement of the ClusterBlockId is controlled by the placer.
    // Goal: Make sure each cluster has a ClusterBlockId and a placement.
    //  - Re-clustering API can create new ClusterBlockIds
    //  - I do not think it can remove a ClusterBlockId (double check)
    //      - This may be an issue if the number of clusters is reduced by us.
    // From inspection, I do not think clustering_data is used when during
    // packing is set to false.

    // Strategy:
    //  - For each tile that has nodes within it:
    //      - Get the first block's ClusterBlockId (this will be the cluster
    //        for this tile.
    //      - Swap / move all other molecules out of this block that do not
    //        belong; while getting the molecules that do belong in.
    //          - If there is nowhere to put a block, create an orphan cluster
    //      - Set the location of this cluster in the placement.
    
    // FIXME: This was stolen from above. Make common. Perhaps make its own custom data structure since the bins may not match the clusters.
    // Initialize the bins with the node positions.
    // FIXME: THIS SHOULD USE GRID BLOCKS FROM THE PLACEMENT CONTEXT.
    // TODO: Encapsulate this into a class.
    const DeviceContext& device_ctx = g_vpr_ctx.device();
    size_t grid_width = device_ctx.grid.width();
    size_t grid_height = device_ctx.grid.height();
    VTR_LOG("\tGrid size = (%zu, %zu)\n", grid_width, grid_height);
    std::vector<Bin> bins(grid_width * grid_height);
    for (size_t i = 0; i < grid_width; i++) {
        for (size_t j = 0; j < grid_height; j++) {
            size_t bin_id = getBinId(i, j);
            bins[bin_id].x_pos = i;
            bins[bin_id].y_pos = j;
        }
    }
    // FIXME: What about the fixed blocks? Shouldnt they be put into the bins?
    for (size_t node_id = 0; node_id < p_placement.num_moveable_nodes; node_id++) {
        double node_x = p_placement.node_loc_x[node_id];
        double node_y = p_placement.node_loc_y[node_id];
        // FIXME: This conversion is not necessarily correct. The bin position
        //        should be the center of the tile, not the corner.
        size_t bin_id = getBinId(node_x, node_y);
        bins[bin_id].contained_node_ids.insert(node_id);
    }

    // See place/place_util.cpp:alloc_and_load_legal_placement_locations
    // Shows how to iterate over the device grid to create locations
    
    // A list of legal locations
    // (tile_type_index, sub_tile_index) -> legal locations at that index
    // NOTE: This is tile TYPE not tile index.
    // TODO: look into equivalent tiles / sites
    //  see: place/cut_spreader.cpp:CutSpreader
    std::vector<std::vector<std::vector<t_pl_loc>>> legal_locs;
    alloc_and_load_legal_placement_locations(legal_locs);

    // A lookup table to get the valid locations at a given tile.
    // (tile_pos_x, tile_pos_y) -> legal locations at that tile
    // FIXME: This is a terrible and inneficient process. Should be fixed.
    //        It is also not very generalizable. Assumes that each tile has
    //        a single sub-tile.
    //  - Should replicate the code in alloc_and_load_legal_placement_locations
    std::vector<std::vector<std::vector<t_pl_loc>>> tile_locs;
    tile_locs.resize(grid_width);
    for (size_t i = 0; i < grid_width; i++)
        tile_locs[i].resize(grid_height);
    for (const auto& tile_type_locs : legal_locs) {
        for (const auto& sub_tile_locs : tile_type_locs) {
            for (const t_pl_loc& loc : sub_tile_locs) {
                tile_locs[loc.x][loc.y].push_back(loc);
            }
        }
    }

    // Initialize the placement context.
    // This resets the locations of the clusters and the grid blocks (reverse lookup).
    // Use set_block_location to set the location fo a cluster.
    init_placement_context();

    // For each bin with moveable blocks in it, assign a cluster to it.
    // FIXME: This will not generalize well.
    std::vector<ClusterBlockId> bin_cluster(bins.size(), ClusterBlockId::INVALID());
    std::vector<ClusterBlockId> moveable_clusters;

    const ClusteringContext& cluster_ctx = g_vpr_ctx.clustering();
    const AtomContext& atom_ctx = g_vpr_ctx.atom();
    ClusterAtomsLookup clusterBlockToAtomBlockLookup;
    clusterBlockToAtomBlockLookup.init_lookup();
    for (ClusterBlockId blk_id : cluster_ctx.clb_nlist.blocks()) {
        // Crude assumption, if the cluster contains any IO blocks, it cannot be used.
        bool has_io_blocks = false;
        for (AtomBlockId atom_blk_id : clusterBlockToAtomBlockLookup.atoms_in_cluster(blk_id)) {
            AtomBlockType atom_blk_ty = atom_ctx.nlist.block_type(atom_blk_id);
            if (atom_blk_ty == AtomBlockType::INPAD || atom_blk_ty == AtomBlockType::OUTPAD)
                has_io_blocks = true;
        }
        if (has_io_blocks)
            continue;
        moveable_clusters.push_back(blk_id);
    }

    VTR_ASSERT(device_ctx.grid.get_num_layers() == 1 && "3D FPGA not supported for AP");
    for (size_t i = 0; i < grid_width; i++) {
        for (size_t j = 0; j < grid_height; j++) {
            size_t bin_id = getBinId(i, j);
            const Bin& bin = bins[bin_id];
            // Does this bin have any moveable nodes?
            std::vector<size_t> moveable_nodes;
            moveable_nodes.reserve(bin.contained_node_ids.size());
            for (size_t node_id : bin.contained_node_ids) {
                if (p_placement.is_moveable_node(node_id))
                    moveable_nodes.push_back(node_id);
            }
            if (moveable_nodes.empty())
                continue;
            
            VTR_ASSERT(moveable_clusters.size() > 0);
            ClusterBlockId new_cluster_id = moveable_clusters.back();
            bin_cluster[bin_id] = new_cluster_id;
            moveable_clusters.pop_back();

            VTR_ASSERT(tile_locs[i][j].size() == 1 && "Only support single-location tiles for now");
            set_block_location(new_cluster_id, tile_locs[i][j][0]);
        }
    }
    VTR_LOG("NUM LEFT: %zu\n", moveable_clusters.size());
    // FIXME: This is awful and hard coded.
    std::vector<size_t> empty_bins;
    for (size_t i = 1; i < grid_width - 1; i++) {
        for (size_t j = 1; j < grid_width - 1; j++) {
            size_t bin_id = getBinId(i, j);
            if (bins[bin_id].contained_node_ids.empty())
                continue;
            empty_bins.push_back(bin_id);
        }
    }
    VTR_ASSERT(empty_bins.size() >= moveable_clusters.size());
    for (size_t i = 0; i < moveable_clusters.size(); i++) {
        // Just shove them in some random empty block
        size_t bin_id = empty_bins.back();
        empty_bins.pop_back();
        bin_cluster[bin_id] = moveable_clusters[i];

        size_t x = bins[bin_id].x_pos;
        size_t y = bins[bin_id].y_pos;
        VTR_ASSERT(tile_locs[x][y].size() == 1 && "Only support single-location tiles for now");
        set_block_location(moveable_clusters[i], tile_locs[x][y][0]);
    }

    // reverse lookup for cluster to bin_id
    // std::map<ClusterBlockId, size_t> cluster_to_bin_id;
    // for (size_t i = 0; i < bin_cluster.size(); i++) {
    //     cluster_to_bin_id[bin_cluster[i]] = i;
    // }

    // Recluster.
    size_t num_immovable = 0;
    for (size_t i = 0; i < p_placement.num_moveable_nodes; i++) {
        VTR_LOG("%zu\n", i);
        // Which cluster do I want to be in?
        size_t target_bin_id = getBinId(p_placement.node_loc_x[i], p_placement.node_loc_y[i]);
        ClusterBlockId target_cluster = bin_cluster[target_bin_id];
        VTR_ASSERT(target_cluster != ClusterBlockId::INVALID());
        // Which cluster am I in?
        t_pack_molecule* mol = p_placement.node_id_to_mol[i];
        ClusterBlockId src_cluster = atom_ctx.lookup.atom_clb(mol->atom_block_ids[0]);
        if (src_cluster == target_cluster)
            continue;
        VTR_LOG("I GOT HERE1\n");
        // Try to move into the target cluster
        t_clustering_data temp_data;
        bool success = false;
        success = move_mol_to_existing_cluster(mol, target_cluster, false, 3, temp_data);
        if (success)
            continue;
        VTR_LOG("I GOT HERE2\n");
        // If that doesnt work try swapping with a block in the target cluster
        const Bin& target_bin = bins[target_bin_id];
        std::vector<t_pack_molecule*> potential_partners;
        for (size_t node_id : target_bin.contained_node_ids) {
            t_pack_molecule* partner_mol = p_placement.node_id_to_mol[node_id];
            ClusterBlockId partner_src = atom_ctx.lookup.atom_clb(partner_mol->atom_block_ids[0]);
            if (partner_src != target_cluster)
                potential_partners.push_back(partner_mol);
        }
        for (t_pack_molecule* partner_mol : potential_partners) {
            success = swap_two_molecules(mol, partner_mol, false, 3, temp_data);
            if (success)
                continue;
        }
        if (!success)
            num_immovable++;
    }
    VTR_LOG("NUM IMMOVABLE: %zu\n", num_immovable);

    // for (size_t i = 0; i < grid_width; i++) {
    //     for (size_t j = 0; j < grid_height; j++) {
    //         size_t bin_id = getBinId(i, j);
    //         const Bin& bin = bins[bin_id];
    //         // Does this bin have any moveable nodes?
    //         std::vector<size_t> moveable_nodes;
    //         moveable_nodes.reserve(bin.contained_node_ids.size());
    //         for (size_t node_id : bin.contained_node_ids) {
    //             if (p_placement.is_moveable_node(node_id))
    //                 moveable_nodes.push_back(node_id);
    //         }
    //         if (moveable_nodes.empty())
    //             continue;
    //         VTR_ASSERT(tile_locs[i][j].size() == 1 && "Only support single-location tiles for now");
    //         VTR_LOG("%zu ", tile_locs[i][j].size());
    //         // t_physical_tile_type_ptr tile = device_ctx.grid.get_physical_type({(int)i, (int)j, 0});
    //         // VTR_ASSERT(tile->sub_tiles.size() == 1 && "Right now assuming that each tile only contains one sub-tile");
    //         // t_sub_tile sub_tile = tile->sub_tiles[0];
    //         // size_t num = 0;
    //         // for (auto loc : legal_locs[tile->index][sub_tile.index])
    //         //     num++;
    //         // VTR_LOG("%zu ", num);
    //         // // VTR_LOG("%zu ", legal_locs[tile->index][sub_tile.index].size());

    //         // For each of the moveable nodes, get them into the same ClusterBlockId
    //         size_t first_node_id = moveable_nodes[0];
    //         t_pack_molecule* first_mol = p_placement.node_id_to_mol[first_node_id];
    //         ClusterBlockId new_cluster_id = atom_ctx.lookup.atom_clb(first_mol->atom_block_ids[0]);
    //     }
    //     VTR_LOG("\n");
    // }

    // PlacementContext& place_ctx = g_vpr_ctx.mutable_placement();
    // for (size_t i = 0; i < grid_width; i++) {
    //     for (size_t j = 0; j < grid_height; j++) {
    //         size_t bin_id = getBinId(i, j);
    //         const Bin& bin = bins[bin_id];
    //         // Does this bin have any moveable nodes?
    //         std::vector<size_t> moveable_nodes;
    //         moveable_nodes.reserve(bin.contained_node_ids.size());
    //         for (size_t node_id : bin.contained_node_ids) {
    //             if (p_placement.is_moveable_node(node_id))
    //                 moveable_nodes.push_back(node_id);
    //         }
    //         if (moveable_nodes.empty())
    //             continue;
    //         // For each of the moveable nodes, get them into the same ClusterBlockId
    //         size_t first_node_id = moveable_nodes[0];
    //         t_pack_molecule* first_mol = p_placement.node_id_to_mol[first_node_id];
    //         ClusterBlockId new_cluster_id = atom_ctx.lookup.atom_clb(first_mol->atom_block_ids[0]);
    //         // Set the position of this ClusterBlockId to be this tile.
    //     }
    // }

    // for (size_t i = 0; i < p_placement.num_moveable_nodes; i++) {
    //     // VTR_LOG("Moveable Node ID: %zu\n", i);
    //     t_pack_molecule* mol = p_placement.node_id_to_mol[i];
    //     VTR_ASSERT(mol->atom_block_ids.size() == 1 && "Molecule expected to have one Atom.");
    // }

    // // There are too many issues with using the reclustering API, namely:
    // //  1) Requires packing to run first, which we hope not to do.
    // //  2) Requires the use of molecules, which we will likely not use the same
    // //     molecules as the packer.
    // // Instead, I will recreate the Clustered netlist from scratch. To make this
    // // easier, I plan to replicate the process done when loading a .net file.
    // ClusteringContext& cluster_ctx = g_vpr_ctx.mutable_clustering();

    // // TODO: Maybe put this in a method to prevent duplicate code.
    // // See /base/vpr_api.cpp:vpr_load_packing
    // /* Ensure we have a clean start with void net remapping information */
    // cluster_ctx.post_routing_clb_pin_nets.clear();
    // cluster_ctx.pre_routing_net_pin_mapping.clear();

    // // See /base/vpr_api.cpp:free_circuit
    // // Free new net structures
    // for (ClusterBlockId blk_id : cluster_ctx.clb_nlist.blocks())
    //     cluster_ctx.clb_nlist.remove_block(blk_id);
    // cluster_ctx.clb_nlist = ClusteredNetlist();

    // // See /base/read_netlist.cpp:read_netlist
    // // NOTE: Not reloading the atom/pb mapping. Assuming that it is accurate.
    // // See /base/read_netlist.cpp:processComplexBlock
    // const DeviceContext& device_ctx = g_vpr_ctx.device();
    // const DeviceGrid& device_grid = device_ctx.grid;
    // size_t device_width = device_grid.width();
    // size_t device_height = device_grid.height();
    // for (size_t i = 0; i < device_width; i++) {
    //     for (size_t j = 0; j < device_height; j++) {
    //         std::vector<size_t> contained_nodes;
    //         // FIXME: This is very inneficient. Precompute which tiles have
    //         //        which block.
    //         for (size_t node_id = 0; node_id < p_placement.num_nodes; node_id++) {
    //             size_t tile_x = static_cast<size_t>(p_placement.node_loc_x[node_id]);
    //             size_t tile_y = static_cast<size_t>(p_placement.node_loc_y[node_id]);
    //             if (tile_x == i && tile_y == j)
    //                 contained_nodes.push_back(node_id);
    //         }
    //         VTR_LOG("%zu ", contained_nodes.size());
    //     }
    //     VTR_LOG("\n");
    // }
}

