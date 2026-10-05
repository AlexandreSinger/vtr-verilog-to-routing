/**
 * @file
 * @author  Alex Singer
 * @date    October 2026
 * @brief   Implementation of methods for estimating the delays of timing arcs
 *          using a flat placement.
 */

#include "ap_timing_estimation.h"
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <utility>
#include "PreClusterDelayCalculator.h"
#include "PreClusterTimingManager.h"
#include "ap_netlist.h"
#include "device_grid.h"
#include "partial_placement.h"
#include "pb_type_graph.h"
#include "physical_types.h"
#include "place_delay_model.h"
#include "router_lookahead_constants.h"
#include "timing_info.h"
#include "vtr_assert.h"
#include "vtr_hash.h"

/**
 * @brief Queries the delay model at a single reference tile pair (grid center
 *        to center+1) to get a representative delay-per-tile estimate.
 *
 * Used as a fallback when the delay model returns ROUTER_LOOKAHEAD_NO_PATH_SENTINEL
 * for a driver/sink pair. Returns 0.0f if the reference point is also missing.
 *
 * TODO: It is possible that the tile at the center of the device has no possible
 *       routes within one tile unit. We should have a more systematic way of
 *       doing this. For now this is better than just returning 0.0.
 */
static float get_delay_per_tile(const PlaceDelayModel& place_delay_model,
                                const DeviceGrid& device_grid) {
    int cx = (int)device_grid.width() / 2;
    int cy = (int)device_grid.height() / 2;
    t_physical_tile_loc from_loc(cx, cy, 0);
    t_physical_tile_loc to_loc(cx + 1, cy, 0);
    float d = place_delay_model.delay(from_loc, 0, to_loc, 0);
    if (d >= ROUTER_LOOKAHEAD_NO_PATH_SENTINEL)
        return 0.0f;
    return d;
}

/**
 * @brief Get the root location of the tile which contains the given block in
 *        the given flat placement.
 *
 * The root location uniquely identifies a tile, even for tiles larger than 1x1.
 */
static t_physical_tile_loc get_containing_tile_root_loc(APBlockId blk_id,
                                                        const PartialPlacement& p_placement,
                                                        const DeviceGrid& device_grid) {
    t_physical_tile_loc tile_loc = p_placement.get_containing_tile_loc(blk_id);
    VTR_ASSERT_SAFE(device_grid.is_valid_tile_loc(tile_loc));
    return device_grid.get_root_location(tile_loc);
}

FlatPlacementArcDelayEstimator::FlatPlacementArcDelayEstimator(const APNetlist& ap_netlist,
                                                               const PreClusterDelayCalculator& delay_calc,
                                                               const PlaceDelayModel& place_delay_model,
                                                               const DeviceGrid& device_grid)
    : ap_netlist_(ap_netlist)
    , place_delay_model_(place_delay_model)
    , device_grid_(device_grid) {
    delay_per_tile_ = get_delay_per_tile(place_delay_model_, device_grid_);

    // Pre-compute the intra-cluster delays of every AP pin. These only depend
    // on the pb_graph pins that the primitive pins are expected to be
    // implemented by, not on the placement.
    // Many primitive pins share the same pb_graph pin, so cache the searches
    // through the pb_graph to save time.
    std::unordered_map<const t_pb_graph_pin*, float> to_boundary_cache;
    std::unordered_map<const t_pb_graph_pin*, float> from_boundary_cache;
    std::unordered_map<std::pair<const t_pb_graph_pin*, const t_pb_graph_pin*>, float, vtr::hash_pair> path_cache;

    pin_cluster_boundary_delay_.resize(ap_netlist_.pins().size(), 0.0f);
    sink_pin_intra_cluster_delay_.resize(ap_netlist_.pins().size(), -1.0f);
    for (APPinId pin_id : ap_netlist_.pins()) {
        const t_pb_graph_pin* gpin = delay_calc.find_pb_graph_pin(ap_netlist_.pin_atom_pin(pin_id));

        // Compute the delay between the pin and the boundary of its cluster.
        // If no path to the boundary is found, assume that the delay is small.
        float boundary_delay;
        if (ap_netlist_.pin_type(pin_id) == PinType::DRIVER) {
            auto [it, inserted] = to_boundary_cache.try_emplace(gpin, 0.0f);
            if (inserted)
                it->second = calc_pb_graph_delay_to_root_pin(gpin);
            boundary_delay = it->second;
        } else {
            auto [it, inserted] = from_boundary_cache.try_emplace(gpin, 0.0f);
            if (inserted)
                it->second = calc_pb_graph_delay_from_root_pin(gpin);
            boundary_delay = it->second;
        }
        pin_cluster_boundary_delay_[pin_id] = std::max(boundary_delay, 0.0f);

        // For sink pins, compute the intra-cluster delay from the driver of
        // the net, in case the driver and sink end up in the same cluster.
        if (ap_netlist_.pin_type(pin_id) != PinType::SINK)
            continue;
        APPinId driver_pin_id = ap_netlist_.net_driver(ap_netlist_.pin_net(pin_id));
        if (!driver_pin_id.is_valid())
            continue;
        const t_pb_graph_pin* driver_gpin = delay_calc.find_pb_graph_pin(ap_netlist_.pin_atom_pin(driver_pin_id));
        auto [it, inserted] = path_cache.try_emplace({driver_gpin, gpin}, -1.0f);
        if (inserted)
            it->second = calc_pb_graph_path_delay(driver_gpin, gpin);
        sink_pin_intra_cluster_delay_[pin_id] = it->second;
    }
}

float FlatPlacementArcDelayEstimator::estimate_inter_cluster_arc_delay_(APPinId driver_pin_id,
                                                                        APPinId sink_pin_id,
                                                                        const PartialPlacement& p_placement) const {
    // Get the root locations of the tiles containing the driver and sink
    // blocks. The placer queries the delay model using the root locations of
    // the clusters, so we do the same here.
    // TODO: The from and to pins are not known since they are cluster-level
    //       pins, and the clusters have not been created yet. The delay models
    //       that we care about do not use the pins yet.
    t_physical_tile_loc driver_loc = get_containing_tile_root_loc(ap_netlist_.pin_block(driver_pin_id),
                                                                  p_placement,
                                                                  device_grid_);
    t_physical_tile_loc sink_loc = get_containing_tile_root_loc(ap_netlist_.pin_block(sink_pin_id),
                                                                p_placement,
                                                                device_grid_);
    float routing_delay = place_delay_model_.delay(driver_loc,
                                                   0 /*from_pin*/,
                                                   sink_loc,
                                                   0 /*to_pin*/);

    // The delay model returns ROUTER_LOOKAHEAD_NO_PATH_SENTINEL when it has
    // no entry for this driver/sink pair (a gap in the model). Use a
    // distance-based estimate so the arc is not treated as free.
    if (routing_delay >= ROUTER_LOOKAHEAD_NO_PATH_SENTINEL) {
        int manhattan_dist = std::abs(driver_loc.x - sink_loc.x)
                             + std::abs(driver_loc.y - sink_loc.y);
        routing_delay = manhattan_dist * delay_per_tile_;
    }

    return pin_cluster_boundary_delay_[driver_pin_id]
           + routing_delay
           + pin_cluster_boundary_delay_[sink_pin_id];
}

float FlatPlacementArcDelayEstimator::estimate_arc_delay(APPinId sink_pin_id,
                                                         const PartialPlacement& p_placement) const {
    VTR_ASSERT_SAFE(ap_netlist_.pin_type(sink_pin_id) == PinType::SINK);
    APNetId net_id = ap_netlist_.pin_net(sink_pin_id);

    // Global nets are not routed through the general routing network.
    // TODO: This is only true for ideal clock modeling. With other clock
    //       modeling options (e.g. route or dedicated_network) clock nets
    //       are routed and will have a delay.
    // NOTE: The AP netlist speculatively marks any net connected to a clock
    //       port or a non-clock global port as global.
    if (ap_netlist_.net_is_global(net_id))
        return 0.0f;

    // Constant nets which are not routed (see --constant_net_method) have no
    // routing delay.
    if (ap_netlist_.net_is_constant(net_id) && ap_netlist_.net_is_ignored(net_id))
        return 0.0f;

    APPinId driver_pin_id = ap_netlist_.net_driver(net_id);
    VTR_ASSERT_SAFE_MSG(driver_pin_id.is_valid(),
                        "Cannot estimate the delay of an arc without a driver");

    // If the driver and sink blocks are in the same tile, assume that they
    // will be packed into the same cluster and use the intra-cluster delay
    // between them.
    // NOTE: This ignores the case where two blocks in the same tile are in
    //       different sub-tiles (for example, two IO blocks placed in the same
    //       IO tile). In that case there is usually no intra-cluster path
    //       between the pins, and the inter-cluster estimate below is used.
    float intra_cluster_delay = sink_pin_intra_cluster_delay_[sink_pin_id];
    if (intra_cluster_delay >= 0.0f) {
        APBlockId driver_blk_id = ap_netlist_.pin_block(driver_pin_id);
        APBlockId sink_blk_id = ap_netlist_.pin_block(sink_pin_id);
        t_physical_tile_loc driver_loc = get_containing_tile_root_loc(driver_blk_id, p_placement, device_grid_);
        t_physical_tile_loc sink_loc = get_containing_tile_root_loc(sink_blk_id, p_placement, device_grid_);
        if (driver_loc == sink_loc)
            return intra_cluster_delay;
    }

    // Otherwise, the arc must be routed between clusters.
    return estimate_inter_cluster_arc_delay_(driver_pin_id, sink_pin_id, p_placement);
}

void FlatPlacementArcDelayEstimator::update_arc_delays(const PartialPlacement& p_placement,
                                                       PreClusterDelayCalculator& delay_calc) const {
    // For each AP sink pin, update the delay of the timing arc going through it.
    // The delay calculator operates on the Atom netlist; however, by construction
    // of the AP netlist, every atom pin corresponds 1to1 to an AP pin.
    for (APPinId ap_pin_id : ap_netlist_.pins()) {
        // Timing arcs are uniquely identified by the sink pin. Only update
        // timing for sink pins.
        if (ap_netlist_.pin_type(ap_pin_id) != PinType::SINK)
            continue;

        float delay = estimate_arc_delay(ap_pin_id, p_placement);
        delay_calc.set_arc_delay(ap_netlist_.pin_atom_pin(ap_pin_id), delay);
    }
}

void update_timing_info_with_flat_placement(PreClusterTimingManager& pre_cluster_timing_manager,
                                            const FlatPlacementArcDelayEstimator* arc_delay_estimator,
                                            const PartialPlacement& p_placement) {
    // If the timing manager is invalid (i.e. timing analysis is off), do not
    // update.
    if (!pre_cluster_timing_manager.is_valid())
        return;
    VTR_ASSERT_SAFE(arc_delay_estimator != nullptr);

    // Re-estimate the delays of all timing arcs using the flat placement.
    arc_delay_estimator->update_arc_delays(p_placement,
                                           *pre_cluster_timing_manager.get_delay_calculator_ptr());

    // If the timing update type is incremental, we need to invalidate all edges which have changed.
    // We assume here that all edge delays change in some way. We could do a more complicated
    // check for each edge modified and check if the delay has changed; but that may likely
    // take more time than just invalidating all of the edges.
    // Since this loop iterates over all of the edges in the timing graph, we only do this if incremental
    // is selected.
    if (pre_cluster_timing_manager.get_timing_update_type() == e_timing_update_type::INCREMENTAL) {
        for (tatum::EdgeId edge : pre_cluster_timing_manager.get_timing_info().timing_graph()->edges()) {
            pre_cluster_timing_manager.get_timing_info_ptr()->invalidate_delay(edge);
        }
    }

    // Update the timing info. This will run STA to recompute the slacks and
    // the criticalities of all timing arcs.
    pre_cluster_timing_manager.update_timing_info();

    // Do not warn again about unconstrained nodes during placement.
    // Without this line, every GP iteration would see the same warning.
    // Ok to warn once after the first iteration.
    pre_cluster_timing_manager.get_timing_info_ptr()->set_warn_unconstrained(false);
}
