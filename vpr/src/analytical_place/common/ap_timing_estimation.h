#pragma once
/**
 * @file
 * @author  Alex Singer
 * @date    October 2026
 * @brief   Methods for estimating the delays of timing arcs using a flat
 *          placement.
 *
 * These estimates are used both to update the timing information during
 * global placement and to estimate the post-routing timing of a flat
 * placement.
 */

#include <string>
#include "ap_netlist_fwd.h"
#include "vtr_vector.h"

class APNetlist;
class DeviceGrid;
class PlaceDelayModel;
class PreClusterDelayCalculator;
class PreClusterTimingManager;
struct PartialPlacement;
namespace tatum {
class TimingGraph;
}

/**
 * @brief How the arc delay estimator computes the delay of a timing arc.
 */
enum class e_flat_placement_arc_type {
    UNROUTED_GLOBAL,   ///< The arc is on a global net, which is not routed.
    UNROUTED_CONSTANT, ///< The arc is on a constant net which is not routed.
    INTRA_CLUSTER,     ///< The driver and sink are in the same tile with an intra-cluster path between them.
    INTER_CLUSTER      ///< The arc is routed between clusters.
};

/**
 * @brief Estimator for the post-routing delays of the timing arcs between
 *        AP blocks, based on a flat placement.
 *
 * Timing arcs are identified by the AP sink pin that they terminate at. Every
 * AP sink pin corresponds to an atom sink pin, which identifies the timing arc
 * in the pre-cluster delay calculator.
 *
 * The delay of an arc is estimated as follows:
 *  - Arcs on global nets and unrouted constant nets have no delay, since these
 *    nets are not routed through the general routing network.
 *  - If the driver and sink blocks are in the same tile, they are assumed to
 *    be clustered together. The arc is given the delay of the intra-cluster
 *    path between the two primitive pins (if such a path exists).
 *  - Otherwise, the arc is given the intra-cluster delay from the driver
 *    primitive pin to the boundary of its cluster, plus the inter-tile routing
 *    delay from the place delay model, plus the intra-cluster delay from the
 *    boundary of the sink cluster to the sink primitive pin.
 *
 * Arcs between atoms in the same molecule or the same chain are handled by the
 * pre-cluster delay calculator itself and are not affected by these estimates.
 *
 * All placement-independent quantities (the intra-cluster delays) are computed
 * once on construction, so updating the arc delays for a new placement is
 * cheap.
 */
class FlatPlacementArcDelayEstimator {
  public:
    /**
     * @brief Constructor for the arc delay estimator.
     *
     *  @param ap_netlist
     *          The AP netlist that placements will be over.
     *  @param delay_calc
     *          The pre-cluster delay calculator. Used to find the pb_graph pins
     *          that each primitive pin is expected to be implemented by, so
     *          that this estimator agrees with the delay calculator.
     *  @param place_delay_model
     *          The delay model used to estimate inter-tile routing delays.
     *  @param device_grid
     *          The device grid that placements will be over.
     */
    FlatPlacementArcDelayEstimator(const APNetlist& ap_netlist,
                                   const PreClusterDelayCalculator& delay_calc,
                                   const PlaceDelayModel& place_delay_model,
                                   const DeviceGrid& device_grid);

    /**
     * @brief Estimate the delay of the timing arc which terminates at the
     *        given AP sink pin, using the given flat placement.
     *
     *  @param sink_pin_id  The AP sink pin which identifies the timing arc.
     *  @param p_placement  The flat placement to estimate the delay with.
     *
     *  @return The estimated delay of the timing arc in seconds.
     */
    float estimate_arc_delay(APPinId sink_pin_id,
                             const PartialPlacement& p_placement) const;

    /**
     * @brief Get how the delay of the timing arc which terminates at the
     *        given AP sink pin is estimated, using the given flat placement.
     */
    e_flat_placement_arc_type get_arc_type(APPinId sink_pin_id,
                                           const PartialPlacement& p_placement) const;

    /**
     * @brief Re-estimate the delays of all timing arcs in the AP netlist using
     *        the given flat placement and store them in the delay calculator.
     *
     * This does not perform STA; the caller must update the timing info which
     * uses the delay calculator afterwards.
     *
     *  @param p_placement  The flat placement to estimate the delays with.
     *  @param delay_calc   The delay calculator to store the arc delays in.
     */
    void update_arc_delays(const PartialPlacement& p_placement,
                           PreClusterDelayCalculator& delay_calc) const;

    /**
     * @brief Write information on every interconnect timing arc to the given
     *        file, for debugging the accuracy of the arc delay estimates.
     *
     * One line is written per interconnect edge in the timing graph, with the
     * edge ID, how its delay is computed (by the delay calculator or by this
     * estimator), the fanout of its net, and the pb_graph pins expected to
     * implement its source and sink. The edge IDs match the edge IDs in the
     * timing graph echo files.
     *
     *  @param filename     The file to write to.
     *  @param p_placement  The flat placement the delays were estimated with.
     *  @param delay_calc   The delay calculator the delays were stored in.
     *  @param timing_graph The timing graph the delay calculator is over.
     */
    void write_arc_info(const std::string& filename,
                        const PartialPlacement& p_placement,
                        const PreClusterDelayCalculator& delay_calc,
                        const tatum::TimingGraph& timing_graph) const;

  private:
    /**
     * @brief Estimate the delay of an arc between the given driver and sink
     *        pins when they are in different clusters.
     *
     * This is the intra-cluster delay out of the driver cluster, plus the
     * inter-tile routing delay, plus the intra-cluster delay into the sink
     * cluster.
     */
    float estimate_inter_cluster_arc_delay_(APPinId driver_pin_id,
                                            APPinId sink_pin_id,
                                            const PartialPlacement& p_placement) const;

    /// @brief The AP netlist that placements are over.
    const APNetlist& ap_netlist_;

    /// @brief The delay model used to estimate inter-tile routing delays.
    const PlaceDelayModel& place_delay_model_;

    /// @brief The device grid that placements are over.
    const DeviceGrid& device_grid_;

    /// @brief A representative routing delay per tile of distance. Used as a
    ///        fallback when the place delay model has no entry for a pair
    ///        of tiles.
    float delay_per_tile_;

    /// @brief The intra-cluster delay between each AP pin and the boundary of
    ///        the cluster that its block will be packed into. For driver pins,
    ///        this is the delay from the pin to the cluster's output pins; for
    ///        sink pins, this is the delay from the cluster's input pins to
    ///        the pin.
    vtr::vector<APPinId, float> pin_cluster_boundary_delay_;

    /// @brief The intra-cluster delay between the driver of each AP sink pin's
    ///        net and the sink pin, assuming both blocks are packed into the
    ///        same cluster. This is negative if there is no intra-cluster path
    ///        between the pins (or if the pin is not a sink pin).
    vtr::vector<APPinId, float> sink_pin_intra_cluster_delay_;
};

/**
 * @brief Update the timing information in the pre-cluster timing manager
 *        using a flat placement as a hint for where the atoms will be placed.
 *
 * The delays of all timing arcs are re-estimated from the flat placement and
 * STA is performed to recompute the slacks and criticalities.
 *
 * If the timing manager is invalid (i.e. timing analysis is off), this does
 * nothing.
 *
 *  @param pre_cluster_timing_manager
 *      Manager object which computes the slacks of timing edges.
 *  @param arc_delay_estimator
 *      Estimator used to compute the delays of timing arcs from the flat
 *      placement. May be nullptr if the timing manager is invalid.
 *  @param p_placement
 *      The flat placement used to update the timing information.
 */
void update_timing_info_with_flat_placement(PreClusterTimingManager& pre_cluster_timing_manager,
                                            const FlatPlacementArcDelayEstimator* arc_delay_estimator,
                                            const PartialPlacement& p_placement);
