#pragma once

#include "vtr_assert.h"

#include "tatum/Time.hpp"
#include "tatum/delay_calc/DelayCalculator.hpp"
#include "tatum/TimingGraph.hpp"

#include "vpr_error.h"
#include "vpr_utils.h"

#include "atom_netlist.h"
#include "atom_lookup.h"
#include "logic_types.h"
#include "physical_types.h"
#include "prepack.h"

class LogicalModels;

class PreClusterDelayCalculator : public tatum::DelayCalculator {
  public:
    PreClusterDelayCalculator(const AtomNetlist& netlist,
                              const AtomLookup& netlist_lookup,
                              const LogicalModels& models,
                              float intercluster_net_delay,
                              const Prepacker& prepacker) noexcept
        : netlist_(netlist)
        , netlist_lookup_(netlist_lookup)
        , models_(models)
        , inter_cluster_net_delay_(intercluster_net_delay)
        , prepacker_(prepacker) {
        //nop
    }

    tatum::Time max_edge_delay(const tatum::TimingGraph& tg, tatum::EdgeId edge_id) const override {
        tatum::NodeId src_node = tg.edge_src_node(edge_id);
        tatum::NodeId sink_node = tg.edge_sink_node(edge_id);

        auto edge_type = tg.edge_type(edge_id);

        if (edge_type == tatum::EdgeType::PRIMITIVE_COMBINATIONAL) {
            return prim_comb_delay(tg, src_node, sink_node);
        } else if (edge_type == tatum::EdgeType::PRIMITIVE_CLOCK_LAUNCH) {
            return prim_tcq_delay(tg, src_node, sink_node);
        } else {
            VTR_ASSERT(edge_type == tatum::EdgeType::INTERCONNECT);

            //External net delay
            return tatum::Time(inter_cluster_net_delay_);
        }
    }

    tatum::Time setup_time(const tatum::TimingGraph& tg, tatum::EdgeId edge_id) const override {
        tatum::NodeId src_node = tg.edge_src_node(edge_id);
        tatum::NodeId sink_node = tg.edge_sink_node(edge_id);
        auto edge_type = tg.edge_type(edge_id);

        VTR_ASSERT_MSG(tg.node_type(src_node) == tatum::NodeType::CPIN, "Edge setup time only valid if source node is a CPIN");
        VTR_ASSERT_MSG(tg.node_type(sink_node) == tatum::NodeType::SINK, "Edge setup time only valid if sink node is a SINK");
        VTR_ASSERT(edge_type == tatum::EdgeType::PRIMITIVE_CLOCK_CAPTURE);

        AtomPinId sink_pin = netlist_lookup_.tnode_atom_pin(sink_node);
        VTR_ASSERT(sink_pin);

        const t_pb_graph_pin* gpin = find_pb_graph_pin(sink_pin);
        VTR_ASSERT(gpin->type == PB_PIN_SEQUENTIAL);

        return tatum::Time(gpin->tsu);
    }

    tatum::Time min_edge_delay(const tatum::TimingGraph& tg, tatum::EdgeId edge_id) const override {
        //Currently return the same delay
        //TODO: use true min delay
        return max_edge_delay(tg, edge_id);
    }

    tatum::Time hold_time(const tatum::TimingGraph& tg, tatum::EdgeId edge_id) const override {
        //Currently return the same as hold time
        //TODO: use true hold time
        return setup_time(tg, edge_id);
    }

  private:
    //TODO: use generic AtomDelayCalc class to avoid code duplication

    tatum::Time prim_tcq_delay(const tatum::TimingGraph& tg, tatum::NodeId src_node, tatum::NodeId sink_node) const {
        VTR_ASSERT_MSG(tg.node_type(src_node) == tatum::NodeType::CPIN
                           && tg.node_type(sink_node) == tatum::NodeType::SOURCE,
                       "Tcq only defined from CPIN to SOURCE");

        AtomPinId sink_pin = netlist_lookup_.tnode_atom_pin(sink_node);
        VTR_ASSERT(sink_pin);

        const t_pb_graph_pin* gpin = find_pb_graph_pin(sink_pin);
        VTR_ASSERT(gpin->type == PB_PIN_SEQUENTIAL);

        //Clock-to-q delay marked on the SOURCE node (the sink node of this edge)
        auto tco = tatum::Time(gpin->tco_max);

        VTR_ASSERT_MSG(tco.valid(), "Found no primitive clock-to-q delay");

        return tco;
    }

    tatum::Time prim_comb_delay(const tatum::TimingGraph& tg, tatum::NodeId src_node, tatum::NodeId sink_node) const {
        auto src_node_type = tg.node_type(src_node);
        auto sink_node_type = tg.node_type(sink_node);
        VTR_ASSERT_MSG((src_node_type == tatum::NodeType::IPIN && sink_node_type == tatum::NodeType::OPIN)
                           || (src_node_type == tatum::NodeType::SOURCE && sink_node_type == tatum::NodeType::SINK)
                           || (src_node_type == tatum::NodeType::SOURCE && sink_node_type == tatum::NodeType::OPIN)
                           || (src_node_type == tatum::NodeType::CPIN && sink_node_type == tatum::NodeType::OPIN)
                           || (src_node_type == tatum::NodeType::IPIN && sink_node_type == tatum::NodeType::SINK),
                       "Primitive combinational delay must be between {SOURCE, IPIN} and {SINK, OPIN}, or CPIN/OPIN");

        //Primitive internal combinational delay
        AtomPinId input_pin = netlist_lookup_.tnode_atom_pin(src_node);
        VTR_ASSERT(input_pin);
        const t_pb_graph_pin* input_gpin = find_pb_graph_pin(input_pin);

        AtomPinId output_pin = netlist_lookup_.tnode_atom_pin(sink_node);
        VTR_ASSERT(output_pin);
        const t_pb_graph_pin* output_gpin = find_pb_graph_pin(output_pin);

        tatum::Time time;
        for (int i = 0; i < input_gpin->num_pin_timing; ++i) {
            const t_pb_graph_pin* sink_gpin = input_gpin->pin_timing[i];

            if (sink_gpin == output_gpin) {
                time = tatum::Time(input_gpin->pin_timing_del_max[i]);
                break;
            }
        }

        VTR_ASSERT_MSG(time.valid(), "Found no primitive combinational delay for edge");

        return time;
    }

    const t_pb_graph_pin* find_pb_graph_pin(const AtomPinId pin) const {
        AtomBlockId blk = netlist_.pin_block(pin);

        const t_pb_graph_node* pb_gnode = prepacker_.get_expected_lowest_cost_pb_gnode(blk);

        AtomPortId port = netlist_.pin_port(pin);
        const t_model_ports* model_port = netlist_.port_model(port);
        int ipin = netlist_.pin_port_bit(pin);

        const t_pb_graph_pin* gpin = get_pb_graph_node_pin_from_model_port_pin(model_port, ipin, pb_gnode);
        VTR_ASSERT(gpin);

        return gpin;
    }

    const t_pb_graph_pin* find_associated_clock_pin(const AtomPinId io_pin) const {
        const t_pb_graph_pin* io_gpin = find_pb_graph_pin(io_pin);

        const t_pb_graph_pin* clock_gpin = io_gpin->associated_clock_pin;

        if (!clock_gpin) {
            AtomBlockId blk = netlist_.pin_block(io_pin);
            std::string model_name = models_.get_model(netlist_.block_model(blk)).name;
            VPR_FATAL_ERROR(VPR_ERROR_TIMING, "Failed to find clock pin associated with pin '%s' (model '%s')", netlist_.pin_name(io_pin).c_str(), model_name.c_str());
        }
        return clock_gpin;
    }

  private:
    const AtomNetlist& netlist_;
    const AtomLookup& netlist_lookup_;
    const LogicalModels& models_;
    const float inter_cluster_net_delay_;
    const Prepacker& prepacker_;
};
