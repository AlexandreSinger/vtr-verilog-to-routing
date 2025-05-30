/*
 * Jason Luu 2008
 * Print complex block information to a file
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string_view>

#include "cluster_legalizer.h"
#include "clustered_netlist.h"
#include "physical_types.h"
#include "physical_types_util.h"
#include "prepack.h"
#include "vpr_context.h"
#include "vtr_assert.h"
#include "vtr_log.h"

#include "vpr_types.h"
#include "vpr_error.h"

#include "pugixml.hpp"

#include "globals.h"
#include "atom_netlist.h"
#include "pb_type_graph.h"
#include "output_clustering.h"
#include "vpr_utils.h"
#include "pack.h"

static void print_clustering_stats_header();
static void print_clustering_stats(std::string_view block_name, int num_block_type, float num_inputs_clocks, float num_outputs);

/**************** Subroutine definitions ************************************/

static void count_clb_inputs_and_outputs_from_pb_route(const t_pb* pb,
                                                       t_logical_block_type_ptr logical_block,
                                                       int ipin,
                                                       e_pin_type pin_type,
                                                       std::unordered_map<AtomNetId, bool>& nets_absorbed,
                                                       int num_clb_inputs_used[],
                                                       int num_clb_outputs_used[]) {
    VTR_ASSERT_DEBUG(!pb->pb_route.empty());
    int pb_graph_pin_id = get_pb_graph_node_pin_from_pb_graph_node(pb->pb_graph_node, ipin)->pin_count_in_cluster;

    if (pb->pb_route.count(pb_graph_pin_id)) {
        //Pin used
        AtomNetId atom_net_id = pb->pb_route[pb_graph_pin_id].atom_net_id;
        if (atom_net_id) {
            nets_absorbed[atom_net_id] = false;
            if (pin_type == RECEIVER) {
                num_clb_inputs_used[logical_block->index]++;
            } else if (pin_type == DRIVER) {
                num_clb_outputs_used[logical_block->index]++;
            }
        }
    }
}

static void count_stats_from_legalizer(const ClusterLegalizer& cluster_legalizer,
                                       std::unordered_map<AtomNetId, bool>& nets_absorbed,
                                       int num_clb_types[],
                                       int num_clb_inputs_used[],
                                       int num_clb_outputs_used[]) {
    for (LegalizationClusterId cluster_id : cluster_legalizer.clusters()) {
        t_logical_block_type_ptr logical_block = cluster_legalizer.get_cluster_type(cluster_id);
        t_physical_tile_type_ptr physical_tile = pick_physical_type(logical_block);
        for (int ipin = 0; ipin < logical_block->pb_type->num_pins; ipin++) {
            int physical_pin = get_physical_pin(physical_tile, logical_block, ipin);
            e_pin_type pin_type = get_pin_type_from_pin_physical_num(physical_tile, physical_pin);

            const t_pb* pb = cluster_legalizer.get_cluster_pb(cluster_id);
            if (pb->pb_route.empty())
                continue;
            count_clb_inputs_and_outputs_from_pb_route(pb,
                                                       logical_block,
                                                       ipin,
                                                       pin_type,
                                                       nets_absorbed,
                                                       num_clb_inputs_used,
                                                       num_clb_outputs_used);
        }
        num_clb_types[logical_block->index]++;
    }
}

static void count_stats_from_netlist(std::unordered_map<AtomNetId, bool>& nets_absorbed,
                                     int num_clb_types[],
                                     int num_clb_inputs_used[],
                                     int num_clb_outputs_used[]) {
    const AtomContext& atom_ctx = g_vpr_ctx.atom();
    const ClusteredNetlist& clb_nlist = g_vpr_ctx.clustering().clb_nlist;

    for (ClusterBlockId blk_id : clb_nlist.blocks()) {
        t_logical_block_type_ptr logical_block = clb_nlist.block_type(blk_id);
        t_physical_tile_type_ptr physical_tile = pick_physical_type(logical_block);
        for (int ipin = 0; ipin < logical_block->pb_type->num_pins; ipin++) {
            int physical_pin = get_physical_pin(physical_tile, logical_block, ipin);
            e_pin_type pin_type = get_pin_type_from_pin_physical_num(physical_tile, physical_pin);

            if (!clb_nlist.block_pb(blk_id)->pb_route.empty()) {
                count_clb_inputs_and_outputs_from_pb_route(clb_nlist.block_pb(blk_id),
                                                           logical_block,
                                                           ipin,
                                                           pin_type,
                                                           nets_absorbed,
                                                           num_clb_inputs_used,
                                                           num_clb_outputs_used);
            } else {
                ClusterNetId clb_net_id = clb_nlist.block_net(blk_id, ipin);
                if (clb_net_id != ClusterNetId::INVALID()) {
                    AtomNetId net_id = atom_ctx.lookup().atom_net(clb_net_id);
                    VTR_ASSERT(net_id);
                    nets_absorbed[net_id] = false;

                    if (pin_type == RECEIVER) {
                        num_clb_inputs_used[logical_block->index]++;
                    } else if (pin_type == DRIVER) {
                        num_clb_outputs_used[logical_block->index]++;
                    }
                }
            }
        }
        num_clb_types[logical_block->index]++;
    }
}

/* Prints out one cluster (clb).  Both the external pins and the *
 * internal connections are printed out.                         */
static void print_stats(const ClusterLegalizer* cluster_legalizer_ptr, bool from_legalizer) {
    const DeviceContext& device_ctx = g_vpr_ctx.device();
    const AtomNetlist& atom_nlist = g_vpr_ctx.atom().netlist();

    int* num_clb_types = new int[device_ctx.logical_block_types.size()];
    int* num_clb_inputs_used = new int[device_ctx.logical_block_types.size()];
    int* num_clb_outputs_used = new int[device_ctx.logical_block_types.size()];

    for (size_t i = 0; i < device_ctx.logical_block_types.size(); i++) {
        num_clb_types[i] = 0;
        num_clb_inputs_used[i] = 0;
        num_clb_outputs_used[i] = 0;
    }

    std::unordered_map<AtomNetId, bool> nets_absorbed;
    for (AtomNetId net_id : atom_nlist.nets()) {
        nets_absorbed[net_id] = true;
    }

    /* Counters used only for statistics purposes. */
    if (from_legalizer) {
        VTR_ASSERT(cluster_legalizer_ptr != nullptr);
        count_stats_from_legalizer(*cluster_legalizer_ptr, nets_absorbed, num_clb_types, num_clb_inputs_used, num_clb_outputs_used);
    } else {
        VTR_ASSERT(cluster_legalizer_ptr == nullptr);
        count_stats_from_netlist(nets_absorbed, num_clb_types, num_clb_inputs_used, num_clb_outputs_used);
    }

    print_clustering_stats_header();

    for (unsigned int itype = 0; itype < device_ctx.logical_block_types.size(); itype++) {
        if (num_clb_types[itype] == 0) {
            print_clustering_stats(device_ctx.logical_block_types[itype].name, num_clb_types[itype], 0.0, 0.0);
        } else {
            print_clustering_stats(device_ctx.logical_block_types[itype].name, num_clb_types[itype],
                                   (float)num_clb_inputs_used[itype] / (float)num_clb_types[itype],
                                   (float)num_clb_outputs_used[itype] / (float)num_clb_types[itype]);
        }
    }

    int total_nets_absorbed = 0;
    for (AtomNetId net_id : atom_nlist.nets()) {
        if (nets_absorbed[net_id] == true) {
            total_nets_absorbed++;
        }
    }
    VTR_LOG("Absorbed logical nets %d out of %d nets, %d nets not absorbed.\n",
            total_nets_absorbed, (int)atom_nlist.nets().size(), (int)atom_nlist.nets().size() - total_nets_absorbed);
    delete[] num_clb_types;
    delete[] num_clb_inputs_used;
    delete[] num_clb_outputs_used;
    /* TODO: print more stats */
}

static void print_clustering_stats_header() {
    VTR_LOG("Final Clustering Statistics: \n");
    VTR_LOG("----------   --------   ------------------------------------   --------------------------\n");
    VTR_LOG("Block Type   # Blocks   Avg. # of input clocks and pins used   Avg. # of output pins used\n");
    VTR_LOG("----------   --------   ------------------------------------   --------------------------\n");
}

static void print_clustering_stats(std::string_view block_name, int num_block_type, float num_inputs_clocks, float num_outputs) {
    VTR_LOG(
        "%10s   "
        "%8d   "
        "%36g   "
        "%26g   ",
        block_name.data(),
        num_block_type,
        num_inputs_clocks,
        num_outputs);

    VTR_LOG("\n");

    fflush(stdout);
}

static const char* clustering_xml_net_text(AtomNetId net_id) {
    /* This routine prints out the atom_ctx.netlist() net name (or open).
     * net_num is the index of the atom_ctx.netlist() net to be printed
     */
    const AtomNetlist& atom_nlist = g_vpr_ctx.atom().netlist();

    if (!net_id) {
        return "open";
    } else {
        return atom_nlist.net_name(net_id).c_str();
    }
}

static std::string clustering_xml_interconnect_text(t_logical_block_type_ptr type, const IntraLbPbPinLookup& pb_graph_pin_lookup_from_index_by_type, int inode, const t_pb_routes& pb_route) {
    if (!pb_route.count(inode) || !pb_route[inode].atom_net_id) {
        return "open";
    }

    int prev_node = pb_route[inode].driver_pb_pin_id;
    int prev_edge;
    if (prev_node == OPEN) {
        /* No previous driver implies that this is either a top-level input pin or a primitive output pin */
        const t_pb_graph_pin* cur_pin = pb_graph_pin_lookup_from_index_by_type.pb_gpin(type->index, inode);
        VTR_ASSERT(cur_pin->parent_node->pb_type->is_root() || (cur_pin->is_primitive_pin() && cur_pin->port->type == OUT_PORT));
        return clustering_xml_net_text(pb_route[inode].atom_net_id);
    } else {
        const t_pb_graph_pin* cur_pin = pb_graph_pin_lookup_from_index_by_type.pb_gpin(type->index, inode);
        const t_pb_graph_pin* prev_pin = pb_graph_pin_lookup_from_index_by_type.pb_gpin(type->index, prev_node);

        for (prev_edge = 0; prev_edge < prev_pin->num_output_edges; prev_edge++) {
            VTR_ASSERT(prev_pin->output_edges[prev_edge]->num_output_pins == 1);
            if (prev_pin->output_edges[prev_edge]->output_pins[0]->pin_count_in_cluster == inode) {
                break;
            }
        }
        VTR_ASSERT(prev_edge < prev_pin->num_output_edges);

        char* name = prev_pin->output_edges[prev_edge]->interconnect->name;
        if (prev_pin->port->parent_pb_type->depth
            >= cur_pin->port->parent_pb_type->depth) {
            /* Connections from siblings or children should have an explicit index, connections from parent does not need an explicit index */
            return vtr::string_fmt("%s[%d].%s[%d]->%s",
                                   prev_pin->parent_node->pb_type->name,
                                   prev_pin->parent_node->placement_index,
                                   prev_pin->port->name,
                                   prev_pin->pin_number, name);
        } else {
            return vtr::string_fmt("%s.%s[%d]->%s",
                                   prev_pin->parent_node->pb_type->name,
                                   prev_pin->port->name,
                                   prev_pin->pin_number, name);
        }
    }
}

/* outputs a block that is open or unused.
 * In some cases, a block is unused for logic but is used for routing. When that happens, the block
 * cannot simply be marked open as that would lose the routing information. Instead, a block must be
 * output that reflects the routing resources used. This function handles both cases.
 */
static void clustering_xml_open_block(pugi::xml_node& parent_node, t_logical_block_type_ptr type, const IntraLbPbPinLookup& pb_graph_pin_lookup_from_index_by_type, t_pb_graph_node* pb_graph_node, int pb_index, bool is_used, const t_pb_routes& pb_route) {
    int i, j, k, m;
    const t_pb_type *pb_type, *child_pb_type;
    t_mode* mode = nullptr;
    int prev_node;
    int mode_of_edge, port_index, node_index;

    mode_of_edge = UNDEFINED;

    pb_type = pb_graph_node->pb_type;

    pugi::xml_node block_node = parent_node.append_child("block");
    block_node.append_attribute("name") = "open";
    block_node.append_attribute("instance") = vtr::string_fmt("%s[%d]", pb_graph_node->pb_type->name, pb_index).c_str();
    std::vector<std::string> block_modes;

    if (is_used) {
        /* Determine mode if applicable */
        port_index = 0;
        for (i = 0; i < pb_type->num_ports; i++) {
            if (pb_type->ports[i].type == OUT_PORT) {
                VTR_ASSERT(!pb_type->ports[i].is_clock);
                for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                    const t_pb_graph_pin* pin = &pb_graph_node->output_pins[port_index][j];
                    node_index = pin->pin_count_in_cluster;
                    if (!pb_type->is_primitive() && pb_route.count(node_index) && pb_route[node_index].atom_net_id) {
                        prev_node = pb_route[node_index].driver_pb_pin_id;
                        const t_pb_graph_pin* prev_pin = pb_graph_pin_lookup_from_index_by_type.pb_gpin(type->index, prev_node);
                        const t_pb_graph_edge* edge = get_edge_between_pins(prev_pin, pin);

                        VTR_ASSERT(edge != nullptr);
                        mode_of_edge = edge->interconnect->parent_mode_index;
                        if (mode != nullptr && &pb_type->modes[mode_of_edge] != mode) {
                            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                                            "Differing modes for block.  Got %s previously and %s for edge %d (interconnect %s).",
                                            mode->name, pb_type->modes[mode_of_edge].name,
                                            port_index,
                                            edge->interconnect->name);
                        }
                        VTR_ASSERT(mode == nullptr || &pb_type->modes[mode_of_edge] == mode);
                        mode = &pb_type->modes[mode_of_edge];
                    }
                }
                port_index++;
            }
        }

        VTR_ASSERT(mode != nullptr && mode_of_edge != UNDEFINED);

        block_node.append_attribute("mode") = mode->name;
        block_node.append_attribute("pb_type_num_modes") = pb_type->num_modes;

        pugi::xml_node inputs_node = block_node.append_child("inputs");

        port_index = 0;
        for (i = 0; i < pb_type->num_ports; i++) {
            if (!pb_type->ports[i].is_clock && pb_type->ports[i].type == IN_PORT) {
                pugi::xml_node port_node = inputs_node.append_child("port");
                port_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;

                std::vector<std::string> pins;
                for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                    node_index = pb_graph_node->input_pins[port_index][j].pin_count_in_cluster;

                    if (pb_type->is_root()) {
                        pins.push_back(clustering_xml_net_text(pb_route[node_index].atom_net_id));
                    } else {
                        pins.push_back(clustering_xml_interconnect_text(type, pb_graph_pin_lookup_from_index_by_type, node_index, pb_route));
                    }
                }
                port_node.text().set(vtr::join(pins.begin(), pins.end(), " ").c_str());
                port_index++;
            }
        }

        pugi::xml_node outputs_node = block_node.append_child("outputs");

        port_index = 0;
        for (i = 0; i < pb_type->num_ports; i++) {
            if (pb_type->ports[i].type == OUT_PORT) {
                VTR_ASSERT(!pb_type->ports[i].is_clock);

                pugi::xml_node port_node = outputs_node.append_child("port");
                port_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;
                std::vector<std::string> pins;
                for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                    node_index = pb_graph_node->output_pins[port_index][j].pin_count_in_cluster;
                    pins.push_back(clustering_xml_interconnect_text(type, pb_graph_pin_lookup_from_index_by_type, node_index, pb_route));
                }
                port_node.text().set(vtr::join(pins.begin(), pins.end(), " ").c_str());
                port_index++;
            }
        }

        pugi::xml_node clock_node = block_node.append_child("clocks");

        port_index = 0;
        for (i = 0; i < pb_type->num_ports; i++) {
            if (pb_type->ports[i].is_clock && pb_type->ports[i].type == IN_PORT) {
                pugi::xml_node port_node = clock_node.append_child("port");
                port_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;

                std::vector<std::string> pins;
                for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                    node_index = pb_graph_node->clock_pins[port_index][j].pin_count_in_cluster;
                    if (pb_type->is_root()) {
                        pins.push_back(clustering_xml_net_text(pb_route[node_index].atom_net_id));
                    } else {
                        pins.push_back(clustering_xml_interconnect_text(type, pb_graph_pin_lookup_from_index_by_type, node_index, pb_route));
                    }
                }
                port_node.text().set(vtr::join(pins.begin(), pins.end(), " ").c_str());
                port_index++;
            }
        }

        if (!pb_type->is_primitive()) {
            for (i = 0; i < mode->num_pb_type_children; i++) {
                child_pb_type = &mode->pb_type_children[i];
                for (j = 0; j < mode->pb_type_children[i].num_pb; j++) {
                    port_index = 0;
                    is_used = false;
                    for (k = 0; k < child_pb_type->num_ports && !is_used; k++) {
                        if (child_pb_type->ports[k].type == OUT_PORT) {
                            for (m = 0; m < child_pb_type->ports[k].num_pins; m++) {
                                node_index = pb_graph_node->child_pb_graph_nodes[mode_of_edge][i][j].output_pins[port_index][m].pin_count_in_cluster;
                                if (pb_route.count(node_index) && pb_route[node_index].atom_net_id) {
                                    is_used = true;
                                    break;
                                }
                            }
                            port_index++;
                        }
                    }
                    clustering_xml_open_block(block_node, type, pb_graph_pin_lookup_from_index_by_type,
                                              &pb_graph_node->child_pb_graph_nodes[mode_of_edge][i][j],
                                              j, is_used, pb_route);
                }
            }
        }
    }
}

/* outputs a block that is used (i.e. has configuration) and all of its child blocks */
static void clustering_xml_block(pugi::xml_node& parent_node, t_logical_block_type_ptr type, const IntraLbPbPinLookup& pb_graph_pin_lookup_from_index_by_type, t_pb* pb, int pb_index, const t_pb_routes& pb_route) {
    int i, j, k, m;
    const t_pb_type *pb_type, *child_pb_type;
    t_pb_graph_node* pb_graph_node;
    t_mode* mode;
    int port_index, node_index;
    bool is_used;

    pb_type = pb->pb_graph_node->pb_type;
    pb_graph_node = pb->pb_graph_node;
    mode = &pb_type->modes[pb->mode];

    pugi::xml_node block_node = parent_node.append_child("block");
    block_node.append_attribute("name") = pb->name;
    block_node.append_attribute("instance") = vtr::string_fmt("%s[%d]", pb_type->name, pb_index).c_str();

    if (!pb_type->is_primitive()) {
        block_node.append_attribute("mode") = mode->name;
    } else {
        const auto& atom_ctx = g_vpr_ctx.atom();
        AtomBlockId atom_blk = atom_ctx.netlist().find_block(pb->name);
        VTR_ASSERT(atom_blk);

        pugi::xml_node attrs_node = block_node.append_child("attributes");
        for (const auto& attr : atom_ctx.netlist().block_attrs(atom_blk)) {
            pugi::xml_node attr_node = attrs_node.append_child("attribute");
            attr_node.append_attribute("name") = attr.first.c_str();
            attr_node.text().set(attr.second.c_str());
        }

        pugi::xml_node params_node = block_node.append_child("parameters");
        for (const auto& param : atom_ctx.netlist().block_params(atom_blk)) {
            pugi::xml_node param_node = params_node.append_child("parameter");
            param_node.append_attribute("name") = param.first.c_str();
            param_node.text().set(param.second.c_str());
        }
    }

    pugi::xml_node inputs_node = block_node.append_child("inputs");

    port_index = 0;
    for (i = 0; i < pb_type->num_ports; i++) {
        if (!pb_type->ports[i].is_clock && pb_type->ports[i].type == IN_PORT) {
            pugi::xml_node port_node = inputs_node.append_child("port");
            port_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;

            std::vector<std::string> pins;
            for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                node_index = pb->pb_graph_node->input_pins[port_index][j].pin_count_in_cluster;

                if (pb_type->is_root()) {
                    if (pb_route.count(node_index)) {
                        pins.push_back(clustering_xml_net_text(pb_route[node_index].atom_net_id));
                    } else {
                        pins.push_back(clustering_xml_net_text(AtomNetId::INVALID()));
                    }
                } else {
                    pins.push_back(clustering_xml_interconnect_text(type, pb_graph_pin_lookup_from_index_by_type, node_index, pb_route));
                }
            }
            port_node.text().set(vtr::join(pins.begin(), pins.end(), " ").c_str());

            //The cluster router may have rotated equivalent pins (e.g. LUT inputs),
            //record the resulting rotation here so it can be unambigously mapped
            //back to the atom netlist
            if (pb_type->ports[i].equivalent != PortEquivalence::NONE && pb_type->parent_mode != nullptr && pb_type->is_primitive()) {
                //This is a primitive with equivalent inputs

                auto& atom_ctx = g_vpr_ctx.atom();
                AtomBlockId atom_blk = atom_ctx.netlist().find_block(pb->name);
                VTR_ASSERT(atom_blk);

                AtomPortId atom_port = atom_ctx.netlist().find_atom_port(atom_blk, pb_type->ports[i].model_port);

                if (atom_port) { //Port exists (some LUTs may have no input and hence no port in the atom netlist)

                    pugi::xml_node port_rotation_node = inputs_node.append_child("port_rotation_map");
                    port_rotation_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;

                    std::set<AtomPinId> recorded_pins;
                    std::vector<std::string> pin_map_list;

                    for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                        node_index = pb->pb_graph_node->input_pins[port_index][j].pin_count_in_cluster;

                        if (pb_route.count(node_index)) {
                            AtomNetId atom_net = pb_route[node_index].atom_net_id;

                            VTR_ASSERT(atom_net);

                            //This physical pin is in use, find the original pin in the atom netlist
                            AtomPinId orig_pin;
                            for (AtomPinId atom_pin : atom_ctx.netlist().port_pins(atom_port)) {
                                if (recorded_pins.count(atom_pin)) continue; //Don't add pins twice

                                AtomNetId atom_pin_net = atom_ctx.netlist().pin_net(atom_pin);

                                if (atom_pin_net == atom_net) {
                                    recorded_pins.insert(atom_pin);
                                    orig_pin = atom_pin;
                                    break;
                                }
                            }

                            VTR_ASSERT(orig_pin);
                            //The physical pin j, maps to a pin in the atom netlist
                            pin_map_list.push_back(vtr::string_fmt("%d", atom_ctx.netlist().pin_port_bit(orig_pin)));
                        } else {
                            //The physical pin is disconnected
                            pin_map_list.push_back("open");
                        }
                    }
                    port_rotation_node.text().set(vtr::join(pin_map_list.begin(), pin_map_list.end(), " ").c_str());
                }
            }

            port_index++;
        }
    }

    pugi::xml_node outputs_node = block_node.append_child("outputs");

    port_index = 0;
    for (i = 0; i < pb_type->num_ports; i++) {
        if (pb_type->ports[i].type == OUT_PORT) {
            VTR_ASSERT(!pb_type->ports[i].is_clock);

            pugi::xml_node port_node = outputs_node.append_child("port");
            port_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;
            std::vector<std::string> pins;
            for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                node_index = pb->pb_graph_node->output_pins[port_index][j].pin_count_in_cluster;
                pins.push_back(clustering_xml_interconnect_text(type, pb_graph_pin_lookup_from_index_by_type, node_index, pb_route));
            }
            port_node.text().set(vtr::join(pins.begin(), pins.end(), " ").c_str());
            port_index++;
        }
    }

    pugi::xml_node clock_node = block_node.append_child("clocks");

    port_index = 0;
    for (i = 0; i < pb_type->num_ports; i++) {
        if (pb_type->ports[i].is_clock && pb_type->ports[i].type == IN_PORT) {
            pugi::xml_node port_node = clock_node.append_child("port");
            port_node.append_attribute("name") = pb_graph_node->pb_type->ports[i].name;

            std::vector<std::string> pins;
            for (j = 0; j < pb_type->ports[i].num_pins; j++) {
                node_index = pb->pb_graph_node->clock_pins[port_index][j].pin_count_in_cluster;
                if (pb_type->is_root()) {
                    if (pb_route.count(node_index)) {
                        pins.push_back(clustering_xml_net_text(pb_route[node_index].atom_net_id));
                    } else {
                        pins.push_back(clustering_xml_net_text(AtomNetId::INVALID()));
                    }
                } else {
                    pins.push_back(clustering_xml_interconnect_text(type, pb_graph_pin_lookup_from_index_by_type, node_index, pb_route));
                }
            }
            port_node.text().set(vtr::join(pins.begin(), pins.end(), " ").c_str());
            port_index++;
        }
    }

    if (!pb_type->is_primitive()) {
        for (i = 0; i < mode->num_pb_type_children; i++) {
            for (j = 0; j < mode->pb_type_children[i].num_pb; j++) {
                /* If child pb is not used but routing is used, I must print things differently */
                if ((pb->child_pbs[i] != nullptr) && (pb->child_pbs[i][j].name != nullptr)) {
                    clustering_xml_block(block_node, type, pb_graph_pin_lookup_from_index_by_type, &pb->child_pbs[i][j], j, pb_route);
                } else {
                    is_used = false;
                    child_pb_type = &mode->pb_type_children[i];
                    port_index = 0;

                    for (k = 0; k < child_pb_type->num_ports && !is_used; k++) {
                        if (child_pb_type->ports[k].type == OUT_PORT) {
                            for (m = 0; m < child_pb_type->ports[k].num_pins; m++) {
                                node_index = pb_graph_node->child_pb_graph_nodes[pb->mode][i][j].output_pins[port_index][m].pin_count_in_cluster;
                                if (pb_route.count(node_index) && pb_route[node_index].atom_net_id) {
                                    is_used = true;
                                    break;
                                }
                            }
                            port_index++;
                        }
                    }
                    clustering_xml_open_block(block_node, type, pb_graph_pin_lookup_from_index_by_type,
                                              &pb_graph_node->child_pb_graph_nodes[pb->mode][i][j],
                                              j, is_used, pb_route);
                }
            }
        }
    }
}

static void clustering_xml_blocks_from_legalizer(pugi::xml_node& block_node,
                                                 const IntraLbPbPinLookup& pb_graph_pin_lookup_from_index_by_type,
                                                 ClusterLegalizer& cluster_legalizer) {
    // Finalize the cluster legalization by ensuring that each cluster pb has
    // its pb_route calculated.
    cluster_legalizer.finalize();
    for (LegalizationClusterId cluster_id : cluster_legalizer.clusters()) {
        clustering_xml_block(block_node,
                             cluster_legalizer.get_cluster_type(cluster_id),
                             pb_graph_pin_lookup_from_index_by_type,
                             cluster_legalizer.get_cluster_pb(cluster_id),
                             size_t(cluster_id),
                             cluster_legalizer.get_cluster_pb(cluster_id)->pb_route);
    }
}

static void clustering_xml_blocks_from_netlist(pugi::xml_node& block_node,
                                               const IntraLbPbPinLookup& pb_graph_pin_lookup_from_index_by_type) {
    const ClusteredNetlist& clb_nlist = g_vpr_ctx.clustering().clb_nlist;
    for (auto blk_id : clb_nlist.blocks()) {
        /* TODO: Must do check that total CLB pins match top-level pb pins, perhaps check this earlier? */
        clustering_xml_block(block_node,
                             clb_nlist.block_type(blk_id),
                             pb_graph_pin_lookup_from_index_by_type,
                             clb_nlist.block_pb(blk_id),
                             size_t(blk_id),
                             clb_nlist.block_pb(blk_id)->pb_route);
    }
}

/* This routine dumps out the output netlist in a format suitable for  *
 * input to vpr. This routine also dumps out the internal structure of *
 * the cluster, in essentially a graph based format.                   */
void output_clustering(ClusterLegalizer* cluster_legalizer_ptr, const std::unordered_set<AtomNetId>& is_clock, const std::string& architecture_id, const char* out_fname, bool skip_clustering, bool from_legalizer) {
    const DeviceContext& device_ctx = g_vpr_ctx.device();
    const AtomNetlist& atom_nlist = g_vpr_ctx.atom().netlist();

    IntraLbPbPinLookup pb_graph_pin_lookup_from_index_by_type(device_ctx.logical_block_types);

    pugi::xml_document out_xml;

    pugi::xml_node block_node = out_xml.append_child("block");
    block_node.append_attribute("name") = out_fname;
    block_node.append_attribute("instance") = "FPGA_packed_netlist[0]";
    block_node.append_attribute("architecture_id") = architecture_id.c_str();
    block_node.append_attribute("atom_netlist_id") = atom_nlist.netlist_id().c_str();

    std::vector<std::string> inputs;
    std::vector<std::string> outputs;

    for (auto blk_id : atom_nlist.blocks()) {
        auto type = atom_nlist.block_type(blk_id);
        switch (type) {
            case AtomBlockType::INPAD:
                if (skip_clustering) {
                    VTR_ASSERT(0);
                }
                inputs.push_back(atom_nlist.block_name(blk_id));
                break;

            case AtomBlockType::OUTPAD:
                if (skip_clustering) {
                    VTR_ASSERT(0);
                }
                outputs.push_back(atom_nlist.block_name(blk_id));
                break;

            case AtomBlockType::BLOCK:
                if (skip_clustering) {
                    VTR_ASSERT(0);
                }
                break;

            default:
                VTR_LOG_ERROR("in output_netlist: Unexpected type %d for atom block %s.\n",
                              type, atom_nlist.block_name(blk_id).c_str());
        }
    }

    block_node.append_child("inputs").text().set(vtr::join(inputs.begin(), inputs.end(), " ").c_str());
    block_node.append_child("outputs").text().set(vtr::join(outputs.begin(), outputs.end(), " ").c_str());

    std::vector<std::string> clocks;
    for (auto net_id : atom_nlist.nets()) {
        if (is_clock.count(net_id)) {
            clocks.push_back(atom_nlist.net_name(net_id));
        }
    }

    block_node.append_child("clocks").text().set(vtr::join(clocks.begin(), clocks.end(), " ").c_str());

    if (skip_clustering == false) {
        if (from_legalizer) {
            VTR_ASSERT(cluster_legalizer_ptr != nullptr);
            clustering_xml_blocks_from_legalizer(block_node, pb_graph_pin_lookup_from_index_by_type, *cluster_legalizer_ptr);
        } else {
            VTR_ASSERT(cluster_legalizer_ptr == nullptr);
            clustering_xml_blocks_from_netlist(block_node, pb_graph_pin_lookup_from_index_by_type);
        }
    }

    out_xml.save_file(out_fname);

    print_stats(cluster_legalizer_ptr, from_legalizer);
}

/********************************************************************
 * A useful API to output packing results to a XML file
 * This function is a wrapper for the function output_clustering()
 * but remove all the requirements on input data structures that
 * have to be built with other APIs
 *
 * As such, this function is expected to be a standard API
 * which can be called anytime and anywhere after packing is finished.
 ********************************************************************/
void write_packing_results_to_xml(const std::string& architecture_id,
                                  const char* out_fname) {
    std::unordered_set<AtomNetId> is_clock = alloc_and_load_is_clock();

    // Since the cluster legalizer is not being used to output the clustering
    // (from_legalizer is false), passing in nullptr.
    output_clustering(nullptr,
                      is_clock,
                      architecture_id,
                      out_fname,
                      false, /*skip_clustering*/
                      false /*from_legalizer*/);
}
