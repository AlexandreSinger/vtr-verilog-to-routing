
#include "noc_aware_cluster_util.h"
#include "atom_netlist.h"
#include "globals.h"
#include "logic_types.h"
#include "vpr_types.h"
#include "vpr_utils.h"

#include <queue>

std::vector<AtomBlockId> find_noc_router_atoms(const AtomNetlist& atom_netlist, const LogicalModels& models) {
    // NoC router atoms are expected to have a specific blif model
    LogicalModelId noc_route_blif_model_id = models.get_model_by_name("noc_router_adapter_block");

    // stores found NoC router atoms
    std::vector<AtomBlockId> noc_router_atoms;

    // iterate over all atoms and find those whose blif model matches
    for (auto atom_id : atom_netlist.blocks()) {
        LogicalModelId model_id = atom_netlist.block_model(atom_id);
        if (model_id == noc_route_blif_model_id) {
            noc_router_atoms.push_back(atom_id);
        }
    }

    return noc_router_atoms;
}

void update_noc_reachability_partitions(const std::vector<AtomBlockId>& noc_atoms,
                                        const AtomNetlist& atom_netlist,
                                        const t_pack_high_fanout_thresholds& high_fanout_thresholds,
                                        vtr::vector<AtomBlockId, NocGroupId>& atom_noc_grp_id) {
    const auto& grid = g_vpr_ctx.device().grid;

    t_logical_block_type_ptr logic_block_type = infer_logic_block_type(grid);
    const char* logical_block_name = logic_block_type != nullptr ? logic_block_type->name.c_str() : "";
    const size_t high_fanout_threshold = high_fanout_thresholds.get_threshold(logical_block_name);

    // get the total number of atoms
    const size_t n_atoms = atom_netlist.blocks().size();

    vtr::vector<AtomBlockId, bool> atom_visited(n_atoms, false);

    atom_noc_grp_id.resize(n_atoms, NocGroupId::INVALID());

    int noc_grp_id_cnt = 0;

    /*
     * Assume that the atom netlist is represented as an undirected graph
     * with all high fanout nets removed. In this graph, we want to find all
     * connected components that include at least one NoC router. We start a
     * BFS from each NoC router and traverse all nets below the high_fanout_threshold,
     * and mark each atom block with a NoC group ID.
     */

    for (auto noc_atom_id : noc_atoms) {
        // check if this NoC router has already been visited
        if (atom_visited[noc_atom_id]) {
            continue;
        }

        auto noc_grp_id = (NocGroupId)noc_grp_id_cnt;
        noc_grp_id_cnt++;

        std::queue<AtomBlockId> q;
        q.push(noc_atom_id);
        atom_visited[noc_atom_id] = true;

        while (!q.empty()) {
            AtomBlockId current_atom = q.front();
            q.pop();

            atom_noc_grp_id[current_atom] = noc_grp_id;

            for (auto pin : atom_netlist.block_pins(current_atom)) {
                AtomNetId net_id = atom_netlist.pin_net(pin);
                size_t net_fanout = atom_netlist.net_sinks(net_id).size();

                if (net_fanout >= high_fanout_threshold) {
                    continue;
                }

                AtomBlockId driver_atom_id = atom_netlist.net_driver_block(net_id);
                if (!atom_visited[driver_atom_id]) {
                    q.push(driver_atom_id);
                    atom_visited[driver_atom_id] = true;
                }

                for (auto sink_pin : atom_netlist.net_sinks(net_id)) {
                    AtomBlockId sink_atom_id = atom_netlist.pin_block(sink_pin);
                    if (!atom_visited[sink_atom_id]) {
                        q.push(sink_atom_id);
                        atom_visited[sink_atom_id] = true;
                    }
                }
            }
        }
    }
}
