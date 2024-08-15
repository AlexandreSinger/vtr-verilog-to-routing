
#include "analytical_placement_flow.h"
#include <memory>
#include "PartialPlacement.h"
#include "ap_netlist.h"
#include "ap_utils.h"
#include "AnalyticalSolver.h"
#include "PlacementLegalizer.h"
#include "atom_netlist.h"
#include "full_legalizer.h"
#include "globals.h"
#include "prepack.h"
#include "read_atom_netlist.h"
#include "user_place_constraints.h"
#include "vpr_context.h"
#include "vpr_types.h"
#include "vtr_assert.h"
#include "vtr_time.h"

void run_analytical_placement_flow(t_vpr_setup& vpr_setup) {
    vtr::ScopedStartFinishTimer timer("Analytical Placement Flow");

    // The global state used/modified by this flow.
    AtomContext& mutable_atom_ctx = g_vpr_ctx.mutable_atom();
    const AtomNetlist& atom_nlist = g_vpr_ctx.atom().nlist;
    const DeviceContext& device_ctx = g_vpr_ctx.device();
    const UserPlaceConstraints& constraints = g_vpr_ctx.floorplanning().constraints;

    // Run the prepacker
    Prepacker prepacker;
    prepacker.init(atom_nlist, device_ctx.logical_block_types);

    // Create the ap netlist from the atom netlist using the result from the
    // prepacker.
    APNetlist ap_netlist = read_atom_netlist(mutable_atom_ctx, prepacker, constraints);
    print_ap_netlist_stats(ap_netlist);

    // Set up the partial placement object
    PartialPlacement p_placement(ap_netlist);
    // Create the solver
    std::unique_ptr<AnalyticalSolver> solver = make_analytical_solver(e_analytical_solver::QP_HYBRID, ap_netlist);
    // Create the legalizer
    FlowBasedLegalizer legalizer(ap_netlist);
    // This for loop always starts at iteration 0
    for (unsigned iteration = 0; iteration < 100; iteration++) {
        VTR_LOG("iteration: %ld\n", iteration);
        solver->solve(iteration, p_placement);
        VTR_ASSERT(p_placement.verify(ap_netlist, device_ctx) && "placement not valid after solve!");
        double post_solve_hpwl = get_hpwl(p_placement, ap_netlist);
        VTR_LOG("HPWL: %f\n", post_solve_hpwl);
        // Partial legalization using flow-based algorithm
        legalizer.legalize(p_placement);
        VTR_ASSERT(p_placement.verify(ap_netlist, device_ctx) && "placement not valid after legalize!");
        double post_legalize_hpwl = get_hpwl(p_placement, ap_netlist);
        VTR_LOG("Post-Legalized HPWL: %f\n", post_legalize_hpwl);
        if(std::abs(post_solve_hpwl - post_legalize_hpwl) < 20){
            VTR_LOG("ended because of convergence\n");
            break;
        }
        // p_placement.unicode_art();
    }

    // Export to a flat placement file.
    export_to_flat_placement_file(p_placement, ap_netlist, mutable_atom_ctx.nlist, "flat_placement_file.txt");

    // Run the full legalizer
    FullLegalizer full_legalizer(ap_netlist,
                                 vpr_setup,
                                 g_vpr_ctx.device().grid,
                                 g_vpr_ctx.device().arch,
                                 g_vpr_ctx.atom().nlist,
                                 prepacker,
                                 g_vpr_ctx.device().logical_block_types,
                                 vpr_setup.PackerRRGraph,
                                 g_vpr_ctx.device().arch->models,
                                 g_vpr_ctx.device().arch->model_library,
                                 vpr_setup.PackerOpts);
    full_legalizer.legalize(p_placement);
}

