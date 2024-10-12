
#include "global_placer.h"
#include <cstdio>
#include <limits>
#include <memory>
#include "analytical_solver.h"
#include "ap_netlist.h"
#include "partial_legalizer.h"
#include "partial_placement.h"
#include "vpr_error.h"
#include "vtr_log.h"
#include "vtr_time.h"

std::unique_ptr<GlobalPlacer> make_global_placer(e_global_placer placer_type,
                                                 const APNetlist& netlist) {
    // Based on the placer type passed in, build the global placer.
    switch (placer_type) {
        case e_global_placer::SimPL:
            return std::make_unique<SimPLGlobalPlacer>(netlist);
        default:
            VPR_FATAL_ERROR(VPR_ERROR_AP,
                            "Unrecognized global placer type");

    }
}

SimPLGlobalPlacer::SimPLGlobalPlacer(const APNetlist& netlist) : GlobalPlacer(netlist) {
    vtr::ScopedStartFinishTimer global_placer_building_timer("Constructing Global Placer");
    solver_ = make_analytical_solver(e_analytical_solver::QP_HYBRID, netlist);
    partial_legalizer_ = make_partial_legalizer(e_partial_legalizer::FLOW_BASED, netlist);
}

static void print_SimPL_status_header() {
    VTR_LOG("---- ---------------- ---------------- ----------- -------------- ----------\n");
    VTR_LOG("Iter Lower Bound HPWL Upper Bound HPWL Solver Time Legalizer Time Total Time\n");
    VTR_LOG("                                             (sec)          (sec)      (sec)\n");
    VTR_LOG("---- ---------------- ---------------- ----------- -------------- ----------\n");
}

static void print_SimPL_status(size_t iteration,
                               double lb_hpwl,
                               double ub_hpwl,
                               float solver_time,
                               float legalizer_time,
                               float total_time) {
    // Iteration
    VTR_LOG("%4zu", iteration);

    // Lower Bound HPWL
    VTR_LOG(" %16.2f", lb_hpwl);

    // Upper Bound HPWL
    VTR_LOG(" %16.2f", ub_hpwl);

    // Solver runtime
    VTR_LOG(" %11.3f", solver_time);

    // Legalizer runtime
    VTR_LOG(" %14.3f", legalizer_time);

    // Total runtime
    VTR_LOG(" %10.3f", total_time);

    VTR_LOG("\n");

    fflush(stdout);
}

PartialPlacement SimPLGlobalPlacer::place() {
    vtr::ScopedStartFinishTimer global_placer_time("AP Global Placer");
    PartialPlacement p_placement(netlist_);

    vtr::Timer runtime_timer;

    print_SimPL_status_header();
    for (int i = 0; i < 100; i++) {
        float iter_start_time = runtime_timer.elapsed_sec();
        // Run the solver.
        float solver_start_time = runtime_timer.elapsed_sec();
        solver_->solve(i, p_placement);
        float solver_end_time = runtime_timer.elapsed_sec();
        double lb_hpwl = p_placement.get_hpwl(netlist_);
        // Run the legalizer.
        float legalizer_start_time = runtime_timer.elapsed_sec();
        partial_legalizer_->legalize(p_placement);
        float legalizer_end_time = runtime_timer.elapsed_sec();
        double ub_hpwl = p_placement.get_hpwl(netlist_);
        // Print some stats
        float iter_end_time = runtime_timer.elapsed_sec();
        print_SimPL_status(i, lb_hpwl, ub_hpwl,
                           solver_end_time - solver_start_time,
                           legalizer_end_time - legalizer_start_time,
                           iter_end_time - iter_start_time);
        if (ub_hpwl - lb_hpwl < 100.0)
            break;
    }

    return p_placement;
}

