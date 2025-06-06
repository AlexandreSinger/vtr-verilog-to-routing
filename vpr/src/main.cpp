/**
 * VPR is a CAD tool used to conduct FPGA architecture exploration.  It takes, as input, a technology-mapped netlist and a description of the FPGA architecture being investigated.
 * VPR then generates a packed, placed, and routed FPGA (in .net, .place, and .route files respectively) that implements the input netlist.
 *
 * This file is where VPR starts execution.
 *
 * Key files in VPR:
 * 1.  libarchfpga/physical_types.h - Data structures that define the properties of the FPGA architecture
 * 2.  vpr_types.h - Very major file that defines the core data structures used in VPR.  This includes detailed architecture information, user netlist data structures, and data structures that describe the mapping between those two.
 * 3.  globals.h - Defines the global variables used by VPR.
 */

#include <cstdio>
#include <cstring>
#include <ctime>

#include "vtr_error.h"
#include "vtr_memory.h"
#include "vtr_log.h"
#include "vtr_time.h"

#include "tatum/error.hpp"

#include "vpr_exit_codes.h"
#include "vpr_error.h"
#include "vpr_api.h"
#include "vpr_signal_handler.h"
#include "vpr_tatum_error.h"

#include "globals.h"

/**
 * VPR program
 * Generate FPGA architecture given architecture description
 * Pack, place, and route circuit into FPGA architecture
 * Electrical timing analysis on results
 *
 * Overall steps
 * 1.  Initialization
 * 2.  Pack
 * 3.  Place-and-route and timing analysis
 * 4.  Clean up
 */
int main(int argc, const char** argv) {
    vtr::ScopedFinishTimer t("The entire flow of VPR");

    t_options Options = t_options();
    t_arch Arch = t_arch();
    t_vpr_setup vpr_setup = t_vpr_setup();

    try {
        vpr_install_signal_handler();

        /* Read options, architecture, and circuit netlist */
        vpr_init(argc, argv, &Options, &vpr_setup, &Arch);
        if (Options.show_version) {
            vpr_free_all(Arch, vpr_setup);
            return SUCCESS_EXIT_CODE;
        }

        if (Options.show_arch_resources) {
            vpr_print_arch_resources(vpr_setup, Arch);
            vpr_free_all(Arch, vpr_setup);
            return SUCCESS_EXIT_CODE;
        }

        bool flow_succeeded = vpr_flow(vpr_setup, Arch);
        if (!flow_succeeded) {
            VTR_LOG("VPR failed to implement circuit\n");
            vpr_free_all(Arch, vpr_setup);
            return UNIMPLEMENTABLE_EXIT_CODE;
        }

        auto& timing_ctx = g_vpr_ctx.timing();
        print_timing_stats("Flow", timing_ctx.stats);

        /* free data structures */
        vpr_free_all(Arch, vpr_setup);

        VTR_LOG("VPR succeeded\n");

    } catch (const tatum::Error& tatum_error) {
        VTR_LOG_ERROR("%s\n", format_tatum_error(tatum_error).c_str());
        vpr_free_all(Arch, vpr_setup);

        return ERROR_EXIT_CODE;

    } catch (const VprError& vpr_error) {
        vpr_print_error(vpr_error);
        vpr_free_all(Arch, vpr_setup);
        if (vpr_error.type() == VPR_ERROR_INTERRUPTED) {
            return INTERRUPTED_EXIT_CODE;
        } else {
            return ERROR_EXIT_CODE;
        }

    } catch (const vtr::VtrError& vtr_error) {
        VTR_LOG_ERROR("%s:%d %s\n", vtr_error.filename_c_str(), vtr_error.line(), vtr_error.what());
        vpr_free_all(Arch, vpr_setup);

        return ERROR_EXIT_CODE;
    }

    /* Signal success to scripts */
    return SUCCESS_EXIT_CODE;
}
