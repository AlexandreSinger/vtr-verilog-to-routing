#pragma once

#include <vector>
#include "read_options.h"
#include "physical_types.h"
#include "vpr_types.h"

void SetupVPR(const t_options* Options,
              const bool TimingEnabled,
              const bool readArchFile,
              t_file_name_opts* FileNameOpts,
              t_arch* Arch,
              t_netlist_opts* NetlistOpts,
              t_packer_opts* PackerOpts,
              t_placer_opts* PlacerOpts,
              t_ap_opts* APOpts,
              t_router_opts* RouterOpts,
              t_analysis_opts* AnalysisOpts,
              t_noc_opts* NocOpts,
              t_server_opts* ServerOpts,
              t_det_routing_arch& RoutingArch,
              std::vector<t_lb_type_rr_node>** PackerRRGraphs,
              std::vector<t_segment_inf>& Segments,
              t_timing_inf* Timing,
              bool* ShowGraphics,
              int* GraphPause,
              bool* SaveGraphics,
              std::string* GraphicsCommands,
              t_power_opts* PowerOpts,
              t_vpr_setup* vpr_setup);
