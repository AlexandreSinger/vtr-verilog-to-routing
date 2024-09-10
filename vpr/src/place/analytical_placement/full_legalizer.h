
#pragma once

#include <vector>

class APNetlist;
class AtomNetlist;
class DeviceGrid;
class PartialPlacement;
class Prepacker;
struct t_logical_block_type;
struct t_lb_type_rr_node;
struct t_model;
struct t_packer_opts;
struct t_vpr_setup;
struct t_arch;

class FullLegalizer {
public:
    // Bring in all the necessary state here. This is the state needed from the
    // AP Context. the Packer Context, and the Placer Context.
    FullLegalizer(const APNetlist& ap_netlist,
                  t_vpr_setup& vpr_setup,
                  const DeviceGrid& device_grid,
                  const t_arch* arch,
                  const AtomNetlist& atom_netlist,
                  const Prepacker& prepacker,
                  const std::vector<t_logical_block_type>& logical_block_types,
                  std::vector<t_lb_type_rr_node>* lb_type_rr_graphs,
                  const t_model* user_models,
                  const t_model* library_models,
                  const t_packer_opts& packer_opts)
            : ap_netlist_(ap_netlist),
              vpr_setup_(vpr_setup),
              device_grid_(device_grid),
              arch_(arch),
              atom_netlist_(atom_netlist),
              prepacker_(prepacker),
              logical_block_types_(logical_block_types),
              lb_type_rr_graphs_(lb_type_rr_graphs),
              user_models_(user_models),
              library_models_(library_models),
              packer_opts_(packer_opts) {}

    void legalize(PartialPlacement& p_placement);
private:
    // AP Context Info
    const APNetlist& ap_netlist_;
    // Overall Setup Info
    // FIXME: I do not like bringing all of this in. Perhaps clean up the methods
    //        that use it.
    t_vpr_setup& vpr_setup_;
    // Device Context Info
    const DeviceGrid& device_grid_;
    const t_arch* arch_;
    // Packing Context Info
    const AtomNetlist& atom_netlist_;
    const Prepacker& prepacker_;
    const std::vector<t_logical_block_type>& logical_block_types_;
    std::vector<t_lb_type_rr_node>* lb_type_rr_graphs_;
    const t_model* user_models_;
    const t_model* library_models_;
    const t_packer_opts& packer_opts_;
    // Placement Context Info
};

