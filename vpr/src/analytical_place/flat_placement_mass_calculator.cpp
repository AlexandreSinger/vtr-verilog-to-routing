/**
 * @file
 * @author  Alex Singer
 * @date    February 2024
 * @brief   Implementation of the mass calculator used in the AP flow.
 */

#include "flat_placement_mass_calculator.h"
#include <cstring>
#include <vector>
#include "ap_mass_report.h"
#include "ap_netlist.h"
#include "arch_types.h"
#include "atom_netlist.h"
#include "atom_netlist_fwd.h"
#include "device_grid.h"
#include "globals.h"
#include "logic_types.h"
#include "physical_types.h"
#include "prepack.h"
#include "primitive_vector.h"
#include "vpr_utils.h"
#include "vtr_log.h"

/**
 * @brief Get the scalar mass of the given model (primitive type).
 *
 * A model with a higher mass will take up more space in its bin which may force
 * more spreading of that type of primitive.
 *
 * TODO: This will be made more complicated later. Models may be weighted based
 *       on some factors.
 */
static float get_model_mass(LogicalModelId model_id) {
    // If this model is a LUT
    if (model_id == LogicalModels::MODEL_NAMES_ID) {
        // FIXME: Remove global reference!
        const LogicalModels& models = g_vpr_ctx.device().arch->models;
        // Go through each of the ports (I think there should only be 1) and count
        // the max number of pins.
        int lut_size = 0;
        t_model_ports* curr_input = models.get_model(model_id).inputs;
        while (curr_input != nullptr) {
            lut_size += curr_input->size;
            curr_input = curr_input->next;
        }

        // FIXME: Is this possible?
        if (lut_size == 0)
            return 1.f;

        return 1.0f * lut_size;
    }
    // Currently, all models have a mass of one.
    (void)model_id;
    return 1.f;
}

static float calc_pb_type_cost(const t_pb_type* pb_type) {
    if (pb_type->num_modes == 0 && pb_type->class_type == MEMORY_CLASS) {
        int num_address_pins = 0;
        int num_data_out_pins = 0;
        int num_address1_pins = 0;
        int num_data_out1_pins = 0;
        int num_address2_pins = 0;
        int num_data_out2_pins = 0;
        for (int i = 0; i < pb_type->num_ports; i++) {
            if (pb_type->ports[i].port_class == nullptr)
                continue;

            std::string port_class = pb_type->ports[i].port_class;
            if (port_class == "address")
                num_address_pins += pb_type->ports[i].num_pins;
            if (port_class == "data_out")
                num_data_out_pins += pb_type->ports[i].num_pins;
            if (port_class == "address1")
                num_address1_pins += pb_type->ports[i].num_pins;
            if (port_class == "data_out1")
                num_data_out1_pins += pb_type->ports[i].num_pins;
            if (port_class == "address2")
                num_address2_pins += pb_type->ports[i].num_pins;
            if (port_class == "data_out2")
                num_data_out2_pins += pb_type->ports[i].num_pins;
        }

        int single_port_num_bits = 0;
        if (num_address_pins > 0) {
            VTR_ASSERT_MSG(num_address1_pins == 0 && num_address2_pins == 0,
                           "Cannot be a single port and dual port memory");
            single_port_num_bits = (1 << num_address_pins) * num_data_out_pins;
        }
        int dual_port1_num_bits = 0;
        if (num_address1_pins > 0) {
            // FIXME: Found that the Titan architecture fails this test! Raise
            //        an issue on this to maybe update the documentation or fix
            //        the architecture.
            //        See: prim_ram_1Kx9
            //        Perhaps the architecture parsing file should check for this.
            /*
            VTR_ASSERT_MSG(num_address_pins == 0 && num_address2_pins > 0,
                           "Ill-formed dual port memory");
            */
            dual_port1_num_bits = (1 << num_address1_pins) * num_data_out1_pins;
        }
        int dual_port2_num_bits = 0;
        if (num_address2_pins > 0) {
            VTR_ASSERT_MSG(num_address_pins == 0 && num_address1_pins > 0,
                           "Ill-formed dual port memory");
            dual_port2_num_bits = (1 << num_address2_pins) * num_data_out2_pins;
        }
        // Note: We take the max of the two dual port num bits since dual ports
        //       access the same memory but with different address bits.
        int total_num_bits = single_port_num_bits + std::max(dual_port1_num_bits, dual_port2_num_bits);
        VTR_ASSERT_MSG(total_num_bits != 0,
                       "Cannot deduce number of bits in memory primitive.");
        return total_num_bits;
    }
    float pb_cost = pb_type->num_input_pins + pb_type->num_output_pins + pb_type->num_clock_pins;
    return pb_cost;
}

// FIXME: This should really take in a molecule. If we knew the pack pattern,
//        we could recognize carry chains.
static float get_atom_mass(AtomBlockId blk_id, const Prepacker& prepacker, const AtomNetlist& atom_netlist) {
    const t_pb_graph_node* primitive = prepacker.get_expected_lowest_cost_pb_gnode(blk_id);

    return calc_pb_type_cost(primitive->pb_type);
    // return compute_primitive_base_cost(primitive);

    LogicalModelId blk_model = atom_netlist.block_model(blk_id);
    // If this model is a LUT
    if (blk_model == LogicalModels::MODEL_NAMES_ID) {
        // NOTE: See base/read_circuit.cpp
        auto in_ports = atom_netlist.block_input_ports(blk_id);
        VTR_ASSERT_MSG(in_ports.size() <= 1,
                       "Expected number of input ports for LUT to be 0 or 1");
        int lut_size = 0;
        if (in_ports.size() == 1) {
            // Use the number of pins in the input port to determine the
            // size of the LUT.
            auto port_id = *in_ports.begin();
            lut_size = atom_netlist.port_pins(port_id).size();
        }

        // FIXME: This is wicked dangerous! Need to check if inputs exists.
        // lut_size = std::max(lut_size, blk_model->inputs->min_size);
        lut_size = std::max(lut_size, 5);
        // FIXME: I think we need to clamp the lut size between the min and max
        //        number of inputs that its model can have...

        // If this is a zero LUT, just pretend that it is a 1-LUT.
        if (lut_size == 0)
            return 1.f;

        return lut_size;
    }

    // Anything we do not recognize, just set to 1.
    return 1.f;
}

// This method is being forward-declared due to the double recursion below.
// Eventually this should be made into a non-recursive algorithm for performance,
// however this is not in a performance critical part of the code.
static PrimitiveVector calc_pb_type_capacity(const t_pb_type* pb_type, std::set<int>& memory_model_ids);

/**
 * @brief Get the amount of primitives this mode can contain.
 *
 * This is part of a double recursion, since a mode contains primitives which
 * themselves have modes.
 */
static PrimitiveVector calc_mode_capacity(const t_mode& mode, std::set<int>& memory_model_ids) {
    // Accumulate the capacities of all the pbs in this mode.
    PrimitiveVector capacity;
    for (int pb_child_idx = 0; pb_child_idx < mode.num_pb_type_children; pb_child_idx++) {
        const t_pb_type& pb_type = mode.pb_type_children[pb_child_idx];
        PrimitiveVector pb_capacity = calc_pb_type_capacity(&pb_type, memory_model_ids);
        // FIXME: This was moved.
        // A mode may contain multiple pbs of the same type, multiply the
        // capacity.
        // pb_capacity *= pb_type.num_pb;
        capacity += pb_capacity;
    }
    return capacity;
}

/**
 * @brief Get the amount of primitives this pb can contain.
 *
 * This is the other part of the double recursion. A pb may have multiple modes.
 * Modes are made of pbs.
 */
static PrimitiveVector calc_pb_type_capacity(const t_pb_type* pb_type, std::set<int>& memory_model_ids) {
    // Since a pb cannot be multiple modes at the same time, we do not
    // accumulate the capacities of the mode. Instead we need to "mix" the two
    // capacities as if the pb could choose either one.
    PrimitiveVector capacity;
    // If this is a leaf / primitive, create the base PrimitiveVector capacity.
    if (pb_type->num_modes == 0) {
        LogicalModelId model_id = pb_type->model_id;
        VTR_ASSERT_SAFE(model_id.is_valid());
        if (pb_type->class_type == MEMORY_CLASS)
            memory_model_ids.insert((size_t)model_id);
        float pb_type_cost = calc_pb_type_cost(pb_type);
        capacity.add_val_to_dim(pb_type_cost * pb_type->num_pb, (size_t)model_id);
        return capacity;
    }
    // For now, we simply mix the capacities of modes by taking the max of each
    // dimension of the capcities. This provides an upper-bound on the amount of
    // primitives this pb can contain.
    for (int mode = 0; mode < pb_type->num_modes; mode++) {
        PrimitiveVector mode_capacity = calc_mode_capacity(pb_type->modes[mode], memory_model_ids);
        capacity = PrimitiveVector::max(capacity, mode_capacity);
    }

    // The current pb only has a set number of pins. Therefore, each dimension
    // of the primitive pb cannot have a higher number of pins than this.
    // Clamp the capacity by the score of this pb.
    // FIXME: Handle shared pins
    // FIXME: Should use the pb_score method.
    float total_curr_score = pb_type->num_input_pins + pb_type->num_output_pins + pb_type->num_clock_pins;
    for (int dim : capacity.get_non_zero_dims()) {
        // If this dimension corresponds to a memory, do not clamp the capacity.
        // Memories just count bits.
        if (memory_model_ids.count(dim) != 0)
            continue;
        capacity.set_dim_val(dim, std::min(capacity.get_dim_val(dim), total_curr_score));
    }

    capacity *= pb_type->num_pb;
    // This is the pin-sharing scaling factor. If there are multiple pbs in the
    // same pb_type, they are likely to share inputs.
    // FIXME: We may be able to actually calculate this number based on the
    //        architecture file.
    /*
    if (pb_type->num_pb > 1)
        capacity *= 0.75f;
        */

    // FIXME: get the number above by taking the number of pins of the parent
    //        pb / (num_pins in this pb * num_pb)
    /*
    if (pb_type->parent_mode != nullptr && pb_type->parent_mode->parent_pb_type != nullptr) {
        t_pb_type* parent_pb_type = pb_type->parent_mode->parent_pb_type;

        float curr_pb_score = pb_type->num_input_pins + pb_type->num_output_pins + pb_type->num_clock_pins;
        curr_pb_score *= pb_type->num_pb;

        float parent_pb_score = parent_pb_type->num_input_pins + parent_pb_type->num_output_pins + parent_pb_type->num_clock_pins;
        // FIXME: Need to handle siblings. If there are other children of the parent,
        //        they too will need pins.
        //        I think putting this calculation in the parent instead of the child
        //        is a much much better idea.
    }
    */

    return capacity;
}

/**
 * @brief Calculate the cpacity of the given logical block type.
 */
static PrimitiveVector calc_logical_block_type_capacity(const t_logical_block_type& logical_block_type) {
    // If this logical block is empty, it cannot contain any primitives.
    if (logical_block_type.is_empty())
        return PrimitiveVector();

    std::set<int> memory_model_ids;
    PrimitiveVector capacity = calc_pb_type_capacity(logical_block_type.pb_type, memory_model_ids);

    // The current pb only has a set number of pins. Therefore, each dimension
    // of the primitive pb cannot have a higher number of pins than this.
    // Clamp the capacity by the score of this pb.
    // FIXME: Handle shared pins
    // t_pb_type* pb_type = logical_block_type.pb_type;
    // float total_curr_score = pb_type->num_input_pins + pb_type->num_output_pins + pb_type->num_clock_pins;
    // for (int dim : capacity.get_non_zero_dims()) {
    //     capacity.set_dim_val(dim, std::min(capacity.get_dim_val(dim), total_curr_score));
    // }

    // The primitive capacity of a logical block is the primitive capacity of
    // its root pb.
    return capacity;
}

/**
 * @brief Get the primitive capacity of the given sub_tile.
 *
 * Sub_tiles may reuse logical blocks between one another, therefore this method
 * requires that the capacities of all of the logical blocks have been
 * pre-calculated and stored in the given vector.
 *
 *  @param sub_tile                         The sub_tile to get the capacity of.
 *  @param logical_block_type_capacities    The capacities of all logical block
 *                                          types.
 */
static PrimitiveVector calc_sub_tile_capacity(const t_sub_tile& sub_tile,
                                              const std::vector<PrimitiveVector>& logical_block_type_capacities) {
    // Similar to getting the primitive capacity of the pb, sub_tiles have many
    // equivalent sites, but it can only be one of them at a time. Need to "mix"
    // the capacities of the different sites this sub_tile may be.
    PrimitiveVector capacity;
    for (t_logical_block_type_ptr block_type : sub_tile.equivalent_sites) {
        const PrimitiveVector& block_capacity = logical_block_type_capacities[block_type->index];
        // Currently, we take the max of each primitive dimension as an upper
        // bound on the capacity of the sub_tile.
        capacity = PrimitiveVector::max(capacity, block_capacity);
    }
    return capacity;
}

/**
 * @brief Get the primitive capacity of a tile of the given type.
 *
 * Tiles may reuse logical blocks between one another, therefore this method
 * requires that the capacities of all of the logical blocks have been
 * pre-calculated and stored in the given vector.
 *
 *  @param tile_type                        The tile type to get the capacity of.
 *  @param logical_block_type_capacities    The capacities of all logical block
 *                                          types.
 */
static PrimitiveVector calc_physical_tile_type_capacity(const t_physical_tile_type& tile_type,
                                                        const std::vector<PrimitiveVector>& logical_block_type_capacities) {
    // Accumulate the capacities of all the sub_tiles in the given tile type.
    PrimitiveVector capacity;
    for (const t_sub_tile& sub_tile : tile_type.sub_tiles) {
        PrimitiveVector sub_tile_capacity = calc_sub_tile_capacity(sub_tile, logical_block_type_capacities);
        // A tile may contain many sub_tiles of the same type. Multiply by the
        // number of sub_tiles of this type.
        sub_tile_capacity *= sub_tile.capacity.total();
        capacity += sub_tile_capacity;
    }
    return capacity;
}

/**
 * @brief Get the primitive mass of the given block.
 *
 * This returns an M-dimensional vector with each entry indicating the mass of
 * that primitive type in this block. M is the number of unique models
 * (primitive types) in the architecture.
 */
static PrimitiveVector calc_block_mass(APBlockId blk_id,
                                       const APNetlist& netlist,
                                       const Prepacker& prepacker,
                                       const AtomNetlist& atom_netlist) {
    PrimitiveVector mass;
    PackMoleculeId mol_id = netlist.block_molecule(blk_id);
    const t_pack_molecule& mol = prepacker.get_molecule(mol_id);
    for (AtomBlockId atom_blk_id : mol.atom_block_ids) {
        // See issue #2791, some of the atom_block_ids may be invalid. They can
        // safely be ignored.
        if (!atom_blk_id.is_valid())
            continue;
        LogicalModelId model_id = atom_netlist.block_model(atom_blk_id);
        VTR_ASSERT_SAFE(model_id.is_valid());
        // FIXME: Get atom mass should not return a float.
        float atom_mass = get_atom_mass(atom_blk_id, prepacker, atom_netlist);
        // Go through each of the atoms pins and check if any of them are fully
        // absorbed into the molecule. These will not contribute to the cost.
        // FIXME: This can be computed better.
        // float absorbed_pins = 0.0f;
        // for (AtomPinId atom_pin_id : atom_netlist.block_pins(atom_blk_id)) {
        //     AtomNetId atom_net_id = atom_netlist.pin_net(atom_pin_id);
        //     bool all_pins_absorbed = true;
        //     for (AtomPinId atom_net_pin_id : atom_netlist.net_pins(atom_net_id)) {
        //         AtomBlockId pin_blk_id = atom_netlist.pin_block(atom_net_pin_id);
        //         if (prepacker.get_atom_molecule(pin_blk_id) != mol_id) {
        //             all_pins_absorbed = false;
        //             break;
        //         }
        //     }
        //     if (all_pins_absorbed)
        //         absorbed_pins++;
        // }
        // float total_cost = atom_mass - absorbed_pins;
        // FIXME: HACK! Found that for VTR architecture, just dividing by 6
        //        puts the values in range. Should instead actually count how
        //        much data is stored within the memory. This can be stored in
        //        the get_atom_mass class.
        if (prepacker.get_expected_lowest_cost_pb_gnode(atom_blk_id)->pb_type->class_type == MEMORY_CLASS) {
            // atom_mass /= 6.0f;
        }
        // FIXME: Currently fractional costs cause problems due to numerical
        //        instability. Round up to nearest whole number. We round up here
        //        since we never want the cost to become zero.
        float total_cost = std::ceil(atom_mass);
        mass.add_val_to_dim(total_cost, (size_t)model_id);
    }
    return mass;
}

FlatPlacementMassCalculator::FlatPlacementMassCalculator(const APNetlist& ap_netlist,
                                                         const Prepacker& prepacker,
                                                         const AtomNetlist& atom_netlist,
                                                         const std::vector<t_logical_block_type>& logical_block_types,
                                                         const std::vector<t_physical_tile_type>& physical_tile_types,
                                                         int log_verbosity)
    : physical_tile_type_capacity_(physical_tile_types.size())
    , logical_block_type_capacity_(logical_block_types.size())
    , block_mass_(ap_netlist.blocks().size())
    , log_verbosity_(log_verbosity) {

    // Precompute the capacity of each logical block type.
    for (const t_logical_block_type& logical_block_type : logical_block_types) {
        logical_block_type_capacity_[logical_block_type.index] = calc_logical_block_type_capacity(logical_block_type);
    }

    // Precompute the capacity of each physical tile type.
    for (const t_physical_tile_type& physical_tile_type : physical_tile_types) {
        physical_tile_type_capacity_[physical_tile_type.index] = calc_physical_tile_type_capacity(physical_tile_type, logical_block_type_capacity_);
    }

    // Precompute the mass of each block in the APNetlist
    VTR_LOGV(log_verbosity_ >= 10, "Pre-computing the block masses...\n");
    for (APBlockId ap_block_id : ap_netlist.blocks()) {
        block_mass_[ap_block_id] = calc_block_mass(ap_block_id,
                                                   ap_netlist,
                                                   prepacker,
                                                   atom_netlist);
    }
    VTR_LOGV(log_verbosity_ >= 10, "Finished pre-computing the block masses.\n");
}

void FlatPlacementMassCalculator::generate_mass_report(const APNetlist& ap_netlist) const {
    generate_ap_mass_report(logical_block_type_capacity_,
                            physical_tile_type_capacity_,
                            block_mass_,
                            ap_netlist);
}
