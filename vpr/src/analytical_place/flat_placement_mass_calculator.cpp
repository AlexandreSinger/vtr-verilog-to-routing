/**
 * @file
 * @author  Alex Singer
 * @date    February 2024
 * @brief   Implementation of the mass calculator used in the AP flow.
 */

#include "flat_placement_mass_calculator.h"
#include <cstring>
#include <queue>
#include <vector>
#include "ap_mass_report.h"
#include "ap_netlist.h"
#include "atom_netlist.h"
#include "atom_netlist_fwd.h"
#include "logic_types.h"
#include "physical_types.h"
#include "prepack.h"
#include "primitive_dim_manager.h"
#include "primitive_vector.h"
#include "vtr_log.h"
#include "vtr_vector.h"

static bool is_primitive_memory_pb_type(const t_pb_type* pb_type) {
    VTR_ASSERT_SAFE(pb_type != nullptr);

    if (pb_type->num_modes != 0)
        return false;

    if (pb_type->class_type != MEMORY_CLASS)
        return false;

    return true;
}

static float calc_pb_type_cost(const t_pb_type* pb_type) {
    if (is_primitive_memory_pb_type(pb_type)) {
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

static float get_atom_mass(AtomBlockId blk_id, const Prepacker& prepacker, const AtomNetlist& atom_netlist) {

    // Get the mass of the most likely pb_type which this atom will be packed into.
    const t_pb_graph_node* primitive = prepacker.get_expected_lowest_cost_pb_gnode(blk_id);
    float mass = calc_pb_type_cost(primitive->pb_type);

    // Since this atom may be part of a molecule, we can use some information
    // about the rest of the molecule to create a more accurate cost.
    // If this primitive is not a memory, then the cost is equal to the number
    // of pins used by the primitive. Some of these pins will be absorbed fully
    // within the cluster due to being part of a molecule. The cost should be
    // adjusted to account for this.
    if (!is_primitive_memory_pb_type(primitive->pb_type)) {
        // Go through each of the atoms pins and check if any of them are fully
        // absorbed into the molecule. These will not contribute to the cost.
        PackMoleculeId mol_id = prepacker.get_atom_molecule(blk_id);
        unsigned absorbed_pins = 0;
        for (AtomPinId atom_pin_id : atom_netlist.block_pins(blk_id)) {
            AtomNetId atom_net_id = atom_netlist.pin_net(atom_pin_id);
            bool all_pins_absorbed = true;
            for (AtomPinId atom_net_pin_id : atom_netlist.net_pins(atom_net_id)) {
                AtomBlockId pin_blk_id = atom_netlist.pin_block(atom_net_pin_id);
                if (prepacker.get_atom_molecule(pin_blk_id) != mol_id) {
                    all_pins_absorbed = false;
                    break;
                }
            }
            if (all_pins_absorbed)
                absorbed_pins++;
        }
        mass -= static_cast<float>(absorbed_pins);
    }

    // Currently, the code does not handle well with fractional masses. Rounding
    // up to the nearest whole number.
    mass = std::ceil(mass);

    return mass;
}

// This method is being forward-declared due to the double recursion below.
// Eventually this should be made into a non-recursive algorithm for performance,
// however this is not in a performance critical part of the code.
static PrimitiveVector calc_pb_type_capacity(const t_pb_type* pb_type,
                                             std::set<PrimitiveVectorDim>& memory_model_dims,
                                             const PrimitiveDimManager& dim_manager);

/**
 * @brief Get the amount of primitives this mode can contain.
 *
 * This is part of a double recursion, since a mode contains primitives which
 * themselves have modes.
 */
static PrimitiveVector calc_mode_capacity(const t_mode& mode,
                                          std::set<PrimitiveVectorDim>& memory_model_dims,
                                          const PrimitiveDimManager& dim_manager) {
    // Accumulate the capacities of all the pbs in this mode.
    PrimitiveVector capacity;
    for (int pb_child_idx = 0; pb_child_idx < mode.num_pb_type_children; pb_child_idx++) {
        const t_pb_type& pb_type = mode.pb_type_children[pb_child_idx];
        PrimitiveVector pb_capacity = calc_pb_type_capacity(&pb_type, memory_model_dims, dim_manager);
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
static PrimitiveVector calc_pb_type_capacity(const t_pb_type* pb_type,
                                             std::set<PrimitiveVectorDim>& memory_model_dims,
                                             const PrimitiveDimManager& dim_manager) {
    // Since a pb cannot be multiple modes at the same time, we do not
    // accumulate the capacities of the mode. Instead we need to "mix" the two
    // capacities as if the pb could choose either one.
    PrimitiveVector capacity;
    // If this is a leaf / primitive, create the base PrimitiveVector capacity.
    if (pb_type->num_modes == 0) {
        LogicalModelId model_id = pb_type->model_id;
        VTR_ASSERT_SAFE(model_id.is_valid());
        PrimitiveVectorDim dim = dim_manager.get_model_dim(model_id);
        VTR_ASSERT(dim.is_valid());
        if (pb_type->class_type == MEMORY_CLASS)
            memory_model_dims.insert(dim);
        float pb_type_cost = calc_pb_type_cost(pb_type);
        capacity.add_val_to_dim(pb_type_cost * pb_type->num_pb, dim);
        return capacity;
    }
    // For now, we simply mix the capacities of modes by taking the max of each
    // dimension of the capcities. This provides an upper-bound on the amount of
    // primitives this pb can contain.
    for (int mode = 0; mode < pb_type->num_modes; mode++) {
        PrimitiveVector mode_capacity = calc_mode_capacity(pb_type->modes[mode], memory_model_dims, dim_manager);
        capacity = PrimitiveVector::max(capacity, mode_capacity);
    }

    // The current pb only has a set number of pins. Therefore, each dimension
    // of the primitive pb cannot have a higher number of pins than this.
    // Clamp the capacity by the score of this pb.
    // FIXME: Handle shared pins
    // FIXME: Should use the pb_score method.
    float total_curr_score = pb_type->num_input_pins + pb_type->num_output_pins + pb_type->num_clock_pins;
    for (PrimitiveVectorDim dim : capacity.get_non_zero_dims()) {
        // If this dimension corresponds to a memory, do not clamp the capacity.
        // Memories just count bits.
        if (memory_model_dims.count(dim) != 0)
            continue;
        capacity.set_dim_val(dim, std::min(capacity.get_dim_val(dim), total_curr_score));
    }

    capacity *= pb_type->num_pb;

    return capacity;
}

/**
 * @brief Calculate the cpacity of the given logical block type.
 */
static PrimitiveVector calc_logical_block_type_capacity(const t_logical_block_type& logical_block_type,
                                                        const PrimitiveDimManager& dim_manager) {
    // If this logical block is empty, it cannot contain any primitives.
    if (logical_block_type.is_empty())
        return PrimitiveVector();

    std::set<PrimitiveVectorDim> memory_model_dims;
    PrimitiveVector capacity = calc_pb_type_capacity(logical_block_type.pb_type,
                                                     memory_model_dims,
                                                     dim_manager);

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
                                       const PrimitiveDimManager& dim_manager,
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

        // Get the dimension in the vector to add value to.
        LogicalModelId model_id = atom_netlist.block_model(atom_blk_id);
        VTR_ASSERT_SAFE(model_id.is_valid());
        PrimitiveVectorDim dim = dim_manager.get_model_dim(model_id);
        VTR_ASSERT(dim.is_valid());

        // Get the amount of mass in this dimension to add.
        float atom_mass = get_atom_mass(atom_blk_id, prepacker, atom_netlist);

        // Add mass to the dimension.
        mass.add_val_to_dim(atom_mass, dim);
    }
    return mass;
}

/**
 * @brief Initialize the dim manager such that every model in the architecture
 *        has a valid dimension in the primitive vector.
 */
static void initialize_dim_manager(PrimitiveDimManager& dim_manager,
                                   const LogicalModels& models,
                                   const std::vector<t_logical_block_type>& logical_block_types,
                                   const AtomNetlist& atom_netlist) {
    // Set the mapping between model IDs and Primitive Vector IDs

    // Pattern-match for "one-hot" primitive modes (i.e. shared primitive dimensions).

    // Populate a lookup between each model and the primitives that implement them.
    // FIXME: Outline to another function.
    VTR_LOG("Looking up the model primitives...\n");
    vtr::vector<LogicalModelId, std::set<t_pb_type*>> model_primitives(models.all_models().size());
    std::queue<t_pb_type*> pb_type_queue;
    for (const t_logical_block_type& block_type : logical_block_types) {
        pb_type_queue.push(block_type.pb_type);
    }
    while (!pb_type_queue.empty()) {
        t_pb_type* pb_type = pb_type_queue.front();
        pb_type_queue.pop();
        if (pb_type == nullptr)
            continue;
        if (pb_type->is_primitive()) {
            model_primitives[pb_type->model_id].insert(pb_type);
            continue;
        }

        for (int mode_idx = 0; mode_idx < pb_type->num_modes; mode_idx++) {
            const t_mode& mode = pb_type->modes[mode_idx];
            for (int pb_child_idx = 0; pb_child_idx < mode.num_pb_type_children; pb_child_idx++) {
                pb_type_queue.push(&mode.pb_type_children[pb_child_idx]);
            }
        }
    }

    // Search for the one-hot primitives.
    VTR_LOG("Searching for one-hot primitives...\n");
    // FIXME: Outline to another function.
    struct OneHotPbType {
        t_pb_type* pb_type;
        std::set<LogicalModelId> shared_models;
    };
    std::vector<OneHotPbType> shared_models;
    VTR_ASSERT(pb_type_queue.empty());
    for (const t_logical_block_type& block_type : logical_block_types) {
        pb_type_queue.push(block_type.pb_type);
    }
    while (!pb_type_queue.empty()) {
        t_pb_type* pb_type = pb_type_queue.front();
        pb_type_queue.pop();
        if (pb_type == nullptr)
            continue;
        if (pb_type->is_primitive()) {
            continue;
        }

        if (pb_type->num_modes > 1) {
            bool is_one_hot = true;
            std::set<LogicalModelId> contained_models;
            // Check if every pb_type in each mode is single and unique.
            for (int mode_idx = 0; mode_idx < pb_type->num_modes; mode_idx++) {
                const t_mode& mode = pb_type->modes[mode_idx];
                if (mode.num_pb_type_children != 1) {
                    is_one_hot = false;
                    break;
                }
                t_pb_type* mode_child_pb = &mode.pb_type_children[0];
                if (mode_child_pb->num_pb > 1 || !mode_child_pb->is_primitive()) {
                    is_one_hot = false;
                    break;
                }
                contained_models.insert(mode_child_pb->model_id);
            }

            if (is_one_hot) {
                OneHotPbType one_hot_pb_type_info;
                one_hot_pb_type_info.pb_type = pb_type;
                one_hot_pb_type_info.shared_models = std::move(contained_models);
                shared_models.push_back(std::move(one_hot_pb_type_info));

                // Do not explore the children of this pb_type.
                continue;
            }
        }
    
        for (int mode_idx = 0; mode_idx < pb_type->num_modes; mode_idx++) {
            const t_mode& mode = pb_type->modes[mode_idx];
            for (int pb_child_idx = 0; pb_child_idx < mode.num_pb_type_children; pb_child_idx++) {
                pb_type_queue.push(&mode.pb_type_children[pb_child_idx]);
            }
        }
    }

    // Log the shared models.
    VTR_LOG("SHARED MODELS:\n");
    for (const OneHotPbType& one_hot_pb_type_info : shared_models) {
        VTR_LOG("\t%s:\n", one_hot_pb_type_info.pb_type->name);
        for (LogicalModelId model_id : one_hot_pb_type_info.shared_models) {
            VTR_LOG("\t\t%s\n", models.model_name(model_id).c_str());
        }
    }

    vtr::vector<LogicalModelId, bool> is_shared_model(models.all_models().size(), false);
    for (const OneHotPbType& one_hot_pb_type_info : shared_models) {
        for (LogicalModelId model_id : one_hot_pb_type_info.shared_models) {
            is_shared_model[model_id] = true;
        }
    }


    // Count the number of occurences of each model in the netlist.
    vtr::vector<LogicalModelId, unsigned> num_model_occurence(models.all_models().size(), 0);
    for (AtomBlockId blk_id : atom_netlist.blocks()) {
        num_model_occurence[atom_netlist.block_model(blk_id)]++;
    }

    // Create a list of models, sorted by their frequency in the netlist.
    // By sorting by frequency, we make the early dimensions more common,
    // which can reduce the overall size of the sparse vector.
    // NOTE: We use stable sort here to keep the order of models the same
    //       as what the user provided in the arch file in the event of a tie.
    std::vector<LogicalModelId> logical_models(models.all_models().begin(), models.all_models().end());
    std::stable_sort(logical_models.begin(), logical_models.end(), [&](LogicalModelId a, LogicalModelId b) {
        // FIXME: Shared models should accumulate the occurences...
        return num_model_occurence[a] > num_model_occurence[b];
    });

    // Create a primitive vector dim for each model.
    for (LogicalModelId model_id : logical_models) {
        if (num_model_occurence[model_id] == 0)
            continue;
        if (is_shared_model[model_id])
            continue;
        dim_manager.create_dim(model_id, models.model_name(model_id));
    }
    for (const OneHotPbType& one_hot_pb_type_info : shared_models) {
        // Create a unique name for the dim. This is used for debugging.
        std::string dim_name = one_hot_pb_type_info.pb_type->name;
        dim_name += "_ap_shared[";
        size_t index = 0;
        for (LogicalModelId model_id : one_hot_pb_type_info.shared_models) {
            dim_name += models.model_name(model_id);
            dim_name += "]";
            if (index != one_hot_pb_type_info.shared_models.size() - 1)
                dim_name += "[";
            index++;
        }
        PrimitiveVectorDim new_dim = dim_manager.create_empty_dim(dim_name);
        for (LogicalModelId model_id : one_hot_pb_type_info.shared_models) {
            VTR_ASSERT(!dim_manager.get_model_dim(model_id).is_valid());
            dim_manager.add_model_to_dim(model_id, new_dim);
        }
    }
    // Now add all unused and non-shared models
    for (LogicalModelId model_id : logical_models) {
        if (num_model_occurence[model_id] != 0)
            continue;
        if (is_shared_model[model_id])
            continue;

        VTR_ASSERT(!dim_manager.get_model_dim(model_id).is_valid());
        dim_manager.create_dim(model_id, models.model_name(model_id));
    }
}

FlatPlacementMassCalculator::FlatPlacementMassCalculator(const APNetlist& ap_netlist,
                                                         const Prepacker& prepacker,
                                                         const AtomNetlist& atom_netlist,
                                                         const std::vector<t_logical_block_type>& logical_block_types,
                                                         const std::vector<t_physical_tile_type>& physical_tile_types,
                                                         const LogicalModels& models,
                                                         int log_verbosity)
    : physical_tile_type_capacity_(physical_tile_types.size())
    , logical_block_type_capacity_(logical_block_types.size())
    , block_mass_(ap_netlist.blocks().size())
    , log_verbosity_(log_verbosity) {

    // Initialize the mapping between model IDs and Primitive Vector dims
    initialize_dim_manager(primitive_dim_manager_,
                           models,
                           logical_block_types,
                           atom_netlist);

    // Precompute the capacity of each logical block type.
    for (const t_logical_block_type& logical_block_type : logical_block_types) {
        logical_block_type_capacity_[logical_block_type.index] = calc_logical_block_type_capacity(logical_block_type, primitive_dim_manager_);
    }

    // Precompute the capacity of each physical tile type.
    for (const t_physical_tile_type& physical_tile_type : physical_tile_types) {
        physical_tile_type_capacity_[physical_tile_type.index] = calc_physical_tile_type_capacity(physical_tile_type, logical_block_type_capacity_);
    }

    // Precompute the mass of each block in the APNetlist
    VTR_LOGV(log_verbosity_ >= 10, "Pre-computing the block masses...\n");
    for (APBlockId ap_block_id : ap_netlist.blocks()) {
        block_mass_[ap_block_id] = calc_block_mass(ap_block_id,
                                                   primitive_dim_manager_,
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
                            primitive_dim_manager_,
                            ap_netlist);
}
