
#include "read_atom_netlist.h"
#include "ap_netlist.h"
#include "atom_netlist.h"
#include "atom_netlist_fwd.h"
#include "netlist_fwd.h"
#include "partition.h"
#include "partition_region.h"
#include "prepack.h"
#include "vpr_context.h"
#include "vpr_types.h"
#include "vtr_assert.h"
#include "vtr_geometry.h"
#include "vtr_time.h"
#include <unordered_set>
#include <vector>

APNetlist read_atom_netlist(const AtomContext& atom_ctx, const Prepacker& prepacker, const UserPlaceConstraints& constraints) {
    vtr::ScopedStartFinishTimer timer("Read Atom Netlist to AP Netlist");

    // FIXME: What to do about the name and ID in this context? For now just
    //        using empty strings.
    APNetlist ap_netlist;

    // Add the APBlocks based on the atom block molecules. This essentially
    // creates supernodes.
    // Each AP block has the name of the first atom block in the molecule.
    // Each port is named "<atom_blk_name>_<atom_port_name>"
    // Each net has the exact same name at the atom netlist
    const AtomNetlist& atom_netlist = atom_ctx.nlist;
    for (AtomBlockId atom_blk_id : atom_netlist.blocks()) {
        // Get the molecule of this block
        t_pack_molecule* mol = prepacker.get_atom_molecule(atom_blk_id);
        // Create the AP block (if not already done)
        const std::string& first_blk_name = atom_netlist.block_name(mol->atom_block_ids[0]);
        APBlockId ap_blk_id = ap_netlist.create_block(first_blk_name, mol);
        // Add the ports and pins of this block to the supernode
        for (AtomPortId atom_port_id : atom_netlist.block_ports(atom_blk_id)) {
            BitIndex port_width = atom_netlist.port_width(atom_port_id);
            PortType port_type = atom_netlist.port_type(atom_port_id);
            const std::string& port_name = atom_netlist.port_name(atom_port_id);
            const std::string& block_name = atom_netlist.block_name(atom_blk_id);
            // The port name needs to be made unique for the supernode (two
            // joined blocks may have the same port name)
            std::string ap_port_name = block_name + "_" + port_name;
            APPortId ap_port_id = ap_netlist.create_port(ap_blk_id, ap_port_name, port_width, port_type);
            for (AtomPinId atom_pin_id : atom_netlist.port_pins(atom_port_id)) {
                BitIndex port_bit = atom_netlist.pin_port_bit(atom_pin_id);
                PinType pin_type = atom_netlist.pin_type(atom_pin_id);
                bool pin_is_const = atom_netlist.pin_is_constant(atom_pin_id);
                AtomNetId pin_atom_net_id = atom_netlist.pin_net(atom_pin_id);
                const std::string& pin_atom_net_name = atom_netlist.net_name(pin_atom_net_id);
                APNetId pin_ap_net_id = ap_netlist.create_net(pin_atom_net_name);
                ap_netlist.create_pin(ap_port_id, port_bit, pin_ap_net_id, pin_type, pin_is_const);
            }
        }
    }

    // Fix the block locations given by the VPR constraints
    for (APBlockId ap_blk_id : ap_netlist.blocks()) {
        const t_pack_molecule* mol = ap_netlist.block_molecule(ap_blk_id);
        for (AtomBlockId mol_atom_blk_id : mol->atom_block_ids) {
            PartitionId part_id = constraints.get_atom_partition(mol_atom_blk_id);
            if (!part_id.is_valid())
                continue;
            // We should not fix a block twice. This would imply that a molecule
            // contains two fixed blocks. This would only make sense if the blocks
            // were fixed to the same location. I am not sure if that is even
            // possible.
            VTR_ASSERT(ap_netlist.block_type(ap_blk_id) == APBlockType::MOVEABLE);
            const PartitionRegion& partition_pr = constraints.get_partition_pr(part_id);
            VTR_ASSERT(partition_pr.get_regions().size() == 1 && "Each partition should contain only one region!");
            const vtr::Rect<int>& region_rect = partition_pr.get_regions()[0].get_rect();
            VTR_ASSERT(region_rect.xmin() == region_rect.xmax() && "Expect each region to be a single point in x!");
            VTR_ASSERT(region_rect.ymin() == region_rect.ymax() && "Expect each region to be a single point in y!");
            // FIXME: Here we are assuming the user is not using the optional
            //        arguments to fix the sub_tile and layer_num. This should be
            //        expressed! Do not merge unless this is done.
            APFixedBlockLoc loc = {region_rect.xmin(), region_rect.ymin(), -1, -1};
            ap_netlist.set_block_loc(ap_blk_id, loc);
        }
    }

    // Cleanup the netlist by removing undesirable nets.
    // Currently undesirable nets are nets that are:
    //  - ignored for placement
    //  - a global net
    //  - connected to 1 or fewer unique blocks
    //  - connected to only fixed blocks
    // TODO: Allow the user to pass a threshold so we can remove high-fanout nets.
    auto remove_net = [&](APNetId net_id) {
        // Remove all pins associated with this net
        for (APPinId pin_id : ap_netlist.net_pins(net_id))
            ap_netlist.remove_pin(pin_id);
        // Remove the net
        ap_netlist.remove_net(net_id);
    };
    for (APNetId ap_net_id : ap_netlist.nets()) {
        // Is the net ignored for placement, if so remove
        const std::string& net_name = ap_netlist.net_name(ap_net_id);
        AtomNetId atom_net_id = atom_netlist.find_net(net_name);
        VTR_ASSERT(atom_net_id.is_valid());
        if (atom_netlist.net_is_ignored(atom_net_id)) {
            remove_net(ap_net_id);
            continue;
        }
        // Is the net global, if so remove
        if (atom_netlist.net_is_global(atom_net_id)) {
            remove_net(ap_net_id);
            continue;
        }
        // Get the unique blocks connectioned to this net
        std::unordered_set<APBlockId> net_blocks;
        for (APPinId ap_pin_id : ap_netlist.net_pins(ap_net_id)) {
            net_blocks.insert(ap_netlist.pin_block(ap_pin_id));
        }
        // If connected to 1 or fewer unique blocks, remove
        if (net_blocks.size() <= 1) {
            remove_net(ap_net_id);
            continue;
        }
        // If all the connected blocks are fixed, remove
        bool is_all_fixed = true;
        for (APBlockId ap_blk_id : net_blocks) {
            if (ap_netlist.block_type(ap_blk_id) == APBlockType::MOVEABLE) {
                is_all_fixed = false;
                break;
            }
        }
        if (is_all_fixed) {
            remove_net(ap_net_id);
            continue;
        }
    }
    ap_netlist.remove_and_compress();

    // TODO: Some nets may go to the same block mutliple times. Should that be
    //       fixed? Maybe we want nets to be strongly attracted towards some
    //       blocks over others.

    // TODO: Should we cleanup the blocks? For example if there is no path
    //       from a fixed block to a given moveable block, then that moveable
    //       block can be removed (since it can literally go anywhere).
    //  - This would be useful to detect and use throughout; but may cause some
    //    issues if we just remove them. When and where will they eventually
    //    be placed?
    //  - Perhaps we can add a flag in the netlist to each of these blocks and
    //    during the solving stage we can ignore them.

    // FIXME: Should this verify even be here?
    VTR_ASSERT(ap_netlist.verify());
    return ap_netlist;
}

