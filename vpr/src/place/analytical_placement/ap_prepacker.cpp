
#include "ap_prepacker.h"
#include <map>
#include "prepack.h"
#include "vpr_context.h"
#include "vpr_types.h"
#include "vtr_assert.h"
#include "vtr_log.h"

APPrepacker::APPrepacker(AtomContext& mutable_atom_ctx) {
    VTR_LOG("Pre-Packing Atom Netlist\n");
    // Run the prepacker.
    list_of_pack_patterns = alloc_and_load_pack_patterns();
    mutable_atom_ctx.list_of_pack_molecules.reset(alloc_and_load_pack_molecules(list_of_pack_patterns.data(), expected_lowest_cost_pb_gnode, list_of_pack_patterns.size()));
}

APPrepacker::~APPrepacker() {
    // When the APPrepacker object is destroyed, clean up the list of pack patterns.
    free_list_of_pack_patterns(list_of_pack_patterns);
}

t_pack_molecule* APPrepacker::get_atom_molecule(AtomBlockId atom_blk_id, const AtomContext& atom_ctx) const {
    // The prepacking performed in the constructor populates the atom_molecules
    // multimap within the atom context.
    const std::multimap<AtomBlockId, t_pack_molecule*>& atom_molecules = atom_ctx.atom_molecules;
    // Get the molecule associated with this atom.
    // This assumes that the map only has one entry for each block. Have not hit
    // a situation where this was not true.
    VTR_ASSERT(atom_molecules.count(atom_blk_id) == 1 && "Only expect one molecule for a given block.");
    auto mol_range = atom_molecules.equal_range(atom_blk_id);
    return mol_range.first->second;
}

