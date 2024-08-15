
#pragma once

#include <unordered_map>
#include <vector>
#include "atom_netlist_fwd.h"

class AtomContext;
class t_pack_molecule;
class t_pack_patterns;
class t_pb_graph_node;

// The APPrepacker encapsulates the prepacker functionality and controls its
// state (so we can avoid messing with global state).
// This will prepack the atoms into molecules, where it will create "forced"
// packs (LUT+FF) and carry chains.
// See /pack/prepack.h and /pack/pack.cpp:try_pack
// Two APPrepacker objects should NOT exist at the same time in a program (since
// they maintain global state).
// FIXME: This is an annoying limitation. Try to find a way to better enforce
//        it.
struct APPrepacker {
    // This class is fiddling with the global state, and as such should not be
    // copied or moved.
    APPrepacker() = delete;
    APPrepacker(const APPrepacker&) = delete;
    APPrepacker& operator=(const APPrepacker&) = delete;
    // Constructor of the prepacker. Performs pre-packing.
    // Populates internal structures within the atom context.
    APPrepacker(AtomContext& mutable_atom_ctx);
    // Destructor. Frees the list of pack patterns. Does not delete the internal
    // state of the atom context.
    ~APPrepacker();
    // Method to get the molecule an atom has be prepacked into.
    t_pack_molecule* get_atom_molecule(AtomBlockId atom_blk_id, const AtomContext& atom_ctx) const;

    std::vector<t_pack_patterns> list_of_pack_patterns;
    // This is a map from the AtomBlock to its expected lowest cost primitive.
    std::unordered_map<AtomBlockId, t_pb_graph_node*> expected_lowest_cost_pb_gnode;
};

