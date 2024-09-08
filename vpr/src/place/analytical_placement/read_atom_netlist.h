
#pragma once

class AtomContext;
class Prepacker;
class UserPlaceConstraints;
class APNetlist;

APNetlist read_atom_netlist(const AtomContext& atom_ctx, const Prepacker& prepacker, const UserPlaceConstraints& constraints);

