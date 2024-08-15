
#pragma once

class AtomContext;
class APPrepacker;
class UserPlaceConstraints;
class APNetlist;

APNetlist read_atom_netlist(const AtomContext& atom_ctx, const APPrepacker& prepacker, const UserPlaceConstraints& constraints);

