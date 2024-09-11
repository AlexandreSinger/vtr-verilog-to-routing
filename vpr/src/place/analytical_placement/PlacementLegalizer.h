
#pragma once

#include "PartialPlacement.h"

class APNetlist;

class PlacementLegalizer {
public:
    PlacementLegalizer(const APNetlist& inetlist) : netlist(inetlist) {}
    virtual ~PlacementLegalizer() {}
    virtual void legalize(PartialPlacement &p_placement) = 0;
protected:
    const APNetlist& netlist;
};

class FlowBasedLegalizer : public PlacementLegalizer {
    using PlacementLegalizer::PlacementLegalizer;
public:
    void legalize(PartialPlacement &p_placement) final;
};

