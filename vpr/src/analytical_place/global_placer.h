
#pragma once

#include <memory>

class APNetlist;
class AnalyticalSolver;
class PartialPlacement;
class PartialLegalizer;

enum class e_global_placer {
    SimPL
};

class GlobalPlacer {
public:
    virtual ~GlobalPlacer() {}

    GlobalPlacer(const APNetlist& netlist, int log_verbosity = 1)
                    : netlist_(netlist),
                      log_verbosity_(log_verbosity) {}

    virtual PartialPlacement place() = 0;

protected:

    const APNetlist& netlist_;

    int log_verbosity_;
};

std::unique_ptr<GlobalPlacer> make_global_placer(e_global_placer placer_type,
                                                 const APNetlist& netlist);

class SimPLGlobalPlacer : public GlobalPlacer {
    std::unique_ptr<AnalyticalSolver> solver_;
    std::unique_ptr<PartialLegalizer> partial_legalizer_;
public:
    SimPLGlobalPlacer(const APNetlist& netlist);

    PartialPlacement place() final;
};

