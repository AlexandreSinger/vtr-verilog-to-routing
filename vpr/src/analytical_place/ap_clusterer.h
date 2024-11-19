

#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "ap_netlist_fwd.h"
#include "cluster_legalizer.h"
#include "vpr_types.h"
#include "vtr_strong_id.h"
#include "vtr_vector.h"

class APNetlist;
class PartialPlacement;
struct t_arch;
struct t_packer_opts;
struct t_vpr_setup;

class APClusterer {
public:
    virtual ~APClusterer() {}

    APClusterer(const APNetlist& netlist, int log_verbosity = 0)
                    : netlist_(netlist),
                      log_verbosity_(log_verbosity) {}

    virtual void create_clusters(const PartialPlacement& p_placement) = 0;

protected:

    const APNetlist& netlist_;

    int log_verbosity_;
};

struct ap_gain_cluster_id_tag;
typedef vtr::StrongId<ap_gain_cluster_id_tag, size_t> APGainClusterId;

struct APGainCluster {
    // Valid will signify if the given cluster is no longer in use (i.e. it was
    // merged into another cluster). This is useful for debugging since every
    // molecule will be in its own cluster at the beginning.
    // FIXME: It may be cleaner and more obvious to just have a map in the
    //        calculator class of booleans. "is contained in cluster".
    bool valid = true;
    // FIXME: The external pins are likely going to be unused. We only need
    //        the internal pins.
    std::unordered_set<AtomPinId> external_pins;
    // Internal pins are pins conencted to blocks within the cluster.
    std::unordered_set<AtomPinId> internal_pins;

    // Nets that have at least one pin outside of the cluster.
    std::unordered_set<AtomNetId> external_nets;
    // Nets that have all pins inside the cluster.
    std::unordered_set<AtomNetId> internal_nets;

    std::vector<const t_pack_molecule*> molecules;

    APGainCluster() {}

    APGainCluster(const t_pack_molecule* mol,
                  const AtomNetlist& netlist);
};

class APClusterGainCalculator {
private:
    static constexpr float pin_gain_weight_ = 1.f;

    // FIXME: SET THIS BACK FROM 0!!!
    static constexpr float wl_gain_weight_ = 0.f;

    const AtomNetlist& atom_netlist_;
    const Prepacker& prepacker_;
    const PartialPlacement& p_placement_;

    vtr::vector<APGainClusterId, APGainCluster> gain_clusters_;

    std::unordered_map<const t_pack_molecule*, APGainClusterId> molecule_gain_cluster_;

public:

    APClusterGainCalculator(const AtomNetlist& atom_netlist,
                            const Prepacker& prepacker,
                            const PartialPlacement& p_placement);

    APGainClusterId create_gain_cluster(const t_pack_molecule* mol);

    // Get the gain of adding molecule mol to the gain cluster.
    float get_gain(APGainClusterId gain_cluster_id, const t_pack_molecule* mol);

    // Add molecule mol to the gain cluster.
    void add_mol_to_gain_cluster(const t_pack_molecule* mol, APGainClusterId gain_cluster_id);

    // Destroy the cluster. Implies that the molecules inside will be clustered
    // into other clusters later.
    void destroy_gain_cluster(APGainClusterId gain_cluster_id);

    // Cleans the given gain cluster assuming that it will never be used again
    // and the molecules inside will never be used again.
    void clean_gain_cluster(APGainClusterId gain_cluster_id);
};

class GreedyAPClusterer : public APClusterer {
private:

    t_pack_high_fanout_thresholds high_fanout_thresholds_;

    // VPR setup used to load the packing.
    // FIXME: This should be removed. Should only bring in what we need.
    t_vpr_setup& vpr_setup_;

    // Arch used to output the clustering.
    const t_arch* arch_;

    // Atom netlist used for getting the model of atoms for starting new clusters.
    const AtomNetlist& atom_netlist_;

    // Prepacker used for the gain calculator.
    const Prepacker& prepacker_;

    // Packer opts used to create the legalizer and output the clustering.
    const t_packer_opts& packer_opts_;

    ClusterLegalizer cluster_legalizer_;

    std::unique_ptr<APClusterGainCalculator> gain_calculator_;

    std::map<const t_model*, std::vector<t_logical_block_type_ptr>> primitive_candidate_block_types_;

    // TODO: Make a APClustererGainCalculator class to abstract the gain
    //       computation in a more clean way.

    // Map from molecule to APBlock. This is necessary to identify which APBlocks
    // are in each cluster.
    // FIXME: This should be maintained by the APNetlist class.
    // FIXME: We can save a lot of time if each molecule had an ID.
    std::unordered_map<const t_pack_molecule*, APBlockId> mol_block_;

    APBlockId select_seed_block();

    LegalizationClusterId start_new_cluster(APBlockId seed);

    bool add_block_to_cluster(APBlockId block, LegalizationClusterId cluster_id);

    bool grow_cluster(APBlockId seed, ClusterLegalizationStrategy strategy);

    // Compatible: unclustered + has a spot in the cluster to be placed.
    APBlockId get_highest_gain_compatible_neighbor(APGainClusterId gain_cluster_id,
                                                   LegalizationClusterId leg_cluster_id);

public:

    GreedyAPClusterer(const APNetlist& netlist,
                      t_vpr_setup& vpr_setup,
                      const t_arch* arch,
                      const AtomNetlist& atom_netlist,
                      const Prepacker& prepacker,
                      const std::vector<t_logical_block_type>& logical_block_types,
                      std::vector<t_lb_type_rr_node>* lb_type_rr_graphs,
                      const t_model* user_models,
                      const t_model* library_models,
                      const t_packer_opts& packer_opts);

    void create_clusters(const PartialPlacement& p_placement) final;
};

