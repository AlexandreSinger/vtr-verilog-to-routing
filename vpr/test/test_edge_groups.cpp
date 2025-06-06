#include <vector>
#include <utility>
#include <set>
#include <random>
#include <algorithm>

#include "catch2/catch_test_macros.hpp"

#include "edge_groups.h"

namespace {

TEST_CASE("edge_groups_create_sets", "[vpr]") {
    // Construct a set of edges that result in these connected sets
    std::vector<std::set<int>> connected_sets{{{1, 2, 3, 4, 5, 6, 7, 8},
                                               {9, 0}}};

    // Build chains from the given connected sets
    int max_node_id = 0;
    std::vector<std::pair<int, int>> edges;
    for (const auto& set : connected_sets) {
        int last = *set.cbegin();
        std::for_each(std::next(set.cbegin()),
                      set.cend(),
                      [&](int node) {
                          edges.emplace_back(last, node);
                          last = node;
                          max_node_id = std::max(max_node_id, node);
                      });
    }

    // Build the id map for node IDs
    std::vector<int> nodes(max_node_id + 1);

    // Initialize nodes to [0, 1, ..., max_node_id]
    std::iota(nodes.begin(), nodes.end(), 0);

    // Create a Mersenne Twister pseudo-random number generator with seed 1
    std::mt19937 g(1);

    // Run the test many times, the PRNG will give differently shuffled inputs
    for (int i = 0; i < 1000; i++) {
        // Shuffle node IDs
        auto random_nodes = nodes;
        std::shuffle(random_nodes.begin(), random_nodes.end(), g);

        // Apply shuffled IDs to edges
        auto random_edges = edges;
        for (auto& edge : random_edges) {
            edge.first = random_nodes[edge.first];
            edge.second = random_nodes[edge.second];
        }

        // Shuffle edges
        std::shuffle(random_edges.begin(), random_edges.end(), g);

        // Add edges to the EdgeGroups object
        EdgeGroups groups;
        for (auto edge : random_edges) {
            groups.add_non_config_edge(RRNodeId(edge.first), RRNodeId(edge.second));
        }

        // The algorithm to test
        groups.create_sets();
        t_non_configurable_rr_sets sets = groups.output_sets();

        // Check for the expected sets
        for (const auto& set : connected_sets) {
            std::set<RRNodeId> random_set;
            for (auto elem : set) {
                random_set.insert(RRNodeId(random_nodes[elem]));
            }
            REQUIRE(std::find(sets.node_sets.begin(), sets.node_sets.end(), random_set) != sets.node_sets.end());
        }
    }
}

} // namespace
