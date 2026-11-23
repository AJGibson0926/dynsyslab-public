#pragma once
#include <cstddef>
#include <optional>
#include <vector>
namespace dynsyslab::math {
struct GraphEdge {std::size_t from{},to{};double length{1};};
struct FiniteGraph {std::size_t vertices{};std::vector<GraphEdge> edges;};
struct GraphTopology {
    std::size_t connectedComponents{},firstBettiNumber{};
    std::vector<std::size_t> component,degree,spanningForestEdges,cycleClosingEdges;
};
// Undirected multigraph, including self-loops and isolated vertices. Graphs
// must be supplied explicitly: this is not a certified Reeb-graph extractor.
GraphTopology graphTopology(const FiniteGraph& graph);
std::vector<std::optional<double>> graphDistances(const FiniteGraph& graph,std::size_t source);
struct PersistenceInterval {
    std::size_t birthVertex{};
    double birth{};
    std::optional<double> death;
};
// H0 lower-star persistence on a graph: edges enter at max(vertex values).
// Equal births use vertex-index tie breaking. Zero-length bars are retained.
std::vector<PersistenceInterval> lowerStarPersistence(const FiniteGraph& graph,
                                                     const std::vector<double>& vertexValues);
}
