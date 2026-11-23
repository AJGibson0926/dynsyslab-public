#pragma once
#include <dynsyslab/math/Topology.hpp>
#include <cstddef>
#include <optional>
#include <vector>
namespace dynsyslab::math {
// Scalar topology on a regular 2D grid, triangulated by the Freudenthal (Kuhn)
// split. Ties break by vertex index (simulation of simplicity), so plateaus get
// a consistent but arbitrary resolution and may yield zero-persistence pairs.
// Persistence: Edelsbrunner, Letscher & Zomorodian, Discrete Comput. Geom. 28
// (2002), doi:10.1007/s00454-002-2885-2. Join/split sweep: Carr, Snoeyink &
// Axen, Comput. Geom. 24 (2003), doi:10.1016/S0925-7721(02)00093-7. Because a
// periodic domain can have loops, arcs come from slab components rather than
// from the contour-tree merge, which is valid only for simply connected domains.
struct ScalarGrid2D {
    std::size_t nx{},ny{};
    std::vector<double> values; // Row-major: index=i+nx*j.
    bool periodicX{false},periodicY{false};
};
enum class GridCriticalKind { Regular, Minimum, Saddle, Maximum };
const char* gridCriticalKindName(GridCriticalKind kind);
struct GridCriticalPoint {
    std::size_t vertex{},i{},j{};
    double value{};
    GridCriticalKind kind{GridCriticalKind::Regular};
    unsigned lowerLinkComponents{},upperLinkComponents{},multiplicity{};
    bool onBoundary{false};
    int eulerContribution{}; // 1-chi(lower link); sums to chi(domain).
};
struct GridPersistencePair {
    unsigned dimension{};
    std::size_t birthVertex{};
    double birth{};
    std::optional<std::size_t> deathVertex; // Empty for an essential class.
    std::optional<double> death;
};
struct LevelGraphNode {
    std::size_t vertex{};
    double value{};
    GridCriticalKind kind{GridCriticalKind::Regular};
    bool onBoundary{false};
    std::size_t lowerDegree{},upperDegree{};
};
struct LevelGraphArc { std::size_t lowerNode{},upperNode{}; };
struct ScalarLevelGraph {
    std::vector<LevelGraphNode> nodes; // Ascending filtration order.
    std::vector<LevelGraphArc> arcs;
};
// Edge lengths are value differences, ready for graphTopology and friends.
FiniteGraph toFiniteGraph(const ScalarLevelGraph& graph);
struct SeparatrixLevel {
    std::size_t vertex{};
    double value{};
    bool onBoundary{false};
    std::size_t lowerArcs{},upperArcs{};
};
struct ScalarTopology {
    std::vector<GridCriticalPoint> criticalPoints;
    int eulerCharacteristic{};
    std::vector<GridPersistencePair> persistence;
    ScalarLevelGraph mergeTree,reebGraph;
    std::size_t reebFirstBetti{};
    std::vector<SeparatrixLevel> separatrices;
};
// Pairs with persistence below minimumPersistence are dropped from the report;
// the graphs are not simplified. At most 1,000,000 vertices.
ScalarTopology analyzeScalarGrid(const ScalarGrid2D& grid,double minimumPersistence=0);
}
