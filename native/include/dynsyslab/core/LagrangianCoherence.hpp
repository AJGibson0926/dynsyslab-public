#pragma once
// Lagrangian coherent structures from the finite-time flow map (docs/architecture/LCS_AND_SSM.md).
//
//   Flow map       Phi_{t0}^{t0+T}(x0) of x' = f(t, x), T > 0 forward, T < 0 backward, integrated
//                  from every node of a DomainGrid of initial conditions.
//   Cauchy-Green   C = (D Phi)^T D Phi restricted to the gridded coordinates (a slice when the grid
//                  anchors the other coordinates), eigenvalues 0 < lambda_1 <= ... <= lambda_k.
//   FTLE           sigma = log(lambda_k) / (2 |T|), the finite-time Lyapunov exponent
//                  (Haller 2001, Physica D 149:248-277; Shadden, Lekien and Marsden 2005,
//                  Physica D 212:271-304).
//   Ridges         second-derivative ridges of sigma: nodes where sigma is a strict maximum along the
//                  Hessian direction of most negative curvature (Shadden et al. 2005, section 4).
//   Variational    in two dimensions, repelling LCS are strainlines r' = xi_1(r) of the forward
//   LCS            tensor on which lambda_1 != lambda_2 > 1 and <xi_2, D^2 lambda_2 xi_2> < 0, kept
//                  where their mean lambda_2 is largest among nearby candidates (Haller 2011,
//                  Physica D 240:574-598; Farazmand and Haller 2012, Chaos 22:013128); attracting
//                  LCS are the repelling LCS of the backward flow (Haller 2015, Annu. Rev. Fluid
//                  Mech. 47:137-162; Haller and Yuan 2000, Physica D 147:352-370).
#include <dynsyslab/core/DomainGrid.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace dynsyslab::core {

struct FlowMapGridOptions {
    double initialTime{0};
    double integrationTime{1};  // signed: forward when positive, backward when negative
    double maximumStep{1e-2};   // classical RK4 steps no longer than this
    // 0: the flow-map gradient by central differences between neighbouring grid nodes (Shadden et
    // al. 2005). > 0: an auxiliary pair of initial conditions +-spacing about each node along each
    // gridded axis, more accurate where the flow map varies on scales below the grid spacing
    // (Farazmand and Haller 2012, section 4.1).
    double auxiliarySpacing{0};
};

// The Cauchy-Green tensor at every node of a grid.
struct CauchyGreenField {
    DomainGrid grid;
    double integrationTime{};
    std::size_t axes{};                 // k, the gridded coordinates
    std::vector<double> eigenvalues;    // points x k, ascending
    std::vector<double> eigenvectors;   // points x k x k; eigenvector j is column j (gridded coordinates)
    std::vector<GridStatus> status;     // undefined where an orbit left the valid state set

    [[nodiscard]] std::size_t size() const noexcept { return status.size(); }
    // log(lambda_k) / (2|T|) at each node (NaN where undefined).
    [[nodiscard]] std::vector<double> ftle() const;
};

// Phi at every node (points x state dimension); a NaN row where the orbit left the valid set.
[[nodiscard]] std::vector<double> flowMapOnGrid(const ODESystem& system, const DomainGrid& grid,
                                                const FlowMapGridOptions& options);
[[nodiscard]] CauchyGreenField cauchyGreenField(const ODESystem& system, const DomainGrid& grid,
                                                const FlowMapGridOptions& options);

// Eigen-decomposition of a symmetric k x k matrix by cyclic Jacobi rotations (Golub and Van Loan,
// Matrix Computations, 4th ed. 2013, section 8.5): ascending eigenvalues, eigenvectors in columns.
void symmetricEigen(std::vector<double> matrix, std::size_t k, std::vector<double>& values,
                    std::vector<double>& vectors);

// 1 at second-derivative ridge nodes of `field` (a scalar on a two-axis grid) whose value is at least
// `threshold` times the largest finite value, else 0.
[[nodiscard]] std::vector<std::uint8_t> ridgeNodes(const DomainGrid& grid, const std::vector<double>& field,
                                                   double threshold = 0.3);

struct StrainlineOptions {
    double threshold{0.3};       // seeds: ridge nodes of lambda_2 above this share of its maximum
    std::size_t seedStride{3};   // every n-th admissible seed in each grid direction
    double step{0.25};           // integration step, in grid spacings
    double maximumLength{0};     // 0: the domain diagonal
    double minimumLength{0};     // 0: a fifth of the shorter domain side
    double failureLength{0};     // stop once the conditions fail this long; 0: two grid spacings
    double separation{0};        // a curve within this of a stronger one is dropped; 0: two spacings
};

struct LagrangianCurve {
    std::vector<double> points;  // samples x 2, in the grid's two coordinates
    double length{};
    double meanStretching{};     // mean lambda_2 along the curve
};

// Variational LCS of a two-axis Cauchy-Green field: strainlines (repelling for T > 0, attracting
// for T < 0), strongest first.
[[nodiscard]] std::vector<LagrangianCurve> strainlines(const CauchyGreenField& field,
                                                       const StrainlineOptions& options = {});

// 1 at grid nodes within half a grid cell of a curve sample, else 0 (for display on the grid).
[[nodiscard]] std::vector<std::uint8_t> rasterizeCurves(const DomainGrid& grid,
                                                        const std::vector<LagrangianCurve>& curves);

}  // namespace dynsyslab::core
