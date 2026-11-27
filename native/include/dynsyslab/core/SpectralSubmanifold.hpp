#pragma once
// Spectral submanifolds (SSMs) and their reduced dynamics (docs/architecture/LCS_AND_SSM.md).
//
// An SSM W(E) of a hyperbolic equilibrium x0 of x' = f(x) is the smoothest invariant manifold
// tangent to a spectral subspace E of A = Df(x0) (Haller and Ponsioen 2016, Nonlinear Dyn.
// 86:1493-1534). The parameterization method writes it as x = x0 + W(p), p in C^d, with reduced
// dynamics p' = R(p), and solves the invariance equation f(x0 + W(p)) = DW(p) R(p) order by order
// (Cabre, Fontich and de la Llave 2003, Indiana Univ. Math. J. 52:283-328; Haro et al., The
// Parameterization Method for Invariant Manifolds, Springer 2016). At order k each monomial p^m
// gives the cohomological equation
//     (A - (m . lambda) I) W_m = V_E R_m + [DW R]_m - [f(x0 + W)]_m
// solved in the eigenbasis of A. In normal-form style R_m keeps exactly the inner (near-)resonant
// terms, lambda_j ~ m . lambda (Jain and Haller 2022, Nonlinear Dyn. 107:1417-1450; for lightly
// damped modes the comparison is of imaginary parts, Ponsioen, Pedergnana and Haller 2018,
// J. Sound Vib. 420:269-295); in graph style R_m takes the whole master-mode projection.
//
// For a two-dimensional SSM of a complex pair the normal form is p' = lambda p + sum_j g_j p^{j+1}
// conj(p)^j, so with p = rho e^{i theta}: rho' = Re(lambda) rho + sum Re(g_j) rho^{2j+1} and
// theta' = Im(lambda) + sum Im(g_j) rho^{2j}, the backbone curve (frequency and damping against
// amplitude; Szalai, Ehrhardt and Haller 2017, Proc. R. Soc. A 473:20160759).
//
// Data-driven SSMs (Cenedese, Axas, Bauerlein, Avila and Haller 2022, Nat. Commun. 13:872; SSMLearn):
// the tangent space from the leading singular vectors of trajectory data, a polynomial graph
// y = y0 + V eta + w(eta) by least squares, polynomial reduced dynamics eta' = r(eta) by regression
// on numerical derivatives, and its normal form by the parameterization method above.
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/NumericalDifferentiation.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <complex>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class SsmStyle { normalForm, graph };

struct SpectralSubmanifoldOptions {
    std::size_t order{5};
    // Indices into the eigenvalues sorted by decreasing real part (slowest first). Empty: mode 0 and,
    // when complex, its conjugate. A complex master mode always brings its conjugate.
    std::vector<std::size_t> masterModes;
    SsmStyle style{SsmStyle::normalForm};
    // Inner resonance lambda_j ~ m . lambda when the difference is within this share of max |lambda_E|.
    double resonanceTolerance{1e-6};
    // Compare imaginary parts only (near-resonance of lightly damped modes).
    bool imaginaryResonance{true};
};

struct SpectralSubmanifold {
    std::size_t stateDimension{}, dimension{}, order{};
    std::vector<double> equilibrium;
    std::vector<std::complex<double>> eigenvalues;  // all, by decreasing real part
    std::vector<std::size_t> masterModes;           // into eigenvalues, in the order of p
    std::vector<std::size_t> exponents;             // monomials x dimension (degree 0..order, graded)
    std::vector<std::complex<double>> W;            // monomials x stateDimension: x - x0
    std::vector<std::complex<double>> R;            // monomials x dimension: p'
    // Real chart q of the master coordinates: p = chart q (p_j = q_a + i q_b, conj p_j = q_a - i q_b
    // for a complex pair; p_j = q_j for a real mode).
    std::vector<std::complex<double>> chart;        // dimension x dimension
    double spectralQuotient{};                      // Int(min Re lambda outside E / max Re lambda in E)
    std::vector<std::string> warnings;              // near outer resonances, undamped modes

    [[nodiscard]] std::size_t monomials() const noexcept { return dimension ? exponents.size() / dimension : 0; }
    [[nodiscard]] std::vector<std::complex<double>> parameterization(const std::complex<double>* p) const;  // x
    [[nodiscard]] std::vector<std::complex<double>> reducedDynamics(const std::complex<double>* p) const;   // p'
    // Real polynomial coefficients in q over polynomialExponentTable(dimension, order, true):
    // features x stateDimension for x (with x0 as the constant), features x dimension for q'.
    [[nodiscard]] std::vector<double> realParameterization() const;
    [[nodiscard]] std::vector<double> realReducedDynamics() const;
    // Largest |f(W(p)) - DW(p) R(p)| / |f| scale over p on the circle |p_j| = radius (invariance check).
    [[nodiscard]] double invarianceError(const std::vector<AnalyticExpression>& field, double radius) const;
};

// The SSM of `field` (expressions in the state) about `equilibrium`.
[[nodiscard]] SpectralSubmanifold spectralSubmanifold(const std::vector<AnalyticExpression>& field,
                                                      const std::vector<double>& equilibrium,
                                                      const SpectralSubmanifoldOptions& options = {});

// Newton's method for f(x) = 0 from `guess` with the symbolic Jacobian.
[[nodiscard]] std::vector<double> equilibriumNear(const std::vector<AnalyticExpression>& field, std::vector<double> guess);

// A polynomial vector field q' = sum_m c_m q^m as expressions (coefficients features x dimension
// over polynomialExponentTable(dimension, degree, true)).
[[nodiscard]] std::vector<AnalyticExpression> polynomialField(std::size_t dimension, std::size_t degree,
                                                              const std::vector<double>& coefficients);

struct BackboneCurve {
    std::vector<double> radius;     // rho
    std::vector<double> frequency;  // theta'(rho), rad per time
    std::vector<double> damping;    // -rho'/(rho |theta'|)
    std::vector<double> amplitude;  // radius x outputs: max over theta of |output(x0 + W) - output(x0)|
    std::size_t outputs{};
};

using StateOutputs = std::function<std::vector<double>(const double* state)>;

// The backbone of a two-dimensional normal-form SSM of a complex pair, at the given radii.
[[nodiscard]] BackboneCurve backboneCurve(const SpectralSubmanifold& ssm, const std::vector<double>& radii,
                                          const StateOutputs& outputs, std::size_t outputCount,
                                          std::size_t angles = 128);

struct SsmFitOptions {
    std::size_t dimension{2};
    std::size_t order{3};            // manifold polynomial degree
    std::size_t dynamicsOrder{3};    // reduced-dynamics polynomial degree
    std::size_t normalFormOrder{0};  // 0: dynamicsOrder
    double ridge{0};
    double skipTime{0};              // discard each trajectory's first samples (off-manifold transient)
    std::vector<double> equilibrium; // empty: the mean of the trajectories' final samples
    DerivativeSpec derivative;       // time derivatives of the reduced coordinates
};

struct SsmFit {
    std::size_t observableDimension{}, dimension{}, order{}, dynamicsOrder{};
    std::vector<double> equilibrium;            // y0
    std::vector<double> tangent;                // observables x dimension, orthonormal columns V
    std::vector<double> manifold;               // features x observables over polynomialExponentTable(d, order, true)
    std::vector<double> dynamics;               // features x dimension over polynomialExponentTable(d, dynamicsOrder, true)
    double manifoldError{};                     // relative RMS of the graph fit
    double dynamicsError{};                     // relative RMS of the reduced-dynamics regression
    std::vector<double> singularValues;         // of the data, for choosing the dimension
    SpectralSubmanifold normalForm;             // of the reduced dynamics (state: eta)

    [[nodiscard]] std::vector<double> lift(const double* eta) const;     // y
    [[nodiscard]] std::vector<double> reduce(const double* y) const;     // eta = V^T (y - y0)
    [[nodiscard]] std::vector<double> rate(const double* eta) const;     // eta'
    // The reduced model from the first sample of `trajectory` lifted to the observables, on its clock.
    [[nodiscard]] TimeSeries predict(const TimeSeries& trajectory, double maximumStep = 0) const;
};

[[nodiscard]] SsmFit fitSpectralSubmanifold(const std::vector<TimeSeries>& trajectories, const SsmFitOptions& options = {});

}  // namespace dynsyslab::core
