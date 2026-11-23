#pragma once
#include <dynsyslab/core/Json.hpp>
#include <string>
#include <vector>

namespace dynsyslab::core::lie {
// Exact symbolic proof data, distinct from the floating-point runtime
// dictionary. Coefficients are decimal integers or p/q strings; no silent
// rationalization of doubles.
struct Term {
  std::vector<unsigned> powers;
  std::string coefficient;
};
struct Polynomial {
  std::size_t dimension{};
  std::vector<Term> terms;
};
using VectorField = std::vector<Polynomial>;
struct Limits {
  std::size_t maximumBasis{32}, maximumOperations{4096};
  unsigned maximumDegree{24};
};
struct AlgebraResult {
  std::string status, reason; // closed, limit_reached, proven_infinite
  std::vector<VectorField> basis;
  std::vector<std::string> words;
  std::vector<std::size_t> growth;
  // [i,j] = sum_k structureConstants[(i*r+j)*r+k] basis[k].
  std::vector<std::string> structureConstants;
  std::vector<std::size_t> derivedDimensions, lowerCentralDimensions;
  std::vector<std::vector<std::string>> center;
  std::vector<std::string> killingForm;
  bool solvable{}, nilpotent{};
  std::size_t operations{};
};
struct ModuleResult {
  std::string status, reason;
  std::vector<Polynomial> basis;
  // Operator -> output observable -> predictor observable, L_X Psi = A_X Psi.
  // Matrices acting on coefficient columns are A_X^T (representation sign
  // convention).
  std::vector<std::vector<std::string>> matrices;
  std::vector<std::vector<Polynomial>> residuals;
  std::vector<std::string> linearStateReadout;
  bool linearStateReconstruction{};
};
struct BoxBound {
  std::string uniformL1, lipschitzL1;
  bool vanishesAtOrigin{};
};
// Rigorous coefficient bounds on |x_j| <= radii_j, centered at the origin.
// If residual(0)=0: ||residual(x)||_1 <= lipschitzL1 * ||x||_1.
BoxBound residualBoxBound(const std::vector<Polynomial> &residual,
                          const std::vector<std::string> &radii);
Polynomial differentiate(const Polynomial &, std::size_t axis);
// Exact value of the stored binary floating coefficient, NOT a claim about
// the exact physical parameter or a fitted function's scientific accuracy.
std::string binaryRational(double value);
std::vector<Polynomial> independentPolynomials(const std::vector<Polynomial>&);
Polynomial multiply(const Polynomial &, const Polynomial &);
Polynomial lieDerivative(const VectorField &, const Polynomial &);
VectorField bracket(const VectorField &x, const VectorField &y);
// Constant antisymmetric Poisson tensor; row-major, dimension^2 exact strings.
Polynomial poissonBracket(const Polynomial &, const Polynomial &,
                          const std::vector<std::string> &tensor);
Polynomial flowBracket(const VectorField &, const Polynomial &,
                       const Polynomial &);
AlgebraResult generatedAlgebra(const std::vector<VectorField> &generators,
                               const Limits &limits = {});
// Smallest ideal of a certified finite ambient algebra containing the seeds.
AlgebraResult generatedIdeal(const AlgebraResult &ambient,
                             const std::vector<VectorField> &seeds,
                             const Limits &limits = {});
ModuleResult invariantModule(const std::vector<Polynomial> &seeds,
                             const std::vector<VectorField> &operators,
                             bool extend = true, const Limits &limits = {});
std::size_t rankAt(const std::vector<VectorField> &fields,
                   const std::vector<std::string> &point);
// Exact defects for supplied (possibly approximate rational) model
// coefficients.
std::vector<Polynomial>
representationDefect(const std::vector<Polynomial> &basis, const VectorField &,
                     const std::vector<std::string> &rowMajorMatrix);
// Read-only versioned gate, callable independently of
// solver/workflow/archive/UI.
json::Value analyze(const json::Value &request);
} // namespace dynsyslab::core::lie
