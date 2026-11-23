#include <cmath>
#include <dynsyslab/core/KromDifferentialAnalysis.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
// Generator/input defects: Strasser et al., arXiv:2312.01441 and SafEDMD
// arXiv:2402.03145. These are evaluated deterministic defects, without claiming
// their statistical, regional, or stability theorems from finite point samples.
KromDifferentialAnalysis
analyzeKromDifferentials(const KoopmanModel &model,
                         const std::vector<double> &states,
                         const std::string &inputMap) {
  if (model.requiresExternalCoordinates() || !model.hasPhysicalSystem())
    throw std::invalid_argument("Physical differential analysis requires an "
                                "instantaneous physical-state provider");
  const auto d = model.stateDimension(), n = model.splitDimension();
  if (states.empty() || states.size() % d || model.inputDimension() != 1)
    throw std::invalid_argument(
        "Require state rows and a scalar physical input");
  for (double x : states)
    if (!std::isfinite(x))
      throw std::invalid_argument("Differential states must be finite");
  auto system = createSystem(model.systemConfiguration());
  const auto generator = model.generatorSplit();
  KromDifferentialAnalysis out;
  out.samples = states.size() / d;
  out.states = d;
  out.lifted = n;
  for (std::size_t s = 0; s < out.samples; ++s) {
    const auto *x = states.data() + s * d;
    if (system->invalid(0, x))
      throw std::invalid_argument(
          "Differential sample outside system validity domain");
    const auto z = model.evaluateSplit(x), j = model.stateJacobianSplit(x);
    std::vector<double> f(d);
    system->rhsWithInput(0, x, 0, f.data());
    const auto g = system->inputMatrix(x),
               b = model.controlMapSplit(x, inputMap);
    for (std::size_t i = 0; i < n; ++i) {
      double v = 0, u = 0, az = 0;
      for (std::size_t a = 0; a < d; ++a) {
        v += j[i * d + a] * f[a];
        u += j[i * d + a] * g[a];
      }
      for (std::size_t k = 0; k < n; ++k)
        az += generator[i * n + k] * z[k];
      out.drift.push_back(v);
      out.generatorDefect.push_back(v - az);
      out.trueInputMap.push_back(u);
      out.fittedInputMap.push_back(b[i]);
      out.inputDefect.push_back(u - b[i]);
    }
    out.values.insert(out.values.end(), z.begin(), z.end());
    out.jacobians.insert(out.jacobians.end(), j.begin(), j.end());
    const auto sv = numerics::singularValueDecomposition(j, n, d);
    out.smallestSingularValue.push_back(n < d ? 0 : sv.singularValues.back());
    out.largestSingularValue.push_back(sv.singularValues.front());
    if (n == 2 && d == 2)
      out.orientedDeterminant.push_back(j[0] * j[3] - j[1] * j[2]);
  }
  return out;
}
InvariantRegression regressInvariantAction(
    const std::vector<double> &train, const std::vector<double> &target,
    const std::vector<double> &valid, const std::vector<double> &validTarget,
    std::size_t n, std::size_t m, double tol) {
  if (!n || !m || train.empty() || valid.empty() || train.size() % n ||
      valid.size() % n || target.size() != train.size() / n * m ||
      validTarget.size() != valid.size() / n * m || !std::isfinite(tol) ||
      tol <= 0)
    throw std::invalid_argument(
        "Invariant regression data shape/tolerance mismatch");
  for (const auto *values : {&train, &target, &valid, &validTarget})
    for (double x : *values)
      if (!std::isfinite(x))
        throw std::invalid_argument("Regression samples must be finite");
  auto augment = [&](const std::vector<double> &x) {
    std::vector<double> a;
    for (std::size_t s = 0; s < x.size() / n; ++s) {
      a.push_back(1);
      a.insert(a.end(), x.begin() + s * n, x.begin() + (s + 1) * n);
    }
    return a;
  };
  const auto a = augment(train), v = augment(valid);
  numerics::LeastSquaresOptions opt;
  opt.method = numerics::LeastSquaresMethod::divideAndConquerSvd;
  opt.singularValueTolerance = tol;
  InvariantRegression out;
  out.predictors = n + 1;
  out.outputs = m;
  out.coefficients =
      numerics::solveLeastSquares(a, train.size() / n, n + 1, target, m, opt);
  auto error = [&](const std::vector<double> &x, const std::vector<double> &y) {
    double e = 0, den = 0;
    for (std::size_t s = 0; s < y.size() / m; ++s)
      for (std::size_t i = 0; i < m; ++i) {
        double pred = 0;
        for (std::size_t k = 0; k < n + 1; ++k)
          pred += x[s * (n + 1) + k] * out.coefficients[k * m + i];
        e += (pred - y[s * m + i]) * (pred - y[s * m + i]);
        den += y[s * m + i] * y[s * m + i];
      }
    return den > 0 ? std::sqrt(e / den)
                   : (e == 0 ? 0 : std::numeric_limits<double>::infinity());
  };
  out.trainingRelativeError = error(a, target);
  out.validationRelativeError = error(v, validTarget);
  return out;
}
} // namespace dynsyslab::core
