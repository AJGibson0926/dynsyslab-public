#pragma once
#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::core {
struct KromDifferentialAnalysis {
  std::size_t samples{}, states{}, lifted{};
  // Sample-major arrays, split-real coordinates. Jacobians: sample x lifted x
  // state.
  std::vector<double> values, jacobians, drift, generatorDefect;
  std::vector<double> trueInputMap, fittedInputMap, inputDefect;
  std::vector<double> smallestSingularValue, largestSingularValue;
  // Only populated for a 2-by-2 state/lift Jacobian.
  std::vector<double> orientedDeterminant;
};
// Pointwise defects, not regional certificates. Physical state only; the ODE's
// autonomous drift and input field are evaluated on the supplied points.
// Nonzero driver configuration is ignored by using rhsWithInput(...,0).
KromDifferentialAnalysis
analyzeKromDifferentials(const KoopmanModel &model,
                         const std::vector<double> &states,
                         const std::string &inputMap = "exact");
struct InvariantRegression {
  std::size_t predictors{}, outputs{};
  std::vector<double> coefficients; // predictor x output
  double trainingRelativeError{}, validationRelativeError{};
};
// Fixed affine feature span [1,z]. Explicit disjoint training/validation data
// avoid reporting interpolation as a closure proof. No data normalization is
// refit.
InvariantRegression
regressInvariantAction(const std::vector<double> &trainingFeatures,
                       const std::vector<double> &trainingAction,
                       const std::vector<double> &validationFeatures,
                       const std::vector<double> &validationAction,
                       std::size_t features, std::size_t outputs,
                       double relativeTolerance = 1e-10);
} // namespace dynsyslab::core
