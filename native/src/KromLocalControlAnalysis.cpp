#include <cmath>
#include <dynsyslab/core/KoopmanStateSpace.hpp>
#include <dynsyslab/core/KromLocalControlAnalysis.hpp>
#include <stdexcept>

namespace dynsyslab::core {
// Apply the SAME compensator to the learned and physical local plants.
// Doyle (1978), Guaranteed Margins for LQG Regulators; Strasser et al.,
// arXiv:2509.02839: model accuracy and closed-loop guarantees are separate.
KromLocalControlAnalysis
analyzeKromLocalControl(const KoopmanModel &model,
                        const KromLocalControlOptions &o) {
  if (model.requiresExternalCoordinates() || !model.hasPhysicalSystem() ||
      o.point.size() != model.stateDimension() || model.inputDimension() != 1 ||
      !std::isfinite(o.holdingInput) || o.frequencyPoints < 100 ||
      o.frequencyPoints > 100000)
    throw std::invalid_argument(
        "Local control analysis requires a scalar-input instantaneous model "
        "and valid point/grid");
  for (double v : {o.statePenalty, o.inputPenalty, o.processCovariance,
                   o.measurementCovariance})
    if (!std::isfinite(v) || v <= 0)
      throw std::invalid_argument(
          "Local design weights must be positive and finite");
  for (double x : o.point)
    if (!std::isfinite(x))
      throw std::invalid_argument("Operating point must be finite");
  if (o.holdingInput != 0 && o.controlMap != "constant")
    throw std::invalid_argument(
        "Nonzero holding input requires constant-map local audit; "
        "state-dependent input Jacobian not frozen consistently");
  auto system = createSystem(model.systemConfiguration());
  if (system->invalid(0, o.point.data()))
    throw std::invalid_argument("Invalid local operating point");
  const auto d = model.stateDimension(), n = model.splitDimension();
  auto drift = system->driftExpressions(), input = system->inputExpressions();
  if (drift.size() != d || input.size() != d)
    throw std::invalid_argument(
        "Physical local analysis requires analytic drift/input Jacobians");
  std::vector<double> a(d * d), c;
  for (std::size_t i = 0; i < d; ++i)
    for (std::size_t j = 0; j < d; ++j)
      a[i * d + j] =
          (drift[i].derivative(j) + o.holdingInput * input[i].derivative(j))
              .evaluate(o.point.data());
  const auto catalog = system->observableCatalog();
  const auto &measurement = catalog.at(o.measurement);
  if (measurement.outputDimension() != 1)
    throw std::invalid_argument(
        "Local margin comparison requires a scalar measurement");
  for (auto v : measurement.jacobian(o.point.data())) {
    if (std::abs(v.imag()) > 1e-12)
      throw std::invalid_argument("Physical sensor must be real");
    c.push_back(v.real());
  }
  KromLocalControlAnalysis out;
  out.physical =
      makeStateSpace(a, system->inputMatrix(o.point.data()), c, {}, d, 1, 1);
  out.physicalEquilibriumResidual.resize(d);
  system->rhsWithInput(0, o.point.data(), o.holdingInput,
                       out.physicalEquilibriumResidual.data());
  KoopmanStateSpaceOptions config;
  config.measurements = {o.measurement};
  config.controlMap = o.controlMap;
  config.operatingState = o.point;
  out.surrogate = koopmanStateSpace(model, config);
  out.liftedAtPoint = model.evaluateSplit(o.point.data());
  out.generatorAtPoint.assign(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j)
      out.generatorAtPoint[i] +=
          out.surrogate.A[i * n + j] * out.liftedAtPoint[j];
    out.generatorAtPoint[i] += out.surrogate.B[i] * o.holdingInput;
  }
  out.measuredAtPoint = measurement.evaluateReal(o.point.data());
  for (auto z : model.reconstruct(o.measurement, o.point.data()))
    out.decodedMeasurementAtPoint.push_back(z.real());
  if (model.hasReconstruction("state"))
    for (auto z : model.reconstruct("state", o.point.data()))
      out.decodedStateAtPoint.push_back(z.real());
  auto diagonal = [](std::size_t count, double value) {
    std::vector<double> v(count * count);
    for (std::size_t i = 0; i < count; ++i)
      v[i * count + i] = value;
    return v;
  };
  const auto design =
      designLQG(out.surrogate, diagonal(n, o.statePenalty), {o.inputPenalty},
                diagonal(n, o.processCovariance), {o.measurementCovariance});
  out.gain = design.regulator.gain;
  out.estimatorGain = design.estimator.gain;
  out.compensator = design.compensator;
  const auto j = model.stateJacobianSplit(o.point.data());
  out.physicalGain.assign(d, 0);
  for (std::size_t k = 0; k < d; ++k)
    for (std::size_t i = 0; i < n; ++i)
      out.physicalGain[k] += out.gain[i] * j[i * d + k];
  out.surrogateLqr = loopRobustness(
      stateFeedbackInputLoop(out.surrogate, out.gain), o.frequencyPoints);
  out.physicalLqr =
      loopRobustness(stateFeedbackInputLoop(out.physical, out.physicalGain),
                     o.frequencyPoints);
  out.surrogateLqg =
      loopRobustness(outputFeedbackInputLoop(out.surrogate, out.compensator),
                     o.frequencyPoints);
  out.physicalLqg =
      loopRobustness(outputFeedbackInputLoop(out.physical, out.compensator),
                     o.frequencyPoints);
  return out;
}
} // namespace dynsyslab::core
