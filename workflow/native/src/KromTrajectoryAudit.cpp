#include "MeasurementMap.hpp"
#include "KromTrajectoryAudit.hpp"
#include "PostprocessInternal.hpp"
#include <cmath>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanStateSpace.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {
using J = core::json::Value;
J array(const std::vector<double> &v) {
  J::Array a;
  for (double x : v)
    a.emplace_back(std::isfinite(x) ? J(x) : J());
  return a;
}
} // namespace
// Compose the existing causal history, frozen KROM, exact ZOH discretization,
// linear simulation and reconstruction services. No fitting or future-state
// correction is used by the replay. The archived input is an exogenous input.
core::json::Value auditKromTrajectory(const core::KoopmanModel &model,
                                      const J &q) {
  auto source = postprocess_detail::loadSource(
      std::filesystem::u8path(q.at("source").string()));
  postprocess_detail::normalizeSourceCoordinates(source);
  core::requireCompatibleSystems(
      model.systemConfiguration(),
      postprocess_detail::systemConfiguration(source));
  const double dt = postprocess_detail::uniformDeltaTime(source);
  if (source.dimension() != model.stateDimension() || source.samples() < 2)
    throw std::invalid_argument(
        "Trajectory must contain physical states on a uniform clock");
  const double strideValue = q.at("stride").number(1);
  if (strideValue < 1 || strideValue > source.samples() ||
      std::floor(strideValue) != strideValue)
    throw std::invalid_argument("Invalid output stride");
  const auto stride = std::size_t(strideValue), n = model.splitDimension(),
             d = source.dimension();
  const auto map = q.at("input_map").string("constant");
  if (map != "constant" && map != "none")
    throw std::invalid_argument(
        "Replay requires a constant map, or none for unforced data");
  if (map == "none")
    for (auto u : source.inputs())
      if (u != 0)
        throw std::invalid_argument(
            "Cannot ignore nonzero input in autonomous replay");
  if (!source.inputs().empty() && source.inputDimension() != 1)
    throw std::invalid_argument("Replay currently requires scalar input");
  std::optional<core::CausalHistoryProjector> history;
  std::optional<core::Observable> measurement;
  if (model.requiresExternalCoordinates()) {
    auto transform = postprocess_detail::loadMeasurementMap(
        model, std::filesystem::u8path(q.at("krom").string()), q).transform;
    if (transform.retainedCoordinates().size() != model.coordinateDimension() ||
        std::abs(transform.sampleInterval() - dt) > 1e-8 * dt)
      throw std::invalid_argument(
          "Frozen history coordinates or clock differ from trajectory");
    auto system = core::createSystem(model.systemConfiguration());
    measurement = system->observableCatalog().at(
        q.at("measurement").string(transform.measurement()));
    if (measurement->outputDimension() != 1)
      throw std::invalid_argument("History audit requires scalar measurement");
    history.emplace(std::move(transform));
  }
  std::vector<double> time, states, lifted, input;
  std::vector<std::complex<double>> complexLifted;
  for (std::size_t s = 0; s < source.samples(); ++s) {
    const auto *x = source.values().data() + s * d;
    std::vector<std::complex<double>> z;
    if (history) {
      const auto c =
          history->push(source.time()[s], measurement->evaluateReal(x).front());
      if (!c)
        continue;
      z = model.evaluateCoordinates(c->data());
    } else
      z = model.evaluate(x);
    const auto real = model.splitValues(z);
    time.push_back(source.time()[s]);
    states.insert(states.end(), x, x + d);
    lifted.insert(lifted.end(), real.begin(), real.end());
    complexLifted.insert(complexLifted.end(), z.begin(), z.end());
    input.push_back(source.inputs().empty() ? 0 : source.inputs()[s]);
  }
  if (time.size() < 2)
    throw std::invalid_argument("Insufficient complete causal history windows");
  const auto b =
      map == "none"
          ? std::vector<double>(n)
          : model.realControlMap(model.controlMapFromLifted(
                std::vector<std::complex<double>>(model.liftedDimension()),
                map));
  const auto linear =
      core::makeStateSpace(model.generatorSplit(), b, {}, {}, n, 1, 0);
  const auto discrete = core::discretize(linear, dt);
  const auto replay = core::simulateLinear(
      linear, time, input,
      std::vector<double>(lifted.begin(), lifted.begin() + n));
  std::vector<double> oneStep, defect;
  double error = 0, change = 0;
  for (std::size_t s = 0; s < time.size(); ++s)
    for (std::size_t i = 0; i < n; ++i) {
      double value = lifted[i];
      if (s) {
        value = discrete.B[i] * input[s - 1];
        for (std::size_t j = 0; j < n; ++j)
          value += discrete.A[i * n + j] * lifted[(s - 1) * n + j];
      }
      oneStep.push_back(value);
      const double delta = lifted[s * n + i] - value;
      defect.push_back(delta / dt);
      if (s) {
        error += delta * delta;
        const double c = lifted[s * n + i] - lifted[(s - 1) * n + i];
        change += c * c;
      }
    }
  std::vector<std::complex<double>> replayComplex;
  for (std::size_t s = 0; s < time.size(); ++s) {
    auto c = model.valuesFromSplit(std::vector<double>(
        replay.states.begin() + s * n, replay.states.begin() + (s + 1) * n));
    replayComplex.insert(replayComplex.end(), c.begin(), c.end());
  }
  J::Object out{
      {"samples", double(time.size())},
      {"lifted_dimension", double(n)},
      {"state_dimension", double(d)},
      {"step_defect_relative_to_change",
       change > 0 ? J(std::sqrt(error / change)) : J()},
      {"scope", "Measured lift/readout versus one-step ZOH and input-replay "
                "predictions. Replay uses only initial lift and recorded "
                "external input; it is not closed-loop prediction. Output "
                "stride never changes the history or simulation clock."}};
  auto save = [&](const std::string &key, const std::vector<double> &values,
                  std::size_t width) {
    std::vector<double> selected;
    for (std::size_t s = 0; s < time.size(); s += stride)
      selected.insert(selected.end(), values.begin() + s * width,
                      values.begin() + (s + 1) * width);
    out[key] = array(selected);
  };
  save("time", time, 1);
  save("states", states, d);
  save("input", input, 1);
  save("lifted", lifted, n);
  save("one_step", oneStep, n);
  save("step_defect_per_time", defect, n);
  save("replay_lifted", replay.states, n);
  for (const auto &name :
       {std::string("state"), std::string("far_field_acoustic_echo")})
    if (model.hasReconstruction(name)) {
      const auto &readout = model.reconstruction(name);
      const auto decoded =
          core::reconstructObservable(readout, complexLifted, time.size());
      const auto predicted =
          core::reconstructObservable(readout, replayComplex, time.size());
      std::vector<double> y, p;
      for (auto c : decoded)
        y.push_back(c.real());
      for (auto c : predicted)
        p.push_back(c.real());
      save("decode_" + name, y, readout.outputs);
      save("replay_" + name, p, readout.outputs);
    }
  return out;
}
} // namespace dynsyslab::workflow
