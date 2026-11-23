#include "KromTrajectoryAudit.hpp"
#include <cmath>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KromDifferentialAnalysis.hpp>
#include <dynsyslab/core/KromLocalControlAnalysis.hpp>
#include <dynsyslab/core/LieAlgebra.hpp>
#include <dynsyslab/workflow/KromAuditQuery.hpp>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {
using namespace core;
using J = json::Value;
J number(double v) { return std::isfinite(v) ? J(v) : J(); }
J array(const std::vector<double> &values) {
  J::Array a;
  for (auto x : values)
    a.push_back(number(x));
  return a;
}
J complexArray(const std::vector<std::complex<double>> &values) {
  J::Array a;
  for (auto x : values)
    a.emplace_back(J::Array{number(x.real()), number(x.imag())});
  return a;
}
std::vector<double> numbers(const J &value) {
  if (!std::holds_alternative<J::Array>(value.storage()))
    throw std::invalid_argument("Expected numeric array");
  std::vector<double> out;
  for (const auto &v : value.array()) {
    if (!std::holds_alternative<double>(v.storage()) ||
        !std::isfinite(v.number()))
      throw std::invalid_argument("Expected finite numbers");
    out.push_back(v.number());
  }
  return out;
}
J polynomial(const lie::Polynomial &p) {
  J::Array a;
  for (const auto &t : p.terms) {
    J::Array powers;
    for (auto k : t.powers)
      powers.emplace_back(double(k));
    a.emplace_back(
        J::Object{{"coefficient", t.coefficient}, {"powers", powers}});
  }
  return a;
}
J structure(const StateSpaceModel &m) {
  return J::Object{{"states", double(m.states)},
                   {"inputs", double(m.inputs)},
                   {"outputs", double(m.outputs)},
                   {"A", array(m.A)},
                   {"B", array(m.B)},
                   {"C", array(m.C)},
                   {"D", array(m.D)}};
}
J margins(const LoopRobustness &m) {
  return J::Object{
      {"nominal_stable", m.closedLoopStable},
      {"gain_lower", number(m.gainMarginLower)},
      {"gain_upper", number(m.gainMarginUpper)},
      {"phase_margin_deg", number(m.phaseMargin)},
      {"minimum_return_difference", number(m.minimumReturnDifference)},
      {"sensitivity_peak", number(m.sensitivityPeak)},
      {"complementary_sensitivity_peak",
       number(m.complementarySensitivityPeak)},
      {"disk_gain_lower", number(m.diskGainLower)},
      {"disk_gain_upper", number(m.diskGainUpper)},
      {"disk_phase_margin_deg", number(m.diskPhaseMargin)},
      {"frequencies", array(m.frequencies)},
      {"response", complexArray(m.response)}};
}
J derivativeClosure(const KoopmanModel &model, const J &request) {
  const auto &bank = model.polynomialEigenfunctions();
  const auto &dict = bank.dictionary();
  const auto d = dict.coordinateDimension, features = dict.featureCount();
  const auto modes = bank.eigenfunctionCount(), split = model.splitDimension();
  std::vector<lie::Polynomial> seeds(split + 1, lie::Polynomial{d, {}});
  seeds[0].terms.push_back({std::vector<unsigned>(d), "1"});
  for (std::size_t k = 0; k < features; ++k) {
    std::vector<std::complex<double>> coefficient(modes);
    for (std::size_t i = 0; i < modes; ++i)
      coefficient[i] = bank.coefficients()[i * features + k];
    const auto c = model.splitValues(coefficient);
    std::vector<unsigned> powers;
    for (std::size_t j = 0; j < d; ++j)
      powers.push_back(unsigned(dict.exponents[k * d + j]));
    for (std::size_t i = 0; i < split; ++i)
      if (c[i] != 0)
        seeds[i + 1].terms.push_back({powers, lie::binaryRational(c[i])});
  }
  seeds = lie::independentPolynomials(seeds);
  J::Array observables;
  for (const auto &p : seeds)
    observables.push_back(polynomial(p));
  J::Object query{{"schema", "dynsyslab.lie-algebra-query/v1"},
                  {"dimension", double(d)},
                  {"action", "module"},
                  {"observables", observables},
                  {"operators", J::Array{}},
                  {"include_partials", true},
                  {"extend", request.at("extend").boolean(true)}};
  if (!request.at("limits").isNull())
    query["limits"] = request.at("limits");
  auto result = lie::analyze(query);
  result.at("coordinate_scope") =
      "Partial derivatives in standardized fitting coordinates. Equivalent "
      "closure in unstandardized fitting coordinates; not automatically "
      "physical-state/input/drift closure.";
  result.at("coefficient_scope") = "Exact binary floating coefficients of the "
                                   "stored approximate polynomial bank";
  result.at("initial_independent_dimension") = double(seeds.size());
  return result;
}
} // namespace
core::json::Value evaluateKromAuditQuery(const core::json::Value &request) {
  if (request.at("schema").string() != "dynsyslab.krom-audit-query/v1")
    throw std::invalid_argument("Unknown KROM audit schema");
  const auto action = request.at("action").string();
  J::Object out{{"schema", "dynsyslab.krom-audit-result/v1"},
                {"action", action}};
  if (action == "regress_action") {
    const auto n = request.at("features").number(),
               m = request.at("outputs").number();
    if (n < 1 || m < 1 || n > 256 || m > 4096 || std::floor(n) != n ||
        std::floor(m) != m)
      throw std::invalid_argument("Invalid action dimensions");
    auto fit = regressInvariantAction(
        numbers(request.at("training_features")),
        numbers(request.at("training_action")),
        numbers(request.at("validation_features")),
        numbers(request.at("validation_action")), std::size_t(n),
        std::size_t(m), request.at("tolerance").number(1e-10));
    out["coefficients"] = array(fit.coefficients);
    out["predictors"] = double(fit.predictors);
    out["outputs"] = double(fit.outputs);
    out["training_relative_error"] = number(fit.trainingRelativeError);
    out["validation_relative_error"] = number(fit.validationRelativeError);
    out["scope"] =
        "Sampled affine-span least squares; not an exact closure certificate";
    return out;
  }
  const auto path = std::filesystem::u8path(request.at("krom").string());
  const auto model = loadKoopmanModel(path);
  out["krom"] = path.generic_u8string();
  out["fitting_coordinates"] =
      eigenfunctionCoordinatesJson(model.coordinates());
  out["lifted_dimension"] = double(model.splitDimension());
  out["eigenvalues"] = complexArray(model.continuousEigenvalues());
  if (action == "trajectory")
    out["trajectory"] = auditKromTrajectory(model, request);
  else if (action == "derivative_closure")
    out["closure"] = derivativeClosure(model, request);
  else if (action == "differential") {
    const auto data =
        analyzeKromDifferentials(model, numbers(request.at("points")),
                                 request.at("input_map").string("exact"));
    out["samples"] = double(data.samples);
    out["state_dimension"] = double(data.states);
    out["values"] = array(data.values);
    out["jacobians"] = array(data.jacobians);
    out["drift"] = array(data.drift);
    out["generator_defect"] = array(data.generatorDefect);
    out["true_input_map"] = array(data.trueInputMap);
    out["fitted_input_map"] = array(data.fittedInputMap);
    out["input_defect"] = array(data.inputDefect);
    out["sigma_min"] = array(data.smallestSingularValue);
    out["sigma_max"] = array(data.largestSingularValue);
    out["determinant"] = array(data.orientedDeterminant);
    out["scope"] =
        "Pointwise physical-state defects and raw-coordinate Jacobian singular "
        "values; no global injectivity or uniform error certificate";
  } else if (action == "local_control") {
    KromLocalControlOptions o;
    o.point = numbers(request.at("point"));
    o.measurement = request.at("measurement").string();
    o.controlMap = request.at("input_map").string("constant");
    o.holdingInput = request.at("holding_input").number(0);
    const double count = request.at("frequency_points").number(1200);
    if (!std::isfinite(count) || count < 100 || count > 100000 ||
        std::floor(count) != count)
      throw std::invalid_argument("Invalid frequency point count");
    o.frequencyPoints = std::size_t(count);
    o.statePenalty = request.at("Q").number(5000);
    o.inputPenalty = request.at("R").number(.01);
    o.processCovariance = request.at("W").number(1);
    o.measurementCovariance = request.at("V").number(1e-4);
    auto a = analyzeKromLocalControl(model, o);
    out["surrogate"] = structure(a.surrogate);
    out["physical"] = structure(a.physical);
    out["compensator"] = structure(a.compensator);
    out["gain"] = array(a.gain);
    out["physical_gain"] = array(a.physicalGain);
    out["estimator_gain"] = array(a.estimatorGain);
    out["physical_equilibrium_residual"] = array(a.physicalEquilibriumResidual);
    out["lifted_at_point"] = array(a.liftedAtPoint);
    out["generator_at_point"] = array(a.generatorAtPoint);
    out["measured_at_point"] = array(a.measuredAtPoint);
    out["decoded_measurement_at_point"] = array(a.decodedMeasurementAtPoint);
    out["decoded_state_at_point"] = array(a.decodedStateAtPoint);
    out["surrogate_lqr"] = margins(a.surrogateLqr);
    out["physical_lqr"] = margins(a.physicalLqr);
    out["surrogate_lqg"] = margins(a.surrogateLqg);
    out["physical_lqg"] = margins(a.physicalLqg);
    out["scope"] =
        "Frozen continuous local loops; same compensator applied to surrogate "
        "and physical sensor/plant. Not a guarantee for tracking, scheduled "
        "gains, saturation, sampling, delay or noise.";
  } else
    throw std::invalid_argument("Unknown KROM audit action");
  return out;
}
} // namespace dynsyslab::workflow
