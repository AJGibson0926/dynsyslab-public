// Frozen-A bilinear K-ROMs through input_map, the certificate step, KLQE and KLQR.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/BilinearKrom.hpp>
#include <dynsyslab/core/EstimationRuntime.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanBilinearBridge.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/SupportedResidualSector.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
namespace w = dynsyslab::workflow;
namespace json = c::json;
namespace fs = std::filesystem;
using Complex = std::complex<double>;

int checks = 0;
void check(bool ok, const std::string& why) {
    ++checks;
    if (!ok) throw std::runtime_error(why);
}
template <class F> void rejects(F f, const std::string& why) {
    ++checks;
    try { f(); } catch (const std::logic_error&) { return; }
    throw std::runtime_error("accepted: " + why);
}
std::string bytes(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
c::SystemConfiguration duffing(double beta) {
    c::SystemConfiguration config; config.identifier = "duffing";
    config.parameters = {{"alpha", 2}, {"beta", beta}, {"delta", 1}};
    return config;
}

void frozenFit() {
    const std::vector<double> A{0.9, 0.1, -0.2, 0.8}, B0{0.3, -0.1}, B1{0.05, 0.0, 0.02, -0.04};
    c::BilinearKrom truth; truth.stateDimension = 2; truth.inputDimension = 1;
    truth.A = A; truth.B0 = B0; truth.B = {B1};
    std::mt19937 random(5); std::uniform_real_distribution<double> draw(-1, 1);
    c::BilinearDataset data;
    for (int k = 0; k < 40; ++k) {
        const std::vector<double> z{draw(random), draw(random)}, u{draw(random)};
        const auto next = c::stepBilinearKrom(truth, z, u);
        data.z.insert(data.z.end(), z.begin(), z.end());
        data.zNext.insert(data.zNext.end(), next.begin(), next.end());
        data.u.push_back(u[0]); ++data.samples;
    }
    const auto fit = c::fitBilinearKromFrozenA(data, A, 2, 1);
    check(fit.model.A == A, "frozen fit changed A");
    check(fit.identifiable && fit.regressors == 3 && fit.rank == 3, "frozen fit diagnostics");
    for (std::size_t k = 0; k < 2; ++k) check(std::abs(fit.model.B0[k] - B0[k]) < 1e-11, "frozen fit B0");
    for (std::size_t k = 0; k < 4; ++k) check(std::abs(fit.model.B[0][k] - B1[k]) < 1e-11, "frozen fit B1");
    check(fit.stats.maximum < 1e-11, "frozen fit residual");
    // A wrong frozen A is not repaired by the fit: the error stays in the residual.
    auto wrong = A; wrong[0] += 0.05;
    const auto biased = c::fitBilinearKromFrozenA(data, wrong, 2, 1);
    check(biased.model.A == wrong && biased.stats.maximum > 1e-3, "frozen fit absorbed an error in A");
}

// psi = x + i v with lambda = -0.5 + 2i: one complex mode, two split coordinates.
c::KoopmanModel complexModel(const c::SystemConfiguration& config, Complex n, Complex conjugateN, Complex b0) {
    auto model = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1},
        {Complex(1, 0), Complex(0, 1)}, {Complex(-0.5, 2.0)});
    c::BilinearControlMap map;
    map.modes = 1; map.inputs = 1; map.operators = {n}; map.conjugateOperators = {conjugateN}; map.offsets = {b0};
    map.provenance = "analytic test map";
    model.setControlMap("bilinear", map);
    return model;
}

void bridgeConversion() {
    const auto config = duffing(0.0);
    const auto model = complexModel(config, {0.3, -0.2}, {0.1, 0.4}, {1.0, 0.5});
    const double dt = 0.05;
    c::KoopmanBilinearBridgeOptions options; options.sampleTime = dt;
    const auto bridge = c::bridgeKoopmanBilinear(model, options);
    check(bridge.model.A == bridge.frozenTransition, "A differs from the transition");
    check(bridge.frozenTransition == dynsyslab::numerics::matrixExponential(model.generatorSplit(), 2, dt),
          "A is not the exact exponential of the frozen generator");
    check(bridge.targetIsEquilibrium && bridge.constantNorm == 0.0, "zero target is an equilibrium");
    // Real-split operator of n psi + c conj(psi).
    const auto& N = bridge.continuous.operators[0];
    check(std::abs(N[0] - 0.4) < 1e-14 && std::abs(N[1] - (0.2 + 0.4)) < 1e-14
          && std::abs(N[2] - (-0.2 + 0.4)) < 1e-14 && std::abs(N[3] - (0.3 - 0.1)) < 1e-14, "split operator");
    // B_1 is the exact first-order term: the defect shrinks as u^2.
    const auto large = c::conversionDefect(bridge, {0.2}), small = c::conversionDefect(bridge, {0.1});
    check(large.transition > 0 && std::abs(large.transition / small.transition - 4.0) < 0.1, "B_1 not first-order exact");

    // A nonzero target at u = 0 is not an equilibrium: the constant is reported.
    options.targetLifted = {0.2, -0.1};
    const auto offset = c::bridgeKoopmanBilinear(model, options);
    check(!offset.targetIsEquilibrium && offset.constantNorm > 1e-3, "constant term hidden");
    check(offset.model.A == bridge.frozenTransition, "centring changed A");

    // Held u*: A_eff = A + u* B_1, and a constant is exactly A_eff w* + B0 u* - w*.
    options.targetInput = {0.3};
    const auto held = c::bridgeKoopmanBilinear(model, options);
    check(held.heldTargetInput, "held input flag");
    for (std::size_t k = 0; k < 4; ++k)
        check(held.model.A[k] == bridge.frozenTransition[k] + 0.3 * bridge.discreteOperators[0][k], "A_eff");
    for (std::size_t r = 0; r < 2; ++r) {
        double e = 0.3 * bridge.discreteOffsets[r] - options.targetLifted[r];
        for (std::size_t j = 0; j < 2; ++j) e += held.model.A[r * 2 + j] * options.targetLifted[j];
        check(std::abs(e - held.constantTerm[r]) < 1e-15, "constant term");
    }
    c::BilinearDataset empty;
    rejects([&] { (void)c::refitBridgeOperators(held, empty); }, "operator refit with a held target input");
}

// A forced record whose values follow the exact held-input bilinear flow.
fs::path syntheticRecord(const fs::path& root, const fs::path& template_, const c::KoopmanModel& truth, double amplitude) {
    auto archive = c::readNpz(template_);
    const auto time = archive.at("time").doubles();
    const double dt = time[1] - time[0];
    const auto split = c::continuousBilinearSplit(truth, "bilinear");
    std::mt19937 random(11); std::uniform_real_distribution<double> draw(-1, 1);
    std::vector<double> values(time.size() * 2), inputs(time.size());
    std::vector<double> z{0.6, -0.4};
    for (std::size_t k = 0; k < time.size(); ++k) {
        if (k % 25 == 0) { z = {draw(random), draw(random)}; }
        inputs[k] = amplitude * draw(random);
        values[k * 2] = z[0]; values[k * 2 + 1] = z[1];
        auto generator = split.generator;
        for (std::size_t e = 0; e < 4; ++e) generator[e] += inputs[k] * split.operators[0][e];
        const auto held = dynsyslab::numerics::zeroOrderHold(generator, split.offsets, 2, 1, dt);
        const std::vector<double> before = z;
        for (std::size_t r = 0; r < 2; ++r)
            z[r] = held.transition[r * 2] * before[0] + held.transition[r * 2 + 1] * before[1] + held.response[r] * inputs[k];
    }
    // Each block restart is a jump, so those pairs must not enter a one-step regression.
    archive.arrays["values"] = c::NpyArray::fromDoubles({time.size(), 2}, values);
    std::vector<double> interval(inputs.begin(), inputs.end() - 1);
    archive.arrays["control_interval_values"] = c::NpyArray::fromDoubles({time.size() - 1, 1}, interval);
    archive.arrays["control_values"] = c::NpyArray::fromDoubles({time.size(), 1}, inputs);
    const auto path = root / "synthetic-forced.npz";
    c::writeNpzAtomic(path, archive);
    return path;
}

json::Value postJob(const char* method, json::Value settings, const std::vector<fs::path>& sources) {
    json::Value::Array list;
    for (const auto& s : sources) list.emplace_back(s.generic_u8string());
    return json::Value::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"},
        {"label", std::string("bilinear-") + method}, {"sources", std::move(list)},
        {"methods", json::Value::Array{json::Value::Object{{"name", method}, {"settings", std::move(settings)}}}}};
}

fs::path forcedDuffing(const fs::path& root, double beta, double amplitude) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"bilinear-forced",
      "system":{"name":"duffing","parameters":{"alpha":2,"delta":1}},
      "time":{"start":0,"stop":4,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","values":[[0.04,0.02]],"count":1,"clusters":1},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_evaluation":"sampled","control_sample_time":0.02},
      "controller":{"kind":"random_pulse","settings":{"model_pulse_duration":0.1,"seed":7}},
      "archive":{"enabled":true}})");
    job.at("system").at("parameters").object()["beta"] = beta;
    job.at("controller").at("settings").object()["model_amplitude"] = amplitude;
    return w::executeJob(job, root).archives.front();
}

// Short trajectories from a grid of lifted states: a zero-input hold, then random pulses,
// so both residual channels are sampled across the same ||z|| range.
std::vector<fs::path> forcedWithHolds(const fs::path& root) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"bilinear-holds",
      "system":{"name":"duffing","parameters":{"alpha":2,"beta":0.2,"delta":1}},
      "time":{"start":0,"stop":0.6,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","count":128,"clusters":1,"values":[]},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_evaluation":"sampled","control_sample_time":0.02},
      "controller":{"kind":"random_pulse","settings":{"model_amplitude":0.75,"model_pulse_duration":0.06,"seed":7,
        "active_time_window":[0.3,0.6]}},
      "archive":{"enabled":true}})");
    json::Value::Array states;
    for (double radius : {0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5})
        for (int angle = 0; angle < 16; ++angle) {
            const double t = 2 * 3.141592653589793 * angle / 16, a = radius * std::cos(t), b = radius * std::sin(t);
            states.push_back(json::Value::Array{(a - b) / 3, (a + 2 * b) / 3});
        }
    job.at("initial_conditions").at("values") = std::move(states);
    return w::executeJob(job, root).archives;
}

// One exact-enough held step of x'' = 2x - 0.2x^3 - x' + u (RK4), independent of the fit.
std::pair<double, double> heldDuffingStep(double x, double v, double u, double dt) {
    const auto f = [u](double a, double b) { return 2 * a - 0.2 * a * a * a - b + u; };
    const int steps = 40; const double h = dt / steps;
    for (int i = 0; i < steps; ++i) {
        const double k1x = v, k1v = f(x, v);
        const double k2x = v + h / 2 * k1v, k2v = f(x + h / 2 * k1x, v + h / 2 * k1v);
        const double k3x = v + h / 2 * k2v, k3v = f(x + h / 2 * k2x, v + h / 2 * k2v);
        const double k4x = v + h * k3v, k4v = f(x + h * k3x, v + h * k3v);
        x += h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x);
        v += h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);
    }
    return {x, v};
}

// Held-out grid inside the certified ellipsoid and input box, never used by the fit.
struct HeldOutResult { std::size_t points{}, violations{}; double worstRatio{}, worstZeroInputRatio{}; };
HeldOutResult heldOutSectorCheck(const c::NpzRecord& archive, const json::Value& record, double dt) {
    c::BilinearKrom model; model.stateDimension = 2; model.inputDimension = 1; model.sampleTime = dt;
    model.A = archive.at("model_A").doubles(); model.B0 = archive.at("model_B0").doubles();
    model.B = {archive.at("model_B").doubles()};
    const auto X = archive.at("certificate_X").doubles();
    const double l11 = std::sqrt(X[0]), l21 = X[2] / l11, l22 = std::sqrt(X[3] - l21 * l21);
    const double cx = record.at("sector").at("c_x").number(), cu = record.at("sector").at("c_u").number();
    const double ubar = record.at("input_bounds").array().front().number();
    HeldOutResult out;
    for (double scale : {0.2, 0.4, 0.6, 0.8, 1.0})
        for (int angle = 0; angle < 48; ++angle)
            for (double fraction : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
                const double t = 2 * 3.141592653589793 * (angle + 0.37) / 48, a = scale * std::cos(t), b = scale * std::sin(t);
                const std::vector<double> z{l11 * a, l21 * a + l22 * b}, u{fraction * ubar};
                const double x = (z[0] - z[1]) / 3, v = (z[0] + 2 * z[1]) / 3;
                const auto [xp, vp] = heldDuffingStep(x, v, u[0], dt);
                const auto predicted = c::stepBilinearKrom(model, z, u);
                const double r = std::hypot(2 * xp + vp - predicted[0], -xp + vp - predicted[1]);
                const double bound = cx * std::hypot(z[0], z[1]) + cu * std::abs(u[0]);
                const double ratio = r / bound;
                ++out.points;
                if (r > bound) ++out.violations;
                out.worstRatio = std::max(out.worstRatio, ratio);
                if (fraction == 0.0) out.worstZeroInputRatio = std::max(out.worstZeroInputRatio, ratio);
            }
    return out;
}

void supportedSectorRules() {
    c::BilinearKrom model; model.stateDimension = 1; model.inputDimension = 1;
    model.A = {0.5}; model.B0 = {1.0}; model.B = {{0.0}};
    const auto build = [&](bool holds, double cubic) {
        c::BilinearDataset d;
        for (int k = 1; k <= 40; ++k) {
            const double z = 0.025 * k, u = holds && k % 2 ? 0.0 : (k % 4 ? 0.7 : -0.7);
            const std::vector<double> next = c::stepBilinearKrom(model, {z}, {u});
            d.z.push_back(z); d.u.push_back(u); d.zNext.push_back(next[0] + cubic * z * z * z + 0.01 * u); ++d.samples;
        }
        return d;
    };
    const auto inputOnly = c::estimateSupportedResidualSector(model, build(false, 0.1));
    check(!inputOnly.supported && inputOnly.reason.find("c_x is unsupported") != std::string::npos, "input-only data accepted");
    const auto both = c::estimateSupportedResidualSector(model, build(true, 0.1));
    check(both.supported && both.bound.cState > 0.09 && both.bound.cInput > 0.0, "two-stage sector");
    c::SupportedSectorOptions declared; declared.oneSided = c::OneSidedSectorDeclaration{"state", "known input-only error"};
    const auto contradicted = c::estimateSupportedResidualSector(model, build(true, 0.1), declared);
    check(!contradicted.supported, "a contradicted one-sided declaration was accepted");
    const auto honest = c::estimateSupportedResidualSector(model, build(false, 0.0), declared);
    check(honest.supported && honest.bound.cState == 0.0 && honest.declarationReason == "known input-only error", "declared one-sided sector");
    declared.oneSided->reason.clear();
    rejects([&] { (void)c::estimateSupportedResidualSector(model, build(false, 0.0), declared); }, "declaration without a reason");
}

void inputMapKeepsA(const fs::path& root) {
    const auto config = duffing(0.0);
    const Complex n(0.3, -0.2), conjugateN(0.1, 0.4), b0(1.0, 0.5);
    const auto truth = complexModel(config, n, conjugateN, b0);
    auto bare = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1},
        {Complex(1, 0), Complex(0, 1)}, {Complex(-0.5, 2.0)});
    const auto kromPath = root / "complex-krom.npz";
    c::savePolynomialKoopmanModel(kromPath, bare);
    const auto kromBytes = bytes(kromPath);
    const auto record = syntheticRecord(root, forcedDuffing(root / "forced-template", 0.0, 0.2), truth, 0.01);

    json::Value::Object settings{{"krom", kromPath.generic_u8string()}, {"representation", "bilinear"},
        {"fit", "held_input"}, {"conjugate_complete", true}};
    const auto output = w::executePostprocess(postJob("input_map", settings, {record}), root).archives.front();
    check(bytes(kromPath) == kromBytes, "input_map modified the source K-ROM");
    const auto fitted = c::loadKoopmanModel(output);
    check(fitted.generatorSplit() == bare.generatorSplit(), "input_map changed the autonomous generator");
    check(fitted.continuousEigenvalues() == bare.continuousEigenvalues(), "input_map changed the eigenvalues");
    const auto& map = std::get<c::BilinearControlMap>(fitted.storedControlMap("bilinear"));
    check(map.conjugateOperators.size() == 1, "conjugate operator not stored");
    // Blocks restart from random states, so a few pairs are jumps; exclude them with a test of scale only.
}

void cleanRecovery(const fs::path& root) {
    const auto config = duffing(0.0);
    const Complex n(0.3, -0.2), conjugateN(0.1, 0.4), b0(1.0, 0.5);
    const auto truth = complexModel(config, n, conjugateN, b0);
    auto bare = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1},
        {Complex(1, 0), Complex(0, 1)}, {Complex(-0.5, 2.0)});
    const auto kromPath = root / "clean-krom.npz";
    c::savePolynomialKoopmanModel(kromPath, bare);
    auto archive = c::readNpz(forcedDuffing(root / "clean-template", 0.0, 0.2));
    const auto time = archive.at("time").doubles();
    const double dt = time[1] - time[0];
    const auto split = c::continuousBilinearSplit(truth, "bilinear");
    std::mt19937 random(3); std::uniform_real_distribution<double> draw(-1, 1);
    std::vector<double> values(time.size() * 2), inputs(time.size() - 1);
    std::vector<double> z{0.8, -0.5};
    for (std::size_t k = 0; k < time.size(); ++k) {
        values[k * 2] = z[0]; values[k * 2 + 1] = z[1];
        if (k + 1 == time.size()) break;
        // Small inputs make the O(u^2) part of the held flow negligible; a weak restoring
        // excitation keeps the state away from the origin.
        inputs[k] = 1e-3 * draw(random);
        auto generator = split.generator;
        for (std::size_t e = 0; e < 4; ++e) generator[e] += inputs[k] * split.operators[0][e];
        const auto held = dynsyslab::numerics::zeroOrderHold(generator, split.offsets, 2, 1, dt);
        const auto before = z;
        for (std::size_t r = 0; r < 2; ++r)
            z[r] = held.transition[r * 2] * before[0] + held.transition[r * 2 + 1] * before[1] + held.response[r] * inputs[k];
    }
    archive.arrays["values"] = c::NpyArray::fromDoubles({time.size(), 2}, values);
    archive.arrays["control_interval_values"] = c::NpyArray::fromDoubles({inputs.size(), 1}, inputs);
    archive.arrays.erase("control_values");
    const auto path = root / "clean-forced.npz";
    c::writeNpzAtomic(path, archive);
    json::Value::Object settings{{"krom", kromPath.generic_u8string()}, {"representation", "bilinear"},
        {"fit", "held_input"}, {"conjugate_complete", true}};
    const auto output = w::executePostprocess(postJob("input_map", settings, {path}), root).archives.front();
    const auto fitted = c::loadKoopmanModel(output);
    check(fitted.generatorSplit() == bare.generatorSplit(), "input_map changed the generator");
    const auto& map = std::get<c::BilinearControlMap>(fitted.storedControlMap("bilinear"));
    std::cout << "clean n=" << map.operators[0] << " c=" << map.conjugateOperators[0] << " b0=" << map.offsets[0] << '\n';
    check(std::abs(map.offsets[0] - b0) < 1e-3, "held-input b0 not recovered");
    check(std::abs(map.operators[0] - n) < 5e-3, "held-input operator not recovered");
    check(std::abs(map.conjugateOperators[0] - conjugateN) < 5e-3, "held-input conjugate operator not recovered");
}

// Causal-history K-ROM: psi = newest delay of the measured position, with a scalar exact
// bilinear flow. The operator kernel phi(lambda, lambda) = dt exp(lambda dt) differs from
// the offset kernel at order dt, so a 1e-3 tolerance detects the older conversion.
void historyHeldInput(const fs::path& root) {
    const auto config = duffing(0.0);
    const double lambda = -1.0, b0 = 0.8, n = 0.5;
    auto archive = c::readNpz(forcedDuffing(root / "history-template", 0.0, 0.2));
    const auto time = archive.at("time").doubles();
    const double dt = time[1] - time[0];
    c::PooledHankelBasis basis; basis.delays = 2; basis.rank = 2;
    basis.leftSingularVectors = {1, 0, 0, 1}; basis.singularValues = {1, 1};
    c::FrozenHistoryTransform transform(basis, dt, {0, 1}, "position");
    c::EigenfunctionCoordinates q; q.dimension = 2; q.kind = c::EigenfunctionCoordinateKind::causalHistory;
    q.providerIdentity = transform.identity();
    c::PolynomialEigenfunctionSet bank({2, {1, 0, 0, 1}, {0, 0}, {1, 1}}, {0., 1.},
        c::EigenfunctionSpectrum::fromContinuousRates({lambda}, dt), q);
    c::KoopmanModel model(bank, config);
    c::NpzRecord extra; c::writeFrozenHistoryTransform(extra, transform);
    const auto kromPath = root / "history-krom.npz";
    c::saveKoopmanModel(kromPath, model, "{}", extra.arrays);

    std::mt19937 random(9); std::uniform_real_distribution<double> draw(-1, 1);
    std::vector<double> values(time.size() * 2, 0.0), inputs(time.size() - 1);
    double y = 1.0;
    for (std::size_t k = 0; k < time.size(); ++k) {
        values[k * 2] = y;
        if (k + 1 == time.size()) break;
        inputs[k] = 1e-3 * draw(random);
        const double rate = lambda + n * inputs[k], flow = std::exp(rate * dt);
        y = flow * y + b0 * inputs[k] * (flow - 1.0) / rate;
    }
    archive.arrays["values"] = c::NpyArray::fromDoubles({time.size(), 2}, values);
    archive.arrays["control_interval_values"] = c::NpyArray::fromDoubles({inputs.size(), 1}, inputs);
    archive.arrays.erase("control_values");
    const auto record = root / "history-forced.npz";
    c::writeNpzAtomic(record, archive);
    json::Value::Object settings{{"krom", kromPath.generic_u8string()}, {"representation", "bilinear"}};
    const auto output = w::executePostprocess(postJob("input_map", settings, {record}), root).archives.front();
    const auto fitted = c::loadKoopmanModel(output);
    check(fitted.spectrum().continuousRates() == model.spectrum().continuousRates(), "history input_map changed the rates");
    const auto& map = std::get<c::BilinearControlMap>(fitted.storedControlMap("bilinear"));
    std::cout << "history b0=" << map.offsets[0] << " n=" << map.operators[0] << '\n';
    check(std::abs(map.offsets[0] - b0) < 1e-3, "history offset not recovered");
    check(std::abs(map.operators[0] - n) < 1e-3, "history operator not recovered by the first-order kernel");
}

void estimatorAndControllerUseBilinear() {
    c::SystemConfiguration config = duffing(0.0);
    auto model = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {1., -2.});
    model.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
    model.setReconstruction("state", {{1. / 3, -1. / 3, 1. / 3, 2. / 3}, 2, 2, 0, 0, 0});
    c::BilinearControlMap map;
    map.modes = 2; map.inputs = 1; map.offsets = {1.0, 1.0}; map.operators = {0.4, 0.1, -0.2, 0.3};
    map.provenance = "analytic test map";
    model.setControlMap("bilinear", map);

    // KLQE: on a stable spectrum, negligible process noise makes the Kalman gain vanish, so the
    // prediction is exp((G + u N) dt) z + ZOH(b0) u. (An unstable mode would keep a nonzero gain.)
    auto stable = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {-1., -2.});
    stable.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
    stable.setReconstruction("state", {{1. / 3, -1. / 3, 1. / 3, 2. / 3}, 2, 2, 0, 0, 0});
    stable.setControlMap("bilinear", map);
    const auto design = c::designKoopmanLQEFromMeasurements(stable, "sensor", {1e-8, 0, 0, 1e-8}, {1e8});
    c::KoopmanLQEObserver observer(stable, design, "bilinear", {0.5, -0.25});
    (void)observer.observe(0.0, {0.25}, {0.7});
    const auto next = observer.observe(0.05, {0.25}, {0.7});
    auto generator = stable.generatorSplit();
    const auto split = c::continuousBilinearSplit(stable, "bilinear");
    for (std::size_t e = 0; e < 4; ++e) generator[e] += 0.7 * split.operators[0][e];
    const auto held = dynsyslab::numerics::zeroOrderHold(generator, split.offsets, 2, 1, 0.05);
    for (std::size_t r = 0; r < 2; ++r) {
        const double exact = held.transition[r * 2] * 0.5 + held.transition[r * 2 + 1] * -0.25 + held.response[r] * 0.7;
        check(std::abs(next.lifted[r] - exact) < 1e-6, "KLQE does not propagate the bilinear term");
    }

    // KLQR: gain is scheduled on B(z) = b0 + N z.
    auto shared = std::make_shared<const c::KoopmanModel>(model);
    const std::vector<double> Q{1, 0, 0, 1}, R{1};
    c::KoopmanLQRController scheduled(shared, Q, R, {}, "bilinear", c::KoopmanTargetSpace::realLiftedValues);
    const std::vector<Complex> lifted{0.5, -0.25};
    const double u = scheduled.valueFromLifted(0.0, lifted).front();
    const std::vector<double> B{1.0 + 0.4 * 0.5 + 0.1 * -0.25, 1.0 - 0.2 * 0.5 + 0.3 * -0.25};
    const auto solved = c::solveContinuousLQR(model.generatorSplit(), B, Q, R, 2, 1);
    check(std::abs(u + solved.gain[0] * 0.5 + solved.gain[1] * -0.25) < 1e-9, "KLQR ignores the bilinear map");
    check(scheduled.gainSolveCount() == 1, "KLQR gain solve count");

    // Certified gain: u = u* + K (z - z*), no Riccati solve, and it refuses augmentation.
    c::KoopmanLQRController certified(shared, Q, R, {}, "bilinear", c::KoopmanTargetSpace::realLiftedValues);
    certified.setCertifiedFeedback({{-2.0, 0.5}, {0.1, -0.05}, {0.2}, 0.02});
    const double v = certified.valueFromLifted(0.0, lifted).front();
    check(std::abs(v - (0.2 - 2.0 * (0.5 - 0.1) + 0.5 * (-0.25 + 0.05))) < 1e-14, "certified sign or target input");
    check(certified.gainSolveCount() == 0, "certified gain ran a Riccati solve");
    rejects([&] { certified.setBalancedReduction({}); }, "balanced truncation on a certified gain");
    c::KoopmanLQRController moved(shared, Q, R, [](double) { return std::vector<double>{0.9, 0.9}; },
                                  "bilinear", c::KoopmanTargetSpace::realLiftedValues);
    moved.setCertifiedFeedback({{-2.0, 0.5}, {0.1, -0.05}, {0.2}, 0.02});
    rejects([&] { (void)moved.valueFromLifted(0.0, lifted); }, "controller target that leaves the certificate");
}

json::Value closedLoopJob(const fs::path& krom, const fs::path& certificate, double x0) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"certified-klqr",
      "system":{"name":"duffing","parameters":{"alpha":2,"beta":0.2,"delta":1}},
      "time":{"start":0,"stop":10,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","values":[[0,0]],"count":1,"clusters":1},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_evaluation":"sampled","control_sample_time":0.02},
      "controller":{"kind":"klqr","target":{"space":"state","kind":"constant","state":[0,0]},
        "settings":{"control_representation":"bilinear","Q":1,"R":1}},
      "archive":{"enabled":true}})");
    job.at("initial_conditions").at("values") = json::Value::Array{json::Value::Array{x0, 0.0}};
    auto& settings = job.at("controller").at("settings").object();
    settings["krom"] = krom.generic_u8string();
    settings["bilinear_certificate"] = certificate.generic_u8string();
    return job;
}

void endToEnd(const fs::path& root) {
    const auto config = duffing(0.2);
    auto krom = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {1., -2.});
    krom.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
    krom.setReconstruction("state", {{1. / 3, -1. / 3, 1. / 3, 2. / 3}, 2, 2, 0, 0, 0});
    const auto kromPath = root / "duffing-krom.npz";
    c::savePolynomialKoopmanModel(kromPath, krom);
    const auto inputOnly = forcedDuffing(root / "forced", 0.2, 1.0);
    const auto forced = forcedWithHolds(root / "forced-holds");
    check(forced.size() == 128, "forced records with holds");

    json::Value::Object mapSettings{{"krom", kromPath.generic_u8string()}, {"representation", "bilinear"}, {"fit", "held_input"}};
    const auto mapped = w::executePostprocess(postJob("input_map", mapSettings, forced), root).archives.front();
    const auto mappedBytes = bytes(mapped);
    check(c::loadKoopmanModel(mapped).generatorSplit() == krom.generatorSplit(), "input_map changed A");

    const auto certify = [&](json::Value::Object extra, std::vector<fs::path> records = {}) {
        if (records.empty()) records = forced;
        json::Value::Object s{{"krom", mapped.generic_u8string()}, {"target", json::Value::Object{{"state", json::Value::Array{0.0, 0.0}}}},
            {"input_bounds", json::Value::Array{0.75}}, {"state_radius", 0.5}, {"contraction", 0.99}, {"inflation", 1.25}, {"tolerance", 1e-5}};
        for (auto& [k, v] : extra) s[k] = std::move(v);
        return w::executePostprocess(postJob("bilinear_certificate", s, records), root).archives.front();
    };
    const auto path = certify({});
    check(bytes(mapped) == mappedBytes, "certificate modified the K-ROM");
    const auto archive = c::readNpz(path);
    const auto meta = json::parse(archive.metadataJson());
    const auto& record = meta.at("bilinear_certificate");
    std::cout << "certificate status=" << record.at("status").string()
              << " c_x=" << record.at("sector").at("c_x").number()
              << " c_u=" << record.at("sector").at("c_u").number()
              << " inscribed=" << record.at("inscribed_radius").number(0) << '\n';
    check(record.at("status").string() == "certified", "Duffing certificate failed: " + record.at("message").string());
    check(record.at("scope").string().find("not a guarantee for unseen states") != std::string::npos, "scope statement missing");
    const auto exact = dynsyslab::numerics::matrixExponential(krom.generatorSplit(), 2, 0.02);
    check(archive.at("frozen_A").doubles() == exact && archive.at("model_A").doubles() == exact, "certificate altered A");

    const auto heldOut = heldOutSectorCheck(archive, record, 0.02);
    std::cout << "held-out grid: points=" << heldOut.points << " violations=" << heldOut.violations
              << " worst r/bound=" << heldOut.worstRatio << " worst zero-input r/bound=" << heldOut.worstZeroInputRatio << '\n';
    check(record.at("sector").at("c_x").number() > 0.0 && record.at("sector").at("c_u").number() > 0.0, "sector is one-sided");
    check(heldOut.violations == 0, "the sector does not bound the true residual on the held-out grid");

    // An input box the records do not reach is refused rather than extrapolated.
    const auto wideBox = c::readNpz(certify({{"input_bounds", json::Value::Array{1.6}}}));
    check(json::parse(wideBox.metadataJson()).at("bilinear_certificate").at("status").string() == "sector_unsupported"
          && !wideBox.contains("certificate_K"), "input box beyond the data was certified");

    // Input-only forced records cannot identify c_x, so they yield no gain and cannot be deployed.
    const auto unsupportedPath = certify({}, {inputOnly});
    const auto unsupported = c::readNpz(unsupportedPath);
    const auto unsupportedMeta = json::parse(unsupported.metadataJson());
    check(unsupportedMeta.at("bilinear_certificate").at("status").string() == "sector_unsupported", "input-only data was not refused");
    check(!unsupportedMeta.at("bilinear_certificate").at("certified").boolean(true) && !unsupported.contains("certificate_K"), "unsupported archive carries a gain");
    auto unsupportedRun = closedLoopJob(mapped, unsupportedPath, 0.01);
    rejects([&] { (void)w::executeJob(unsupportedRun, root / "unsupported"); }, "KLQR deployed a sector_unsupported archive");

    // Closed loop from inside the certified ellipsoid.
    const double inscribed = record.at("inscribed_radius").number();
    const double x0 = 0.4 * inscribed / std::sqrt(5.0);
    const auto run = c::readNpz(w::executeJob(closedLoopJob(mapped, path, x0), root).archives.front());
    const auto x = run.at("values").doubles();
    check(std::hypot(x[x.size() - 2], x[x.size() - 1]) < 1e-3 * std::max(x0, 1e-3) + 1e-6, "certified KLQR did not settle");
    const auto runMeta = json::parse(run.metadataJson());
    const auto evidence = json::dump(runMeta);
    check(evidence.find("not a guarantee for unseen states") != std::string::npos, "run metadata lacks the certificate scope");
    const auto K = archive.at("certificate_K").doubles();
    const auto u = run.at("control_values").doubles();
    const auto z0 = krom.evaluateSplit(x.data());
    check(std::abs(u[0] - (K[0] * z0[0] + K[1] * z0[1])) < 1e-12, "applied input is not u = K z");

    // Output feedback: KLQE on position, same K-ROM and map.
    auto feedback = closedLoopJob(mapped, path, x0);
    feedback.at("controller").at("settings").object()["estimator"] = json::parse(
        R"({"measurement":"sensor","measurement_observable":"position","output":"state",
            "process_covariance":1,"measurement_covariance":0.001,"initial_lifted":[0,0]})");
    const auto estimated = c::readNpz(w::executeJob(feedback, root / "klqe").archives.front());
    const auto xe = estimated.at("values").doubles();
    check(std::hypot(xe[xe.size() - 2], xe[xe.size() - 1]) < 1e-3, "KLQE + certified gain did not settle");
    check(json::dump(json::parse(estimated.metadataJson())).find("not certified: the design assumes z is known") != std::string::npos,
          "output-feedback limit not reported");

    // Refused deployments.
    auto badClock = closedLoopJob(mapped, path, x0);
    badClock.at("solver").at("control_sample_time") = 0.01;
    rejects([&] { (void)w::executeJob(badClock, root / "bad-clock"); }, "wrong control clock");
    auto badCap = closedLoopJob(mapped, path, x0);
    badCap.at("controller").at("settings").object()["control_limit"] = 0.5;
    rejects([&] { (void)w::executeJob(badCap, root / "bad-cap"); }, "actuator cap inside the certified box");
    auto wrongKrom = closedLoopJob(kromPath, path, x0);
    rejects([&] { (void)w::executeJob(wrongKrom, root / "bad-krom"); }, "K-ROM without the certified map");

    // A target that is not an equilibrium yields an honest archive without a gain.
    const auto off = c::readNpz(certify({{"target", json::Value::Object{{"state", json::Value::Array{0.3, 0.0}}}}}));
    const auto offMeta = json::parse(off.metadataJson());
    check(offMeta.at("bilinear_certificate").at("status").string() == "target_not_equilibrium", "constant term not reported");
    check(!off.contains("certificate_K"), "uncertified archive carries a gain");
    auto uncertified = closedLoopJob(mapped, certify({{"target", json::Value::Object{{"state", json::Value::Array{0.3, 0.0}}}}}), x0);
    rejects([&] { (void)w::executeJob(uncertified, root / "uncertified"); }, "uncertified archive deployed");

    // Fitted operators keep A too.
    const auto fitted = c::readNpz(certify({{"operators", "fit"}}));
    check(fitted.at("model_A").doubles() == exact, "operator fit altered A");
}
}  // namespace

int main(int argc, char** argv) {
    const auto root = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path()
        / ("dynsyslab-bilinear-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        fs::create_directories(root);
        frozenFit();
        supportedSectorRules();
        bridgeConversion();
        inputMapKeepsA(root);
        cleanRecovery(root);
        historyHeldInput(root);
        estimatorAndControllerUseBilinear();
        endToEnd(root);
        std::cout << "bilinear certificate tests passed: " << checks << " checks\n";
        if (argc <= 1) { std::error_code ignored; fs::remove_all(root, ignored); }
        return 0;
    } catch (const std::exception& failure) {
        std::cerr << failure.what() << '\n';
        return 1;
    }
}
