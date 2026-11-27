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

#include "BilinearTestSupport.hpp"

namespace {
using namespace bilinear_test;


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
    const auto& N = bridge.continuous.B[0];
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
        auto generator = split.A;
        for (std::size_t e = 0; e < 4; ++e) generator[e] += inputs[k] * split.B[0][e];
        const auto held = dynsyslab::numerics::zeroOrderHold(generator, split.B0, 2, 1, dt);
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



// Short trajectories from a grid of lifted states: a zero-input hold, then random pulses,
// so both residual channels are sampled across the same ||z|| range.

// One exact-enough held step of x'' = 2x - 0.2x^3 - x' + u (RK4), independent of the fit.

// Held-out grid inside the certified ellipsoid and input box, never used by the fit.

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
    c::saveKoopmanModel(kromPath, bare);
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
    c::saveKoopmanModel(kromPath, bare);
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
        auto generator = split.A;
        for (std::size_t e = 0; e < 4; ++e) generator[e] += inputs[k] * split.B[0][e];
        const auto held = dynsyslab::numerics::zeroOrderHold(generator, split.B0, 2, 1, dt);
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
    c::HankelLeftBasis basis; basis.rows = 2; basis.rank = 2;
    basis.leftSingularVectors = {1, 0, 0, 1}; basis.singularValues = {1, 1};
    c::FrozenHistoryTransform transform(basis, dt, {0, 1}, "position");
    c::FunctionDomain q; q.dimension = 2;
    q.providerIdentity = transform.identity();
    c::EigenfunctionFitResult bank{c::FunctionExpansion({2, {1, 0, 0, 1}, {0, 0}, {1, 1}}, {0., 1.}, q), c::Spectrum::fromContinuousRates({lambda}, dt), {}};
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
        std::cout << "bilinear certificate tests passed: " << checks << " checks\n";
        if (argc <= 1) { std::error_code ignored; fs::remove_all(root, ignored); }
        return 0;
    } catch (const std::exception& failure) {
        std::cerr << failure.what() << '\n';
        return 1;
    }
}
