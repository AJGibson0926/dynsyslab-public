// Field models (docs/architecture/FIELD_MODELS.md): every kind fitted, written as a definition,
// read back and checked against a known field; the damped-drive waveform fit, waveform weights,
// resolved directions and trajectory calibration.
#include <dynsyslab/core/FieldModels.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
void require(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

// x' = -0.1 x + y, y' = -x - 0.1 y, i.e. u' = (-0.1 - i) u for u = x + i y.
c::FieldRecord oscillator(double x0, double y0, std::size_t n, double dt, double rateScale = 1) {
    c::FieldRecord r{{}, {}, {}, 2, 2};
    for (std::size_t k = 0; k < n; ++k) {
        const double t = dt * double(k), e = std::exp(-0.1 * t);
        const double x = e * (x0 * std::cos(t) + y0 * std::sin(t)), y = e * (y0 * std::cos(t) - x0 * std::sin(t));
        r.time.push_back(t);
        r.state.insert(r.state.end(), {x, y});
        r.label.insert(r.label.end(), {rateScale * (-0.1 * x + y), rateScale * (-x - 0.1 * y)});
    }
    return r;
}

double fieldError(const c::Observable& f, const std::vector<c::FieldRecord>& records) {
    double worst = 0;
    for (const auto& r : records)
        for (std::size_t k = 0; k < r.samples(); k += 7) {
            const auto v = f.evaluateReal(r.state.data() + k * r.dimension);
            for (std::size_t j = 0; j < r.outputs; ++j) worst = std::max(worst, std::abs(v[j] - r.label[k * r.outputs + j]));
        }
    return worst;
}

c::Observable reread(const c::json::Value& d) { return c::fieldFromDefinition(c::json::parse(c::json::dump(d))); }

void testDataDrivenKinds(const std::vector<c::FieldRecord>& records) {
    c::CompactRadialOptions co;
    co.stride = 10; co.neighbors = 6; co.affineTail = true; co.ridge = 1e-8; co.threshold = 0;
    const auto compact = reread(c::fitCompactRadialField(records, co));
    require(fieldError(compact, records) < 0.02, "compact radial field with affine tail reproduces the rates");
    c::LocalTangentOptions lo;
    lo.stride = 20; lo.samples = 60; lo.neighbors = 4; lo.ridge = 1e-8; lo.threshold = 0;
    const auto local = reread(c::fitLocalTangentField(records, lo));
    require(fieldError(local, records) < 0.02, "local tangent-chart field reproduces the rates");
    c::PairedPhaseOptions po;
    po.threshold = 1e-5; po.ridge = 1e-10;
    const auto pairedDefinition = c::fitPairedPhaseField(records, po);
    require(pairedDefinition.at("harmonics").array().at(0).number() == 1, "winding number of the reference pair");
    const auto paired = reread(pairedDefinition);
    const double pairedError = fieldError(paired, records);
    require(pairedError < 1e-3, "paired amplitude/phase field reproduces u' = (-0.1 - i) u: " + std::to_string(pairedError));
    c::PairedPhaseOptions weak;
    weak.weakWindows = {5, 10}; weak.instantWeight = 0.5; weak.threshold = 1e-5; weak.ridge = 1e-10;
    require(fieldError(reread(c::fitPairedPhaseField(records, weak)), records) < 5e-3, "weak-form paired field");
    const auto roll = c::rolloutField(paired, {records[0].state[0], records[0].state[1]}, records[0].time);
    require(roll.complete && std::abs(roll.state.back() - records[0].state.back()) < 0.01, "paired rollout follows the record");
}

void testStructured(const std::vector<c::FieldRecord>& records) {
    c::MapFitOptions linear;
    linear.degree = 1; linear.threshold = 0; linear.ridge = 1e-12;
    // A second-order reading of the oscillator: (x, y) with y = x' is not kinematic, so use
    // q'' = -q - 0.2 q' sampled at random states.
    c::FieldRecord second{{}, {}, {}, 2, 2};
    for (int k = 0; k < 60; ++k) {
        const double q = std::sin(0.7 * k), v = std::cos(1.3 * k);
        second.time.push_back(k); second.state.insert(second.state.end(), {q, v}); second.label.insert(second.label.end(), {v, -q - 0.2 * v});
    }
    require(fieldError(reread(c::fitSecondOrderField({second}, linear)), {second}) < 1e-8, "second-order field: kinematic rows and acceleration");
    // Interior field on three taps: outer accelerations -x0 and -x4, interior drive f = x0 + 0.5 x5, gamma 2.
    c::FieldRecord interior{{}, {}, {}, 6, 6};
    const double g = 2;
    for (int k = 0; k < 120; ++k) {
        double z[6];
        for (int i = 0; i < 6; ++i) z[i] = std::sin(0.37 * k * (i + 1) + i);
        interior.time.push_back(k);
        interior.state.insert(interior.state.end(), z, z + 6);
        const double f = z[0] + 0.5 * z[5];
        interior.label.insert(interior.label.end(), {z[1], -z[0], z[3], f - 2 * g * z[3] - g * g * z[2], z[5], -z[4]});
    }
    c::InteriorOptions io;
    io.gamma = {g}; io.outer = linear; io.drive = linear;
    require(fieldError(reread(c::fitInteriorField({interior}, io)), {interior}) < 1e-8, "interior field: outer map and damped interior drive");
    // A two-tap pair network with symmetric weights reproduces kinematic rows and is evaluable.
    c::FieldRecord network{{}, {}, {}, 4, 4};
    for (const auto& r : records)
        for (std::size_t k = 0; k < r.samples(); ++k) {
            network.time.push_back(double(network.time.size()));
            const double* x = r.state.data() + 2 * k;
            network.state.insert(network.state.end(), {x[0], x[1], x[1], -x[0]});
            network.label.insert(network.label.end(), {x[1], -x[0], -x[0], -x[1]});
        }
    const auto net = reread(c::fitPairNetwork({network}, {0, 1}, {0, 1, 1, 0}, linear));
    require(fieldError(net, {network}) < 1e-8, "pair network shares one acceleration map");
}

void testDampedDrive() {
    // x'' + 2 g x' + g^2 x = g^2 u(t), u = sin t, by fine RK4; the drive f(u) = g^2 u is recovered.
    const double g = 3, dt = 0.01;
    c::DriveRecord r{{}, {}, {}, {}, 1, 1};
    double x = 0.2, v = -0.1;
    const auto rate = [g](double t, double p, double q) { return std::pair{q, g * g * std::sin(t) - 2 * g * q - g * g * p}; };
    for (int k = 0; k <= 1000; ++k) {
        const double t = dt * k;
        r.time.push_back(t); r.input.push_back(std::sin(t)); r.position.push_back(x); r.velocity.push_back(v);
        const int sub = 20; const double h = dt / sub;
        for (int s = 0; s < sub; ++s) {
            const double tt = t + s * h;
            const auto [a1, b1] = rate(tt, x, v);
            const auto [a2, b2] = rate(tt + h / 2, x + h / 2 * a1, v + h / 2 * b1);
            const auto [a3, b3] = rate(tt + h / 2, x + h / 2 * a2, v + h / 2 * b2);
            const auto [a4, b4] = rate(tt + h, x + h * a3, v + h * b3);
            x += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4); v += h / 6 * (b1 + 2 * b2 + 2 * b3 + b4);
        }
    }
    c::MapFitOptions o;
    o.degree = 1; o.threshold = 0; o.ridge = 1e-12;
    const auto phi = reread(c::fitDampedDrive({r}, g, o));
    const double u = 0.5;
    const double driven = phi.evaluateReal(&u)[0];
    require(std::abs(driven - g * g * 0.5) < 2e-2, "filtered-library fit recovers the drive: " + std::to_string(driven));
}

void testWaveformObjective() {
    std::vector<double> t, y;
    for (int k = 0; k <= 100; ++k) { t.push_back(0.01 * k); y.insert(y.end(), {2.0, -1.0}); }
    const auto w = c::waveformWeights(t, y, 2, {1}, {0.5}, 0.25, {1e-6});
    for (std::size_t ch = 0; ch < 2; ++ch) {
        double loss = 0;
        for (std::size_t k = 0; k < t.size(); ++k) loss += std::pow(w[k * 2 + ch] * y[k * 2 + ch], 2);
        require(std::abs(loss - 1.5) < 1e-9, "a unit relative error costs global plus local weight");
    }
    std::vector<double> responses;
    const auto dirs = c::resolvedDirections({1, 0, 0, 1e-12, 0, 0}, 3, {1, 0, 0, 1, 0, 0}, 3, 2, responses);
    require(dirs.size() == 3 && responses.size() == 3, "unresolved direction dropped");
    require(dirs[2] == 0, "zero coefficient rows stay zero");
}

void testCalibration() {
    std::vector<c::FieldRecord> biased = {oscillator(1, 0, 200, 0.05, 1.08)}, truth = {oscillator(1, 0, 200, 0.05)};
    c::MapFitOptions linear;
    linear.degree = 1; linear.threshold = 0.01; linear.ridge = 1e-12;
    const auto model = c::fitMap(biased, linear);
    c::CalibrationOptions o;
    o.duration = 4; o.startsPerRecord = 3; o.maximumEvaluations = 30;
    const auto result = c::calibrateField(model, truth, o);
    require(result.accepted && result.candidateLoss < 0.1 * result.initialLoss, "calibration reduces the shooting loss");
    require(fieldError(reread(result.definition), truth) < 0.01, "calibrated field matches the true rates");
}

}  // namespace

int main() {
    try {
        const std::vector<c::FieldRecord> records = {oscillator(1, 0, 400, 0.05), oscillator(0.3, 0.8, 400, 0.05)};
        testDataDrivenKinds(records);
        testStructured(records);
        testDampedDrive();
        testWaveformObjective();
        testCalibration();
        std::cout << "field model tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
