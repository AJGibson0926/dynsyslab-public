// Known linear ARX system: identification, realization, record boundaries, held-input
// decimation, and LQ tracking with the feedforward from the backward recursion.
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/IoDelayModel.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>

namespace c = dynsyslab::core;

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

c::IoDelayModel truth() {
    c::IoDelayModel m;  // poles 0.9 and 0.5, stable, with one input delay
    m.ny = 2; m.nu = 2; m.a = {1.4, -0.45}; m.b = {0.2, 0.1}; m.c = 0.05; m.sampleTime = 0.1;
    return m;
}

double advance(const c::IoDelayModel& m, const std::vector<double>& y, const std::vector<double>& u, std::size_t k) {
    double next = m.c;
    for (std::size_t i = 0; i < m.ny; ++i) next += m.a[i] * (k >= i ? y[k - i] : 0.0);
    for (std::size_t j = 0; j < m.nu; ++j) next += m.b[j] * (k >= j ? u[k - j] : 0.0);
    return next;
}

c::IoDelayRecord simulate(const c::IoDelayModel& m, std::size_t steps, unsigned seed, double start) {
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal;
    c::IoDelayRecord record;
    record.y.assign(steps + 1, start);
    record.u.resize(steps);
    for (auto& value : record.u) value = normal(generator);
    for (std::size_t k = 1; k < steps; ++k) record.y[k + 1] = advance(m, record.y, record.u, k);
    return record;
}

void recoversCoefficientsAcrossRecords() {
    const auto m = truth();
    // Very different levels: a row bridging the two records would corrupt the fit.
    const std::vector<c::IoDelayRecord> records{simulate(m, 400, 1, 0.0), simulate(m, 400, 2, 25.0)};
    c::IoDelayFitOptions options;
    options.ny = 2; options.nu = 2; options.relativeRidge = 1e-14;
    const auto [fit, report] = c::fitIoDelayModel(records, options, m.sampleTime);
    for (std::size_t i = 0; i < 2; ++i) {
        require(std::abs(fit.a[i] - m.a[i]) < 1e-6, "a coefficient");
        require(std::abs(fit.b[i] - m.b[i]) < 1e-6, "b coefficient");
    }
    require(std::abs(fit.c - m.c) < 1e-6, "offset");
    require(report.rows == 2 * 399, "rows stay inside each record");
    const auto validation = c::validateIoDelayModel(fit, simulate(m, 300, 3, 0.0));
    require(validation.oneStepRelative < 1e-6 && validation.freeRunRelative < 1e-5, "held-out errors");
}

void companionMatchesRecursion() {
    const auto m = truth();
    const auto r = c::ioDelayRealization(m);
    require(r.system.states == 3, "non-minimal dimension");
    const auto record = simulate(m, 40, 4, 0.0);
    std::vector<double> z{record.y[1], record.y[0], record.u[0]};
    for (std::size_t k = 1; k < 40; ++k) {
        std::vector<double> next(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) next[i] += r.system.A[i * 3 + j] * z[j];
            next[i] += r.system.B[i] * record.u[k] + r.offset[i];
        }
        z = next;
        require(std::abs(z[0] - record.y[k + 1]) < 1e-10, "companion output");
        require(std::abs(z[2] - record.u[k]) < 1e-14, "input delay state");
    }
}

void heldInputDecimation() {
    std::vector<double> y(21, 1.0), u(20, 1.0);
    for (std::size_t k = 12; k < 20; ++k) u[k] = -1.0;  // switch inside the [10, 15) interval
    const auto record = c::decimateHeldRecord(y, u, 5);
    require(record.y.size() == 5 && record.u.size() == 4, "decimated length");
    require(record.heldValid[0] && record.heldValid[1] && !record.heldValid[2] && record.heldValid[3],
            "unheld interval flagged");
    auto fitRecord = simulate(truth(), 200, 5, 0.0);
    fitRecord.heldValid.assign(fitRecord.u.size(), 1);
    fitRecord.heldValid[50] = 0;
    c::IoDelayFitOptions options; options.ny = 2; options.nu = 2;
    const auto report = c::fitIoDelayModel({fitRecord}, options, 0.1).second;
    require(report.rowsRejectedForHold == 2, "rows touching an unheld interval are rejected");
}

void tracksKnownReference() {
    const auto m = truth();
    const std::size_t steps = 1500;
    std::vector<double> reference(steps + 80);
    for (std::size_t k = 0; k < reference.size(); ++k)
        reference[k] = 0.7 + 0.4 * std::sin(0.08 * static_cast<double>(k));
    auto tracker = std::make_shared<c::IoDelayTracker>(m, 1.0, 1e-6, reference);
    require(tracker->riccatiResidual() < 1e-8, "Riccati residual");
    for (const auto& pole : tracker->closedLoopPoles()) require(std::abs(pole) < 1.0, "closed loop stable");
    c::IoDelayTrackerStream stream(tracker, 0.0);
    std::vector<double> y(steps + 1, 0.0), u(steps, 0.0);
    double squares = 0.0;
    std::size_t counted = 0;
    for (std::size_t k = 0; k < steps; ++k) {
        u[k] = stream.step(k, y[k]);
        require(stream.step(k, y[k]) == u[k], "a repeated tick holds the command");
        y[k + 1] = advance(m, y, u, k);
        if (k > 200) { squares += std::pow(y[k + 1] - reference[k + 1], 2); ++counted; }
    }
    require(std::sqrt(squares / static_cast<double>(counted)) < 1e-3, "tracking error with preview");

    c::IoDelayTrackerStream capped(tracker, 0.0, std::make_pair(-0.05, 0.05));
    std::vector<double> yc(101, 0.0);
    for (std::size_t k = 0; k < 100; ++k) {
        const double value = capped.step(k, yc[k]);
        require(std::abs(value) <= 0.05 + 1e-15, "control limit");
        yc[k + 1] = 0.9 * yc[k] + value;
    }
}

}  // namespace

int main() {
    try {
        recoversCoefficientsAcrossRecords();
        companionMatchesRecursion();
        heldInputDecimation();
        tracksKnownReference();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "io_delay_model tests passed\n";
    return 0;
}
