// Path representation error of an approximate input map and the rank/condition audit of a
// controlled training design (ControlMapAudit.hpp). Ported from
// scripts/generate_rpe_control_identification_data.py (_representation_error, _dataset_diagnostics).
// Percentiles interpolate linearly between order statistics (Hyndman and Fan 1996, definition 7).
#include <dynsyslab/core/ControlMapAudit.hpp>

#include <dynsyslab/numerics/Svd.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

double percentile(std::vector<double> v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double position = q / 100.0 * double(v.size() - 1);
    const auto low = std::size_t(std::floor(position));
    const auto high = std::min(low + 1, v.size() - 1);
    return v[low] + (position - double(low)) * (v[high] - v[low]);
}

// Rank by the numpy default tolerance: s_max * max(rows, columns) * epsilon.
std::size_t rank(const std::vector<double>& s, std::size_t rows, std::size_t columns) {
    if (s.empty()) return 0;
    const double tolerance = s.front() * double(std::max(rows, columns)) * std::numeric_limits<double>::epsilon();
    return std::size_t(std::count_if(s.begin(), s.end(), [&](double v) { return v > tolerance; }));
}

}  // namespace

ControlMapPathError controlMapPathError(const KoopmanModel& model, const std::string& approximation, const std::vector<double>& states,
                                        std::size_t samples, double relativeFloor) {
    const auto n = model.stateDimension();
    if (!samples || states.size() != samples * n) throw std::invalid_argument("control-map path error: samples x state dimension states");
    ControlMapPathError out;
    out.samples = samples;
    std::vector<double> relative(samples);
    double error = 0, energy = 0;
    for (std::size_t k = 0; k < samples; ++k) {
        const double* x = states.data() + k * n;
        const auto exact = model.controlMap(x, "exact");
        const auto fitted = model.controlMap(x, approximation);
        if (exact.size() != fitted.size()) throw std::invalid_argument("control-map path error: map shapes differ");
        double e = 0, m = 0;
        for (std::size_t i = 0; i < exact.size(); ++i) { e += std::norm(fitted[i] - exact[i]); m += std::norm(exact[i]); }
        error += e; energy += m;
        relative[k] = std::sqrt(e) / std::max({std::sqrt(m), relativeFloor, std::numeric_limits<double>::min()});
    }
    out.rmse = std::sqrt(error / double(samples));
    out.relativeFrobenius = energy > 0 ? std::sqrt(error / energy) : std::numeric_limits<double>::infinity();
    out.medianRelative = percentile(relative, 50);
    out.p95Relative = percentile(relative, 95);
    out.maximumRelative = *std::max_element(relative.begin(), relative.end());
    return out;
}

ControlDesignAudit controlDesignAudit(const KoopmanModel& model, const std::vector<double>& states, const std::vector<double>& inputs,
                                      std::size_t samples, std::size_t m, std::size_t stride) {
    const auto n = model.stateDimension();
    if (!samples || !m || !stride || states.size() != samples * n || inputs.size() != samples * m)
        throw std::invalid_argument("control design audit: aligned states and held inputs");
    ControlDesignAudit out;
    out.inputs = m;
    out.inputMinimum = *std::min_element(inputs.begin(), inputs.end());
    out.inputMaximum = *std::max_element(inputs.begin(), inputs.end());
    double sum = 0, square = 0;
    for (const auto u : inputs) { sum += u; square += u * u; }
    out.inputMean = sum / double(inputs.size());
    out.inputRms = std::sqrt(square / double(inputs.size()));
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> column(samples);
        for (std::size_t k = 0; k < samples; ++k) column[k] = states[k * n + j];
        out.stateMinimum.push_back(*std::min_element(column.begin(), column.end()));
        out.stateMaximum.push_back(*std::max_element(column.begin(), column.end()));
        out.stateP01.push_back(percentile(column, 1));
        out.stateP99.push_back(percentile(column, 99));
    }
    std::vector<double> constant;
    std::vector<std::complex<double>> bilinear;
    std::size_t modes = 0;
    for (std::size_t k = 0; k < samples; k += stride) {
        const auto psi = model.evaluateFromState(states.data() + k * n);
        modes = psi.size();
        for (std::size_t j = 0; j < m; ++j) {
            const double u = inputs[k * m + j];
            constant.push_back(u);
            for (const auto& p : psi) bilinear.push_back(u * p);
        }
        for (std::size_t j = 0; j < m; ++j) bilinear.push_back(inputs[k * m + j]);
        ++out.rows;
    }
    // Rows were written input-major per sample: [u_1 psi, ..., u_m psi, u_1..u_m].
    out.bilinearColumns = m * modes + m;
    out.constantRank = rank(numerics::singularValues(std::vector<std::complex<double>>(constant.begin(), constant.end()), out.rows, m), out.rows, m);
    const auto s = numerics::singularValues(bilinear, out.rows, out.bilinearColumns);
    out.bilinearRank = rank(s, out.rows, out.bilinearColumns);
    out.bilinearCondition = s.empty() || !(s.back() > 0) ? std::numeric_limits<double>::infinity() : s.front() / s.back();
    return out;
}

}  // namespace dynsyslab::core
