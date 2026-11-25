// ARX identification on a held-input clock (Ljung, System Identification, 2nd ed., 1999) and its
// non-minimal companion realization, whose state is the measured delay window (Korda and Mezic,
// Automatica 93, 2018, use the same delay-embedded output idea for linear predictors).
#include <dynsyslab/core/IoDelayModel.hpp>

#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/StreamingQr.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

void requireFinite(const std::vector<double>& values, const char* name) {
    if (std::any_of(values.begin(), values.end(), [](double v) { return !std::isfinite(v); }))
        throw std::invalid_argument(std::string(name) + " must be finite");
}

std::size_t intervals(const IoDelayRecord& record) {
    return record.y.empty() ? 0 : std::min(record.y.size() - 1, record.u.size());
}

bool heldWindow(const IoDelayRecord& record, std::size_t k, std::size_t nu) {
    if (record.heldValid.empty()) return true;
    for (std::size_t j = 0; j < nu; ++j)
        if (!record.heldValid.at(k - j)) return false;
    return true;
}

double deviation(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    double mean = 0.0, sum = 0.0;
    for (double v : values) mean += v;
    mean /= static_cast<double>(values.size());
    for (double v : values) sum += (v - mean) * (v - mean);
    return std::sqrt(sum / static_cast<double>(values.size()));
}

double predict(const IoDelayModel& m, const double* y, const double* u, std::size_t k) {
    double next = m.c;
    for (std::size_t i = 0; i < m.ny; ++i) next += m.a[i] * y[k - i];
    for (std::size_t j = 0; j < m.nu; ++j) next += m.b[j] * u[k - j];
    return next;
}

}  // namespace

IoDelayRecord decimateHeldRecord(const std::vector<double>& y, const std::vector<double>& u,
                                 std::size_t decimation, double relativeTolerance) {
    if (!decimation || y.size() < 2 || u.size() + 1 < y.size())
        throw std::invalid_argument("io_delay record needs y, u aligned on one clock and a positive decimation");
    double scale = 0.0;
    for (double v : u) scale = std::max(scale, std::abs(v));
    const double tolerance = relativeTolerance * std::max(scale, 1e-300);
    const std::size_t count = (y.size() - 1) / decimation + 1;
    IoDelayRecord out;
    out.y.reserve(count); out.u.reserve(count - 1); out.heldValid.reserve(count - 1);
    for (std::size_t m = 0; m < count; ++m) {
        out.y.push_back(y[m * decimation]);
        if (m + 1 == count) break;
        const double held = u[m * decimation];
        bool valid = true;
        for (std::size_t s = m * decimation + 1; s < (m + 1) * decimation; ++s)
            valid = valid && std::abs(u[s] - held) <= tolerance;
        out.u.push_back(held);
        out.heldValid.push_back(valid ? 1 : 0);
    }
    return out;
}

std::pair<IoDelayModel, IoDelayFitReport> fitIoDelayModel(
    const std::vector<IoDelayRecord>& records, const IoDelayFitOptions& options, double sampleTime) {
    const auto ny = options.ny, nu = options.nu;
    if (!ny || !nu || !(sampleTime > 0) || !(options.relativeRidge >= 0) || records.empty())
        throw std::invalid_argument("io_delay fit needs ny, nu, a sample time, a nonnegative ridge and records");
    std::vector<double> allY, allU;
    for (const auto& record : records) {
        requireFinite(record.y, "io_delay y"); requireFinite(record.u, "io_delay u");
        allY.insert(allY.end(), record.y.begin(), record.y.end());
        allU.insert(allU.end(), record.u.begin(), record.u.end());
    }
    const double scaleY = deviation(allY), scaleU = deviation(allU);
    if (!(scaleY > 0) || !(scaleU > 0))
        throw std::invalid_argument("io_delay fit needs a varying output and input");
    const std::size_t parameters = ny + nu + (options.fitOffset ? 1 : 0);
    const std::size_t columns = parameters + 1;  // the regression target rides along in the QR
    numerics::WeightedQrAccumulator qr(columns);
    IoDelayFitReport report;
    constexpr std::size_t blockRows = 2048;
    std::vector<double> block;
    block.reserve(blockRows * columns);
    const auto flush = [&] {
        if (block.empty()) return;
        qr.append(block, block.size() / columns);
        block.clear();
    };
    const std::size_t first = std::max(ny, nu) - 1;
    for (const auto& record : records) {
        const auto usable = intervals(record);
        for (std::size_t k = first; k < usable; ++k) {
            if (!heldWindow(record, k, nu)) { ++report.rowsRejectedForHold; continue; }
            for (std::size_t i = 0; i < ny; ++i) block.push_back(record.y[k - i] / scaleY);
            for (std::size_t j = 0; j < nu; ++j) block.push_back(record.u[k - j] / scaleU);
            if (options.fitOffset) block.push_back(1.0);
            block.push_back(record.y[k + 1] / scaleY);
            ++report.rows;
            if (block.size() >= blockRows * columns) flush();
        }
    }
    flush();
    if (report.rows < parameters)
        throw std::invalid_argument("io_delay fit has fewer valid rows than parameters");
    // Solve the ridge problem on R: [R1; sqrt(lambda rows) I] theta = [r; 0].
    const auto& factor = qr.factor();
    const std::size_t factorRows = qr.rows();
    const std::size_t penalized = ny + nu;
    const std::size_t stacked = factorRows + penalized;
    std::vector<double> matrix(stacked * parameters, 0.0), rhs(stacked, 0.0);
    for (std::size_t row = 0; row < factorRows; ++row) {
        std::copy_n(factor.data() + row * columns, parameters, matrix.data() + row * parameters);
        rhs[row] = factor[row * columns + parameters];
    }
    const double penalty = std::sqrt(options.relativeRidge * static_cast<double>(report.rows));
    for (std::size_t p = 0; p < penalized; ++p)
        matrix[(factorRows + p) * parameters + p] = penalty;
    const auto theta = numerics::solveLeastSquares(matrix, stacked, parameters, rhs, 1);
    IoDelayModel model;
    model.ny = ny; model.nu = nu; model.sampleTime = sampleTime;
    model.a.assign(theta.begin(), theta.begin() + static_cast<std::ptrdiff_t>(ny));
    model.b.resize(nu);
    for (std::size_t j = 0; j < nu; ++j) model.b[j] = theta[ny + j] * scaleY / scaleU;
    model.c = options.fitOffset ? theta[ny + nu] * scaleY : 0.0;
    double squares = 0.0;
    for (const auto& record : records) {
        const auto usable = intervals(record);
        for (std::size_t k = first; k < usable; ++k) {
            if (!heldWindow(record, k, nu)) continue;
            const double e = record.y[k + 1] - predict(model, record.y.data(), record.u.data(), k);
            squares += e * e;
        }
    }
    report.oneStepRmse = std::sqrt(squares / static_cast<double>(report.rows));
    return {std::move(model), report};
}

IoDelayRealization ioDelayRealization(const IoDelayModel& m) {
    if (!m.ny || !m.nu || m.a.size() != m.ny || m.b.size() != m.nu || !(m.sampleTime > 0))
        throw std::invalid_argument("io_delay model is malformed");
    const std::size_t n = m.stateDimension();
    std::vector<double> A(n * n, 0.0), B(n, 0.0), C(n, 0.0), d(n, 0.0);
    std::copy(m.a.begin(), m.a.end(), A.begin());
    for (std::size_t j = 1; j < m.nu; ++j) A[m.ny + j - 1] = m.b[j];
    for (std::size_t i = 1; i < m.ny; ++i) A[i * n + i - 1] = 1.0;       // y delay chain
    for (std::size_t j = 1; j + 1 < m.nu; ++j)
        A[(m.ny + j) * n + m.ny + j - 1] = 1.0;                          // u delay chain
    B[0] = m.b[0];
    if (m.nu > 1) B[m.ny] = 1.0;
    C[0] = 1.0;
    d[0] = m.c;
    return {makeStateSpace(std::move(A), std::move(B), std::move(C), {0.0}, n, 1, 1,
                           TimeDomain::discrete, m.sampleTime), std::move(d)};
}

IoDelayValidation validateIoDelayModel(const IoDelayModel& m, const IoDelayRecord& record) {
    IoDelayValidation out;
    const std::size_t first = std::max(m.ny, m.nu) - 1, usable = intervals(record);
    if (usable <= first) throw std::invalid_argument("io_delay validation record is shorter than the model window");
    double squares = 0.0;
    for (std::size_t k = first; k < usable; ++k) {
        if (!heldWindow(record, k, m.nu)) continue;
        const double e = record.y[k + 1] - predict(m, record.y.data(), record.u.data(), k);
        squares += e * e; ++out.rows;
    }
    if (!out.rows) throw std::invalid_argument("io_delay validation record has no valid rows");
    const double scale = deviation(record.y);
    out.oneStepRmse = std::sqrt(squares / static_cast<double>(out.rows));
    out.oneStepRelative = scale > 0 ? out.oneStepRmse / scale : std::numeric_limits<double>::quiet_NaN();
    // Free run: only the first window comes from the record; afterwards y is the model's own.
    auto simulated = record.y;
    double freeSquares = 0.0;
    const double bound = 1e6 * std::max(scale, 1e-12);
    for (std::size_t k = first; k < usable; ++k) {
        const double next = predict(m, simulated.data(), record.u.data(), k);
        if (!std::isfinite(next) || std::abs(next) > bound + std::abs(record.y[k + 1])) {
            out.freeRunDiverged = true; break;
        }
        simulated[k + 1] = next;
        const double e = next - record.y[k + 1];
        freeSquares += e * e;
    }
    const auto steps = static_cast<double>(usable - first);
    out.freeRunRmse = out.freeRunDiverged ? std::numeric_limits<double>::infinity() : std::sqrt(freeSquares / steps);
    out.freeRunRelative = out.freeRunDiverged || !(scale > 0)
        ? std::numeric_limits<double>::infinity() : out.freeRunRmse / scale;
    return out;
}

void storeIoDelayModel(const IoDelayModel& m, NpzRecord& archive) {
    const auto realization = ioDelayRealization(m);
    const auto n = m.stateDimension();
    archive.arrays["io_delay_a"] = NpyArray::fromDoubles({m.ny}, m.a);
    archive.arrays["io_delay_b"] = NpyArray::fromDoubles({m.nu}, m.b);
    archive.arrays["io_delay_c"] = NpyArray::fromDoubles({1}, std::vector<double>{m.c});
    archive.arrays["io_delay_sample_time"] = NpyArray::fromDoubles({1}, std::vector<double>{m.sampleTime});
    archive.arrays["companion_A"] = NpyArray::fromDoubles({n, n}, realization.system.A);
    archive.arrays["companion_B"] = NpyArray::fromDoubles({n, 1}, realization.system.B);
    archive.arrays["companion_C"] = NpyArray::fromDoubles({1, n}, realization.system.C);
    archive.arrays["companion_d"] = NpyArray::fromDoubles({n}, realization.offset);
}

IoDelayModel loadIoDelayModel(const NpzRecord& archive) {
    IoDelayModel m;
    m.a = archive.at("io_delay_a").doubles();
    m.b = archive.at("io_delay_b").doubles();
    m.c = archive.at("io_delay_c").doubles().at(0);
    m.sampleTime = archive.at("io_delay_sample_time").doubles().at(0);
    m.ny = m.a.size(); m.nu = m.b.size();
    (void)ioDelayRealization(m);  // validates shape and finiteness
    return m;
}

}  // namespace dynsyslab::core
