#include <dynsyslab/core/TransferFunction.hpp>
#include <dynsyslab/core/SystemNorms.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
constexpr double kPi = 3.14159265358979323846;
using Cx = std::complex<double>;
}  // namespace

StateSpaceModel channel(const StateSpaceModel& m, std::size_t output, std::size_t input) {
    if (output >= m.outputs || input >= m.inputs) throw std::out_of_range("transfer-function channel is out of range");
    const std::size_t n = m.states;
    std::vector<double> B(n), C(n);
    for (std::size_t i = 0; i < n; ++i) { B[i] = m.B[i * m.inputs + input]; C[i] = m.C[output * n + i]; }
    const double D = m.D.empty() ? 0.0 : m.D[output * m.inputs + input];
    return makeStateSpace(m.A, std::move(B), std::move(C), {D}, n, 1, 1, m.domain, m.sampleTime);
}

std::vector<double> polynomialFromRoots(const std::vector<Cx>& roots) {
    std::vector<Cx> c{1.0};
    for (const auto& r : roots) {
        std::vector<Cx> next(c.size() + 1, 0.0);
        for (std::size_t k = 0; k < c.size(); ++k) { next[k] += c[k]; next[k + 1] -= r * c[k]; }
        c = std::move(next);
    }
    std::vector<double> out;
    for (const auto& v : c) out.push_back(v.real());
    return out;
}

TransferFunction transferFunction(const StateSpaceModel& model, std::size_t output, std::size_t input, double tol) {
    const auto g = channel(model, output, input);
    TransferFunction tf;
    tf.output = output; tf.input = input;
    tf.directGain = g.D[0];
    tf.poles = systemPoles(g);
    tf.zeros = g.states ? transmissionZeros(g) : std::vector<Cx>{};
    tf.relativeDegree = int(tf.poles.size()) - int(tf.zeros.size());
    // gain from G at a few points away from poles and zeros: G = k prod(s - z) / prod(s - p)
    double scale = 1.0;
    for (const auto& p : tf.poles) scale = std::max(scale, std::abs(p));
    const bool discrete = g.domain == TimeDomain::discrete;
    std::vector<double> w;
    for (double f : {0.37, 0.71, 1.13, 1.71, 2.39}) w.push_back(discrete ? f * kPi / (3.0 * g.sampleTime) : f * scale);
    const auto response = frequencyResponse(g, w).response;
    std::vector<double> gains;
    for (std::size_t i = 0; i < w.size(); ++i) {
        const Cx s = discrete ? std::exp(Cx(0.0, w[i] * g.sampleTime)) : Cx(0.0, w[i]);
        Cx k = response[i];
        for (const auto& p : tf.poles) k *= (s - p);
        for (const auto& z : tf.zeros) k /= (s - z);
        if (std::isfinite(k.real())) gains.push_back(k.real());
    }
    if (!gains.empty()) {
        std::sort(gains.begin(), gains.end());
        tf.gain = gains[gains.size() / 2];
    }
    if (std::abs(tf.gain) < 1e-300) tf.zeros.clear(), tf.relativeDegree = int(tf.poles.size());
    tf.denominator = polynomialFromRoots(tf.poles);
    tf.numerator = polynomialFromRoots(tf.zeros);
    for (auto& c : tf.numerator) c *= tf.gain;
    for (std::size_t i = 0; i < tf.poles.size(); ++i)
        for (std::size_t j = 0; j < tf.zeros.size(); ++j)
            if (std::abs(tf.poles[i] - tf.zeros[j]) <= tol * std::max(1.0, std::abs(tf.poles[i])))
                tf.cancellations.emplace_back(i, j);
    return tf;
}

std::vector<double> unwrappedPhaseDeg(const std::vector<Cx>& response) {
    std::vector<double> out;
    double previous = 0.0, offset = 0.0;
    for (std::size_t i = 0; i < response.size(); ++i) {
        const double raw = std::arg(response[i]) * 180.0 / kPi;
        if (i) {
            double step = raw + offset - previous;
            while (step > 180.0) { offset -= 360.0; step -= 360.0; }
            while (step < -180.0) { offset += 360.0; step += 360.0; }
        }
        previous = raw + offset;
        out.push_back(previous);
    }
    return out;
}

BodeData bode(const StateSpaceModel& model, std::size_t output, std::size_t input, const std::vector<double>& w) {
    const auto g = channel(model, output, input);
    const auto response = frequencyResponse(g, w).response;
    BodeData out;
    out.frequencies = w;
    for (const auto& v : response)
        out.magnitudeDb.push_back(20.0 * std::log10(std::max(std::abs(v), std::numeric_limits<double>::min())));
    out.phaseDeg = unwrappedPhaseDeg(response);
    return out;
}

}  // namespace dynsyslab::core
