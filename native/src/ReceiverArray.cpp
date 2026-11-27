// Synchronized receiver arrays and the linearized delay inverse (ReceiverArray.hpp). Ported from
// the receiver-reconstruction side-car (signals.cpp).
#include <dynsyslab/core/ReceiverArray.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

ReceiverArraySignals rpeReceiverArray(const ODESystem& system, const std::vector<AcousticReceiver>& receivers,
                                      const std::vector<double>& t, const std::vector<double>& x, const ReceiverArrayOptions& options) {
    const auto n = system.descriptor().stateDimension;
    if (receivers.empty()) throw std::invalid_argument("receiver array: no receivers");
    if (t.size() < 2 || x.size() != t.size() * n) throw std::invalid_argument("receiver array: trajectory shape");
    if (!options.stride) throw std::invalid_argument("receiver array: stride must be positive");
    const double T = system.coordinates().secondsPerModelTime;
    const double R0 = system.configuration().parameters.at("equilibrium_radius_m");
    ReceiverArraySignals out;
    out.receivers = receivers.size();
    out.stateDimension = n;
    std::vector<ObservableCatalog> catalogs;
    std::vector<double> delay;
    for (const auto& r : receivers) {
        catalogs.push_back(rpeReceiverPressures(system, r));
        out.ranges.push_back(r.rangeMetres());
        delay.push_back(r.travelTimeSeconds() / T);
    }
    out.nearest = std::size_t(std::min_element(out.ranges.begin(), out.ranges.end()) - out.ranges.begin());
    for (std::size_t i = 0; i < t.size(); ++i)
        if (!(x[i * n] > 0 && x[i * n] * R0 < out.ranges[out.nearest]))
            throw std::invalid_argument("receiver array: a receiver lies inside the bubble");
    for (const auto d : delay) out.lags.push_back(d - delay[out.nearest]);
    const double maxLag = *std::max_element(out.lags.begin(), out.lags.end());
    out.terms = {"monopole_pressure", "spherical_near_field_pressure", "dipole_induction_pressure", "dipole_radiation_pressure"};
    std::vector<std::vector<double>> pressure(out.terms.size()), nearAt(out.terms.size());
    std::vector<std::size_t> left(receivers.size(), 0);
    const auto value = [&](const Observable& term, std::size_t a, std::size_t b, double f) {
        const double ya = term.evaluateReal(x.data() + a * n).front();
        return f == 0 ? ya : ya + f * (term.evaluateReal(x.data() + b * n).front() - ya);
    };
    for (std::size_t i = 0; i < t.size(); i += options.stride) {
        if (t[i] > options.end) break;
        if (t[i] < options.begin || t[i] - maxLag < t.front()) continue;
        out.time.push_back(t[i]);
        out.states.insert(out.states.end(), x.begin() + std::ptrdiff_t(i * n), x.begin() + std::ptrdiff_t((i + 1) * n));
        for (std::size_t j = 0; j < receivers.size(); ++j) {
            const double emitted = t[i] - out.lags[j];
            while (left[j] + 1 < t.size() && t[left[j] + 1] <= emitted) ++left[j];
            const auto right = std::min(left[j] + 1, t.size() - 1);
            const double f = right == left[j] ? 0 : (emitted - t[left[j]]) / (t[right] - t[left[j]]);
            for (std::size_t k = 0; k < out.terms.size(); ++k) {
                const double p = value(catalogs[j].at(out.terms[k]), left[j], right, f);
                pressure[k].push_back(p);
                nearAt[k].push_back(j == out.nearest ? p : value(catalogs[out.nearest].at(out.terms[k]), left[j], right, f));
            }
        }
    }
    if (out.time.empty()) throw std::invalid_argument("receiver array: no samples in the window after propagation warm-up");
    for (std::size_t k = 0; k < out.terms.size(); ++k) {
        out.pressures.insert(out.pressures.end(), pressure[k].begin(), pressure[k].end());
        out.nearAtEmission.insert(out.nearAtEmission.end(), nearAt[k].begin(), nearAt[k].end());
    }
    return out;
}

std::vector<double> linearDelayInverse(const std::vector<double>& A, const std::vector<double>& H, const std::vector<double>& delays,
                                       std::size_t n) {
    const auto r = delays.size();
    if (A.size() != n * n || H.size() != r * n || !r) throw std::invalid_argument("linear delay inverse: shapes");
    std::vector<double> O(r * n, 0);
    for (std::size_t i = 0; i < r; ++i) {
        const auto E = numerics::matrixExponential(A, n, -delays[i]);  // exp(-A delta_i)
        for (std::size_t c = 0; c < n; ++c)
            for (std::size_t k = 0; k < n; ++k) O[i * n + c] += H[i * n + k] * E[k * n + c];
    }
    const auto svd = numerics::singularValueDecomposition(O, r, n);
    const double top = svd.singularValues.empty() ? 0 : svd.singularValues.front();
    std::size_t rank = 0;
    for (const auto s : svd.singularValues) rank += s > 1e-10 * top;
    if (rank < n) throw std::invalid_argument("linear delay inverse: the delayed outputs do not observe every state");
    std::vector<double> identity(r * r, 0);
    for (std::size_t i = 0; i < r; ++i) identity[i * r + i] = 1;
    return numerics::solveLeastSquares(O, r, n, identity, r);
}

}  // namespace dynsyslab::core
