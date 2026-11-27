#include <dynsyslab/core/SpectrumComparison.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace dynsyslab::core {

double symmetricNearestRms(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b) {
    if (a.empty() || b.empty()) return std::numeric_limits<double>::quiet_NaN();
    const auto directed = [](const std::vector<std::complex<double>>& from, const std::vector<std::complex<double>>& to) {
        double sum = 0.0;
        for (const auto& x : from) {
            double best = std::numeric_limits<double>::infinity();
            for (const auto& y : to) best = std::min(best, std::norm(x - y));
            sum += best;
        }
        return sum / static_cast<double>(from.size());
    };
    return std::sqrt(0.5 * (directed(a, b) + directed(b, a)));
}

}  // namespace dynsyslab::core
