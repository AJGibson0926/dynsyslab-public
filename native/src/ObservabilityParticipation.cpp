// Mode participation of observability directions (core/Observability.hpp). The participation
// of a coordinate group in a direction is the squared norm of the direction's components on that
// group, as in modal participation factors (Perez-Arriaga, Verghese & Schweppe, IEEE Trans. Power
// Apparatus and Systems PAS-101(9), 1982).
#include <dynsyslab/core/Observability.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <stdexcept>
#include <vector>

namespace dynsyslab::core {

std::vector<double> observabilityModeParticipation(const ObservabilityResult& result,
                                                   const std::vector<std::vector<std::size_t>>& groups) {
    const std::size_t n = result.stateDimension;
    if (result.rightDirections.size() != n * n)
        throw std::invalid_argument("observability result has no square right-direction matrix");
    std::vector<bool> used(n, false);
    std::vector<double> participation(groups.size() * n, 0.0);
    for (std::size_t g = 0; g < groups.size(); ++g) {
        if (groups[g].empty()) throw std::invalid_argument("mode participation groups may not be empty");
        for (const auto index : groups[g]) {
            if (index >= n) throw std::invalid_argument("mode participation group index is outside the state");
            if (used[index]) throw std::invalid_argument("mode participation groups must not overlap");
            used[index] = true;
            for (std::size_t j = 0; j < n; ++j) {
                const double value = result.rightDirections[index * n + j];
                participation[g * n + j] += value * value;
            }
        }
    }
    return participation;
}

std::vector<std::vector<std::size_t>> splitModeGroups(const KoopmanModel& model) {
    // generatorSplit has a 2x2 rotation block exactly where a complex mode owns two coordinates.
    const auto generator = model.generatorSplit();
    const std::size_t n = model.splitDimension();
    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t row = 0; row < n;) {
        if (row + 1 < n && generator[row * n + row + 1] != 0.0) {
            groups.push_back({row, row + 1});
            row += 2;
        } else {
            groups.push_back({row});
            ++row;
        }
    }
    return groups;
}

}  // namespace dynsyslab::core
