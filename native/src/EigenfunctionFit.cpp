#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <stdexcept>

namespace dynsyslab::core {
EigenfunctionFitDiagnostics
EigenfunctionFitDiagnostics::extract(const std::vector<std::size_t> &modes) const {
    auto result = *this;
    result.residualMaximum.clear();
    for (auto j : modes) {
        if (!residualMaximum.empty())
            result.residualMaximum.push_back(residualMaximum.at(j));
    }
    // A subset is not another fitting run; retain the parent fit count.
    return result;
}
} // namespace dynsyslab::core
