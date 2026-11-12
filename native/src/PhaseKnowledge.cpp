#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::phase_detail {
void applyEquilibriumGuides(const ODESystem &system, PhaseSpaceResult &result) {
    const auto knowledge = system.dynamicalKnowledge();
    if (!knowledge.equilibriumGuides)
        return;
    auto &options = result.options;
    const auto guides = knowledge.equilibriumGuides(options.lower, options.upper);
    if (guides.size() > 4096)
        throw std::invalid_argument("System supplied too many equilibrium guides");
    for (const auto &guide : guides) {
        PhaseKnowledgeCheck check;
        check.factId = guide.factId;
        check.kind = "equilibrium-guide";
        check.state = guide.state;
        check.status = "rejected";
        if (guide.state.size() != options.lower.size() || !inside(guide.state, options))
            check.reason = "Candidate is outside the declared state domain";
        else
            try {
                check.residual = norm(field(system, guide.state), options);
                check.status =
                    check.residual <= options.residualTolerance ? "verified" : "seed-only";
                check.reason = "Checked against the actual zero-input field; Jacobian "
                               "classification is computed by the equilibrium search";
                options.equilibriumSeeds.push_back(guide.state);
            } catch (const std::exception &error) {
                check.reason = error.what();
            }
        result.knowledgeChecks.push_back(std::move(check));
    }
}
void checkKnownInvariants(const ODESystem &system, PhaseSpaceResult &result) {
    const auto knowledge = system.dynamicalKnowledge();
    const auto n = result.options.lower.size();
    for (const auto &invariant : knowledge.invariantGuides) {
        if (!invariant.value)
            continue;
        for (std::size_t index = 0; index < result.manifolds.size(); ++index) {
            const auto &curve = result.manifolds[index];
            if (curve.values.size() < n)
                continue;
            PhaseKnowledgeCheck check;
            check.factId = invariant.factId;
            check.kind = "manifold-invariant";
            check.objectIndex = index;
            check.status = "unresolved";
            try {
                PhaseState x(curve.values.begin(), curve.values.begin() + n);
                const double baseline = invariant.value(x);
                if (!std::isfinite(baseline))
                    throw std::runtime_error("Nonfinite initial invariant");
                for (std::size_t row = 0; row < curve.values.size() / n; ++row) {
                    std::copy_n(curve.values.begin() + row * n, n, x.begin());
                    const double value = invariant.value(x);
                    if (!std::isfinite(value))
                        throw std::runtime_error("Nonfinite invariant along candidate");
                    check.residual = std::max(check.residual, std::abs(value - baseline) /
                                                                  (1 + std::abs(baseline)));
                }
                check.status = check.residual <= 100 * (result.options.relativeTolerance +
                                                        result.options.absoluteTolerance)
                                   ? "consistent"
                                   : "discrepancy";
                check.reason = "Sampled conservation defect relative to the first manifold sample, "
                               "scaled by 1+abs(initial invariant); diagnostic, not a proof";
            } catch (const std::exception &error) {
                check.reason = error.what();
            }
            result.knowledgeChecks.push_back(std::move(check));
        }
    }
}
} // namespace dynsyslab::core::phase_detail
