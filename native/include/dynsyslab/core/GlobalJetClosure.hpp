#pragma once

#include <dynsyslab/core/Observable.hpp>
#include <utility>

namespace dynsyslab::core {

// Global identities refer to the regular real domain of the supplied autonomous
// expressions, with each stored binary64 coefficient interpreted exactly.
struct GlobalJetOptions {
    std::size_t highestOrder{4};
    std::size_t polynomialDegree{2};
    std::vector<std::size_t> positiveAxes;
    bool emitSolverProblems{};
    std::size_t solverTimeoutMs{}; // per decision, zero disables external solving
    // Approximate pairs to validate by exact interval arithmetic.
    std::vector<std::pair<std::vector<double>, std::vector<double>>> pairSeeds;
    bool searchAllOrders{}; // keep checking for smoother equations or a state inverse
};

[[nodiscard]] json::Value globalJetClosure(const std::vector<AnalyticExpression>& drift,
    const Observable& coordinates, const GlobalJetOptions& options = {});

// Numerical evaluation of an explicit certified order on a realizable jet,
// using the real branch for odd roots. Implicit graphs require a fiber solve.
[[nodiscard]] std::vector<double> evaluateJetRhs(const json::Value& order,
    const std::vector<double>& jetState);

} // namespace dynsyslab::core
