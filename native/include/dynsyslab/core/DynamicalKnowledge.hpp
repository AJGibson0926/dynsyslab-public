#pragma once
#include <functional>
#include <string>
#include <vector>
namespace dynsyslab::core {
// Bibliographic information is advisory; adding a paper never installs executable guides.
struct ResearchPaper {
    std::string id, title, venue, doi;
    int year{};
    std::vector<std::string> authors, urls, topics, relatedFactIds;
    std::string findings, applicability, coordinates, verification, notes;
};
struct DynamicalFact {
    std::string id, title, kind, statement, formula, conditions, coordinates;
    std::vector<std::string> references;
    bool applicable{};
    std::string applicabilityReason;
    std::string evidence{"researched relation; verify against configured equations"};
};
// A hand-derived eigenfunction of the configured system, kept as text so that systems carry their
// own eigenfunctions without depending on the Koopman libraries. Expressions use the state names and
// the parameter names; the eigenvalue is the continuous-time rate lambda with psi(x(t)) =
// exp(lambda t) psi(x(0)). kind "exact" holds on the stated region; "principal" is the linearization
// w^T (x - x*) of a principal eigenfunction at an equilibrium (accurate near x* only); "shell" is
// exact on the stated region but its rate text uses the state through an invariant (e.g.
// i n Omega(H)), so it is an eigenfunction along each orbit, with the eigenvalue of that orbit.
// Expressions may use atan2 and the elliptic integrals (AnalyticExpressionText.cpp).
// Loaded by libraryEigenfunctions (KnownEigenfunctions.hpp) as Eigenfunction objects.
struct EigenfunctionEntry {
    std::string id, title, kind{"exact"};
    std::string real, imag{"0"};               // psi(x) = real + i imag
    std::string rate, rateImag{"0"};           // lambda = rate + i rateImag
    std::string conditions, derivation;
    std::vector<std::string> references;
    bool applicable{true};                     // false when the configuration violates `conditions`
};
struct EquilibriumGuide {
    std::string factId;
    std::vector<double> state;
};
struct InvariantGuide {
    std::string factId, label;
    std::function<double(const std::vector<double> &)> value;
};
struct DynamicalKnowledge {
    std::string schema{"dynsyslab.dynamical-knowledge/v1"}, version{"1"},
        researchStatus{"not_reviewed"};
    std::vector<DynamicalFact> facts;
    std::vector<ResearchPaper> papers;
    // Suggestions only. The sidecar always checks f(x)=0 and classifies the actual Jacobian.
    std::function<std::vector<EquilibriumGuide>(const std::vector<double> &lower,
                                                const std::vector<double> &upper)>
        equilibriumGuides;
    std::vector<InvariantGuide> invariantGuides;
    std::vector<EigenfunctionEntry> eigenfunctions;
};
} // namespace dynsyslab::core
