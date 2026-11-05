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
};
} // namespace dynsyslab::core
