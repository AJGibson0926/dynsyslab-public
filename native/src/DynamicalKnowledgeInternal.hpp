#pragma once
#include <dynsyslab/core/System.hpp>
namespace dynsyslab::core::detail {
DynamicalKnowledge hodgkinHuxleyKnowledge(const SystemConfiguration &);
DynamicalKnowledge kuramotoKnowledge(const SystemConfiguration &);
DynamicalKnowledge springPendulumKnowledge(const SystemConfiguration &);
DynamicalKnowledge marmottantKnowledge(const SystemConfiguration &);
DynamicalKnowledge reddySzeriKnowledge(const SystemConfiguration &);
std::vector<ResearchPaper> doublePendulumResearchPapers();
DynamicalKnowledge doublePendulumKnowledge(const SystemConfiguration &configuration);
std::vector<ResearchPaper> additionalResearchPapers(const std::string &identifier);
std::vector<ResearchPaper> builtInResearchPapers(const std::string &identifier);
DynamicalKnowledge builtInDynamicalKnowledge(const SystemDescriptor &descriptor,
                                             const SystemConfiguration &configuration);
DynamicalKnowledge rpeDynamicalKnowledge(const SystemConfiguration &configuration);
} // namespace dynsyslab::core::detail
