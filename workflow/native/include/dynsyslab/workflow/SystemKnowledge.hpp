#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <filesystem>
namespace dynsyslab::workflow {
namespace json = dynsyslab::core::json;
json::Value researchPaperJson(const dynsyslab::core::ResearchPaper &paper);
dynsyslab::core::ResearchPaper researchPaperFromJson(const json::Value &value);
json::Value dynamicalKnowledgeJson(const dynsyslab::core::DynamicalKnowledge &knowledge);
// Merge immutable local bibliography additions into the system-owned built-in metadata.
dynsyslab::core::DynamicalKnowledge loadSystemKnowledge(const dynsyslab::core::ODESystem &system,
                                                        const std::filesystem::path &archiveRoot);
dynsyslab::core::DynamicalKnowledge
loadSystemKnowledge(const dynsyslab::core::SystemDescriptor &descriptor,
                    const std::filesystem::path &archiveRoot);
// Versioned list/add protocol. Additions are advisory; executable guides remain in the core.
json::Value evaluateResearchCatalogueQuery(const json::Value &request,
                                           const std::filesystem::path &archiveRoot);
} // namespace dynsyslab::workflow
