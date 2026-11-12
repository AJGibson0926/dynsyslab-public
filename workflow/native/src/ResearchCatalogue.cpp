#include "MetadataPaths.hpp"
#include "ResearchCatalogueLock.hpp"
#include "WorkflowSupport.hpp"
#include <algorithm>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>
#include <memory>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using dynsyslab::core::ResearchPaper;
std::string identity(const ResearchPaper &paper) {
    return paper.doi.empty() ? "id:" + paper.id : "doi:" + paper.doi;
}
std::filesystem::path catalogueDirectory(const std::filesystem::path &root, const std::string &id) {
    return root / "system_metadata" / "paper_catalogues" /
           dynsyslab::core::wide(dynsyslab::core::sha256Text(id));
}
void merge(std::vector<ResearchPaper> &papers, const ResearchPaper &paper) {
    for (const auto &existing : papers) {
        if (existing.id == paper.id || identity(existing) == identity(paper)) {
            if (json::dump(researchPaperJson(existing)) != json::dump(researchPaperJson(paper)))
                throw std::runtime_error(
                    "Conflicting saved bibliography identities; all evidence was preserved");
            return;
        }
    }
    papers.push_back(paper);
}
} // namespace
dynsyslab::core::DynamicalKnowledge
loadSystemKnowledge(const dynsyslab::core::SystemDescriptor &descriptor,
                    const std::filesystem::path &archiveRoot) {
    using namespace dynsyslab::core;
    auto knowledge = descriptor.research;
    if (archiveRoot.empty())
        return knowledge;
    const auto &id = descriptor.identifier;
    const auto directory = detail::metadataIoPath(catalogueDirectory(archiveRoot, id));
    if (!std::filesystem::exists(directory))
        return knowledge;
    std::vector<std::filesystem::path> files;
    for (const auto &entry : std::filesystem::directory_iterator(directory))
        if (entry.is_regular_file() && entry.path().extension() == ".dslpaper")
            files.push_back(entry.path());
    std::sort(files.begin(), files.end());
    for (const auto &file : files) {
        const auto record = json::parse(readNpz(file).metadataJson());
        if (record.at("schema").string() != "dynsyslab.paper-addition/v1" ||
            record.at("system").string() != id)
            throw std::runtime_error("Bibliography record has a mismatched system or schema");
        const auto paper = researchPaperFromJson(record.at("paper"));
        if (file.stem().u8string() != sha256Text(identity(paper)) ||
            record.at("content_sha256").string() !=
                sha256Text(json::dump(researchPaperJson(paper))))
            throw std::runtime_error("Bibliography record failed its identity or content check");
        merge(knowledge.papers, paper);
    }
    return knowledge;
}
dynsyslab::core::DynamicalKnowledge loadSystemKnowledge(const dynsyslab::core::ODESystem &system,
                                                        const std::filesystem::path &archiveRoot) {
    auto descriptor = system.descriptor();
    descriptor.research = system.dynamicalKnowledge();
    return loadSystemKnowledge(descriptor, archiveRoot);
}
json::Value evaluateResearchCatalogueQuery(const json::Value &request,
                                           const std::filesystem::path &archiveRoot) {
    using namespace dynsyslab::core;
    if (request.at("schema").string() != "dynsyslab.research-query/v1")
        throw std::invalid_argument("Unsupported research query schema");
    const auto action = request.at("action").string("list");
    if (action != "list" && action != "add")
        throw std::invalid_argument("Research action must be list or add");
    SystemConfiguration config;
    const auto &spec = request.at("system");
    config.identifier = spec.at("name").string(spec.at("identifier").string());
    config.preset = spec.at("preset").string();
    const auto definition = resolveSystemDefinition(config);
    const auto &descriptor = definition.descriptor;
    std::unique_ptr<detail::ResearchCatalogueLock> lock;
    if (action == "add") {
        if (archiveRoot.empty())
            throw std::invalid_argument("Adding a paper requires an archive root");
        (void)researchPaperFromJson(request.at("paper")); // Validate before creating anything.
        lock = std::make_unique<detail::ResearchCatalogueLock>(
            catalogueDirectory(archiveRoot, descriptor.identifier));
    }
    auto knowledge = loadSystemKnowledge(descriptor, archiveRoot);
    std::string status = "listed", addedPath;
    if (action == "add") {
        if (archiveRoot.empty())
            throw std::invalid_argument("Adding a paper requires an archive root");
        const auto paper = researchPaperFromJson(request.at("paper"));
        const ResearchPaper *duplicate = nullptr;
        for (const auto &existing : knowledge.papers) {
            if (existing.id == paper.id && identity(existing) != identity(paper))
                throw std::invalid_argument("Paper id already belongs to a different reference");
            if (identity(existing) == identity(paper))
                duplicate = &existing;
        }
        if (duplicate)
            status = "already_present";
        else {
            const auto directory = catalogueDirectory(archiveRoot, descriptor.identifier);
            const auto file = directory / wide(sha256Text(identity(paper)) + ".dslpaper");
            NpzRecord archive;
            archive.arrays["metadata_json"] = NpyArray::fromString(json::dump(json::Value::Object{
                {"schema", "dynsyslab.paper-addition/v1"},
                {"system", descriptor.identifier},
                {"created_utc", detail::timestamp()},
                {"origin", "local addition"},
                {"content_sha256", sha256Text(json::dump(researchPaperJson(paper)))},
                {"paper", researchPaperJson(paper)}}));
            try {
                writeNpzAtomicNew(detail::metadataIoPath(file), archive);
                status = "added";
            } catch (const std::exception &) {
                if (!std::filesystem::exists(detail::metadataIoPath(file)))
                    throw;
                // A concurrent writer may have added the same DOI. Read the winner below.
                status = "already_present";
            }
            knowledge = loadSystemKnowledge(descriptor, archiveRoot);
            addedPath = file.generic_u8string();
        }
    }
    auto result = dynamicalKnowledgeJson(knowledge).at("paper_catalogue");
    result.object()["system"] = descriptor.identifier;
    result.object()["status"] = status;
    result.object()["write_policy"] = "immutable additions; duplicate DOI or id is not overwritten";
    if (!addedPath.empty())
        result.object()["added_path"] = addedPath;
    return result;
}
} // namespace dynsyslab::workflow
