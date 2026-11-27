#include "ObjectPublishing.hpp"
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/LongPath.hpp>

#include <chrono>
#include <optional>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

struct Stage {
    std::string method;
    J settings;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
};
thread_local std::vector<Stage> stages;
thread_local std::vector<std::filesystem::path> additional;

}  // namespace

StageScope::StageScope(std::string method, J settings) { stages.push_back({std::move(method), std::move(settings)}); }
StageScope::~StageScope() { stages.pop_back(); }

std::string currentStageMethod() { return stages.empty() ? std::string() : stages.back().method; }

std::filesystem::path kindFolder(const std::filesystem::path& root, const std::string& kind) {
    return c::objectDirectory(root, kind);
}

std::string withDataObject(const std::string& metadataJson, const std::string& kind,
                           const std::vector<std::filesystem::path>& inputs) {
    c::NpzRecord carrier;
    carrier.arrays["metadata_json"] = c::NpyArray::fromString(metadataJson.empty() ? "{}" : metadataJson);
    stampDataObject(carrier, kind, kind, inputs);
    return carrier.metadataJson();
}

void stampDataObject(c::NpzRecord& archive, const std::string& kind, const std::string& name,
                     const std::vector<std::filesystem::path>& inputs, J properties) {
    auto meta = archive.contains("metadata_json") && !archive.metadataJson().empty()
        ? json::parse(archive.metadataJson()) : J(J::Object{});
    if (!meta.at("data_object").isNull()) return;
    c::DataObjectInfo info;
    info.kind = kind;
    info.producedBy.method = stages.empty() ? name : stages.back().method;
    info.producedBy.settings = stages.empty() ? J() : stages.back().settings;
    auto sources = inputs;  // none given: the archives its metadata names as its sources
    if (sources.empty())
        for (const char* key : {"source_archive", "source_krom", "source_spectrum"})
            if (!meta.at(key).string().empty()) sources.push_back(std::filesystem::u8path(meta.at(key).string()));
    for (const auto& input : sources) {
        if (input.empty() || !std::filesystem::is_regular_file(c::extendedLengthPath(input))) continue;
        const auto absolute = std::filesystem::absolute(input).lexically_normal();
        info.producedBy.inputs.push_back({absolute.generic_u8string(), c::sha256File(absolute), "input"});
    }
    info.properties = std::move(properties);
    meta.object()["data_object"] = c::dataObjectJson(info);
    // Stage timing: seconds from the start of this step (source loading and the algorithm)
    // to the publication of this archive.
    if (!stages.empty() && meta.at("timing").isNull())
        meta.object()["timing"] = J(J::Object{{"schema", J("dynsyslab.stage-timing/v1")},
            {"method", J(stages.back().method)},
            {"stage_seconds", J(std::chrono::duration<double>(
                std::chrono::steady_clock::now() - stages.back().started).count())}});
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
}

c::ModelBindingOptions composedBinding(const J& settings) {
    c::ModelBindingOptions binding;
    // Which of the eigenfunctions' representations the model evaluates (by name; empty: the krom's first).
    binding.representation = settings.at("eigenfunction_representation").string();
    if (!settings.at("reconstructions").isNull())
        throw std::invalid_argument("\"reconstructions\" objects no longer exist: readouts C belong to the krom "
                                    "(fit them with observable_reconstruction and use the krom it publishes)");
    return binding;
}

void registerAdditionalObject(const std::filesystem::path& path) { additional.push_back(path); }

std::vector<std::filesystem::path> takeAdditionalObjects() {
    auto out = std::move(additional);
    additional.clear();
    return out;
}

std::vector<std::filesystem::path> collectStepObjects(const std::filesystem::path& primary, const std::string& method) {
    std::vector<std::filesystem::path> out{primary};
    for (auto& path : takeAdditionalObjects()) out.push_back(std::move(path));
    for (const auto& path : out) {
        if (path.extension() != L".npz" || !std::filesystem::is_regular_file(c::extendedLengthPath(path))) continue;  // e.g. a campaign folder
        if (c::dataObjectKind(c::readNpzIndexStored(path).metadataJson).empty())
            throw std::logic_error(method + " wrote " + path.filename().u8string() + " without declaring its object kind");
    }
    return out;
}

}  // namespace dynsyslab::workflow::postprocess_detail
