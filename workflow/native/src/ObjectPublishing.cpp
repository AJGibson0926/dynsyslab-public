#include "ObjectPublishing.hpp"
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/LongPath.hpp>

#include <optional>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

struct Stage {
    std::string method;
    J settings;
};
thread_local std::vector<Stage> stages;
thread_local std::vector<std::filesystem::path> additional;

bool startsWith(const std::string& text, const char* prefix) { return text.rfind(prefix, 0) == 0; }

}  // namespace

StageScope::StageScope(std::string method, J settings) { stages.push_back({std::move(method), std::move(settings)}); }
StageScope::~StageScope() { stages.pop_back(); }

std::string objectKindFor(const std::string& name) {
    namespace k = c::object_kind;
    if (name == "koopman-hankel-dmd") return k::eigenfunctions;
    if (name == "reconstructions") return k::reconstructions;
    if (name == "sindy-eigenfunctions" || name == "koopman-dmd" || name == "edmd"
        || name == "koopman-eigenfunction-expansion")
        return k::eigenfunctions;
    if (name == "havok") return k::dynamicalModel;
    if (name == "coordinate-transform" || name == "history-projection" || name == "hankel_svd" || name == "hankel_pod"
        || name == "klqe")
        return k::timeSeries;
    if (name == "history-basis" || name == "pca" || name == "svd" || name == "pod" || name == "bagged_pca"
        || name == "robust_pca")
        return k::coordinateMap;
    if (name == "sindy-dynamics") return k::dynamicalModel;
    if (name == "hankel_matrix") return k::matrix;
    if (name == "domain-field") return k::field;
    if (name == "control-target" || name == "reference-signal") return k::timeSeries;
    if (name == "history-transform") return k::coordinateMap;
    if (name == "io-delay-model") return k::dynamicalModel;
    if (name == "animation" || name == "figure-dataset") return "figure";
    if (name == "koopman-validation" || name == "koopman-eigenvalue-optimization" || name == "koopman-bilinear-certificate"
        || name == "koopman-mode-selection")
        return k::report;
    if (startsWith(name, "koopman-") || name == "echo-history-koopman") return k::koopmanModel;
    return k::report;  // analyses, validations, error tables
}

std::filesystem::path objectFolder(const std::filesystem::path& root, const std::string& name) {
    return c::objectDirectory(root, objectKindFor(name));
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
    for (const auto& input : inputs) {
        if (input.empty() || !std::filesystem::is_regular_file(c::extendedLengthPath(input))) continue;
        const auto absolute = std::filesystem::absolute(input).lexically_normal();
        info.producedBy.inputs.push_back({absolute.generic_u8string(), c::sha256File(absolute), "input"});
    }
    info.properties = std::move(properties);
    meta.object()["data_object"] = c::dataObjectJson(info);
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
}

c::ModelBindingOptions composedBinding(const J& settings) {
    c::ModelBindingOptions binding;
    const auto& given = settings.at("reconstructions");
    if (given.isNull()) return binding;
    if (std::holds_alternative<std::string>(given.storage())) binding.reconstructions.push_back(std::filesystem::u8path(given.string()));
    else for (const auto& item : given.array()) binding.reconstructions.push_back(std::filesystem::u8path(item.string()));
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
        if (!c::dataObjectKind(c::readNpzIndexStored(path).metadataJson).empty()) continue;
        auto archive = c::readNpzStored(path);
        const auto folder = path.parent_path();
        const auto kind = folder.parent_path().filename() == L"objects" ? folder.filename().u8string() : objectKindFor(method);
        std::vector<std::filesystem::path> inputs;
        if (!archive.metadataJson().empty()) {
            const auto meta = json::parse(archive.metadataJson());
            for (const char* key : {"source_archive", "source_krom", "source_spectrum"})
                if (!meta.at(key).string().empty()) inputs.push_back(std::filesystem::u8path(meta.at(key).string()));
        }
        stampDataObject(archive, kind, method, inputs);
        c::writeNpzAtomic(path, archive);
    }
    return out;
}

}  // namespace dynsyslab::workflow::postprocess_detail
