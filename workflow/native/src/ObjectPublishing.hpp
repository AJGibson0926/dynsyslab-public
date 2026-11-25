#pragma once
// Processing steps publish data objects (core/DataObject.hpp). publish() files each archive by its
// object kind and stamps the envelope: what the object is, and the method, settings and inputs that
// produced it. A step that yields several objects publishes the others with publishObject(); the
// first object it returns is its primary output.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

// The processing step being run, for the provenance of what it publishes. Set by the dispatcher.
class StageScope {
public:
    StageScope(std::string method, json::Value settings);
    ~StageScope();
    StageScope(const StageScope&) = delete;
    StageScope& operator=(const StageScope&) = delete;
};

// The object kind an archive published under `name` holds (e.g. "koopman-hankel-dmd" -> eigenfunctions).
[[nodiscard]] std::string objectKindFor(const std::string& name);

// The folder for archives published under `name`: <root>/objects/<kind>.
[[nodiscard]] std::filesystem::path objectFolder(const std::filesystem::path& root, const std::string& name);

// Adds the data_object envelope (unless the archive declares one) for an archive of `kind` read from `inputs`.
void stampDataObject(dynsyslab::core::NpzRecord& archive, const std::string& kind, const std::string& name,
                     const std::vector<std::filesystem::path>& inputs, json::Value properties = {});

// Binding options that compose a Koopman model from data objects: the setting "reconstructions"
// (a path or a list) names reconstructions objects to attach as readouts.
[[nodiscard]] dynsyslab::core::ModelBindingOptions composedBinding(const json::Value& settings);

// Objects a step published besides its returned one; the dispatcher collects them after each step.
void registerAdditionalObject(const std::filesystem::path& path);
[[nodiscard]] std::vector<std::filesystem::path> takeAdditionalObjects();

// After a step: its primary output, then the objects it registered. An archive written without the
// envelope is stamped here, once, before anything reads it.
[[nodiscard]] std::vector<std::filesystem::path> collectStepObjects(const std::filesystem::path& primary,
                                                                    const std::string& method);

}  // namespace dynsyslab::workflow::postprocess_detail
