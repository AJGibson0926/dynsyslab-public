#pragma once
// Processing steps publish data objects (core/DataObject.hpp). Every writer declares the object's
// kind (a core::object_kind); publishObject() files the archive by it and stamps the envelope: what
// the object is, and the method, settings and inputs that produced it. A step that yields several
// objects registers the others (registerAdditionalObject); the first it returns is its primary output.
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

// The method of the step being run (its stage name, e.g. "hdmd"); empty outside a step. Archives
// record it as their "method" and as data_object.produced_by.method, so the two never differ.
[[nodiscard]] std::string currentStageMethod();

// The folder for objects of `kind` (a core::object_kind): <root>/objects/<kind>.
[[nodiscard]] std::filesystem::path kindFolder(const std::filesystem::path& root, const std::string& kind);

// Metadata text with the data_object envelope added: for writers that hand metadata to a core saver
// (saveKoopmanModel, makeKoopmanArchive) instead of publishing an archive themselves.
[[nodiscard]] std::string withDataObject(const std::string& metadataJson, const std::string& kind,
                                         const std::vector<std::filesystem::path>& inputs);

// Adds the data_object envelope (unless the archive declares one) for an archive of `kind` read from
// `inputs`; with none, from the archives its metadata names (source_archive, source_krom, source_spectrum).
void stampDataObject(dynsyslab::core::NpzRecord& archive, const std::string& kind, const std::string& name,
                     const std::vector<std::filesystem::path>& inputs, json::Value properties = {});

// Binding options from a step's settings: "eigenfunction_representation" names the representation
// the K-ROM evaluates.
[[nodiscard]] dynsyslab::core::ModelBindingOptions composedBinding(const json::Value& settings);

// Objects a step published besides its returned one; the dispatcher collects them after each step.
void registerAdditionalObject(const std::filesystem::path& path);
[[nodiscard]] std::vector<std::filesystem::path> takeAdditionalObjects();

// After a step: its primary output, then the objects it registered. Every object is stamped when it is
// written (publishObject, stampDataObject, withDataObject); one without its envelope is an error.
[[nodiscard]] std::vector<std::filesystem::path> collectStepObjects(const std::filesystem::path& primary,
                                                                    const std::string& method);

}  // namespace dynsyslab::workflow::postprocess_detail
