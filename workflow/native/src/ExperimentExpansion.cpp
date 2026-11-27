#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/workflow/Experiment.hpp>
#include "PipelineFingerprint.hpp"
#include <dynsyslab/core/OutputPath.hpp>

#include <dynsyslab/workflow/Plan.hpp>

#include <dynsyslab/core/Text.hpp>

#include "WorkflowSupport.hpp"
#include "ExperimentJson.hpp"
#include "ExperimentStages.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

// Experiment expansion: parameters, extends chains and sequence stages, before any stage runs
// (EXPERIMENTS.md). Execution and stage reuse are in Experiment.cpp.
namespace dynsyslab::workflow {
using namespace experiment_detail;
namespace {

constexpr const char* kSchema = "dynsyslab.compute-job/v1";
constexpr std::size_t kMaximumExtendsDepth = 8;


std::string stageId(const json::Value& stage) {
    const auto id = stage.at("id").string();
    if (id.empty() || !std::regex_match(id, std::regex("[A-Za-z][A-Za-z0-9_-]*")))
        throw std::invalid_argument("experiment stage ids must start with a letter and use [A-Za-z0-9_-]");
    return id;
}

// Parameters merge like a JSON merge patch except that null is kept as a value, not
// treated as "delete": null is how a toggle parameter (control_limit, target_ramp)
// is switched off in an extending file or with --param NAME=null.
json::Value mergeParameters(json::Value target, const json::Value& patch) {
    const auto* fields = std::get_if<json::Value::Object>(&patch.storage());
    if (!fields) return patch;
    if (!std::holds_alternative<json::Value::Object>(target.storage())) target = json::Value(json::Value::Object{});
    for (const auto& [key, value] : *fields)
        target.object()[key] = value.isNull() ? json::Value() : mergeParameters(target.at(key), value);
    return target;
}

// "${params.NAME}" as a whole string becomes the parameter's JSON value (any
// type); embedded in a longer string it is spliced in as text. Parameter values
// may themselves contain parameter or stage references. "${params.NAME.FIELD...}"
// reads a field of an object parameter (e.g. ${params.B_psi.representation}).
json::Value substituteParameters(const json::Value& value, const json::Value& parameters,
                                 std::size_t depth) {
    if (depth > 16) throw std::invalid_argument("experiment parameters reference each other too deeply");
    static const std::regex pattern(R"(\$\{params\.([A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*)\})");
    // A parameter (or field) that is present with the value null is defined: null is
    // how toggles such as control_limit or target_ramp are switched off. Only a
    // missing key is an error.
    auto lookup = [&](const std::string& path) -> json::Value {
        const auto child = [&](const json::Value& object, const std::string& key) -> json::Value {
            const auto* fields = std::get_if<json::Value::Object>(&object.storage());
            if (!fields || !fields->count(key)) throw std::invalid_argument("undefined experiment parameter: " + path);
            return fields->at(key);
        };
        std::size_t begin = 0, dot = path.find('.');
        json::Value found = child(parameters, path.substr(0, dot));
        while (dot != std::string::npos) {
            begin = dot + 1; dot = path.find('.', begin);
            // A field of a parameter that is itself a reference is read after that reference resolves.
            if (std::holds_alternative<std::string>(found.storage()))
                found = substituteParameters(found, parameters, depth + 1);
            json::Value next = child(found, path.substr(begin, dot == std::string::npos ? std::string::npos : dot - begin));
            found = std::move(next);
        }
        return found;
    };
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        std::smatch match;
        if (std::regex_match(*text, match, pattern))
            return substituteParameters(lookup(match[1].str()), parameters, depth + 1);
        std::string out; auto begin = text->cbegin();
        while (std::regex_search(begin, text->cend(), match, pattern)) {
            out.append(begin, match[0].first);
            const auto replaced = substituteParameters(lookup(match[1].str()), parameters, depth + 1);
            const auto* asText = std::get_if<std::string>(&replaced.storage());
            out += asText ? *asText : json::dump(replaced);
            begin = match[0].second;
        }
        out.append(begin, text->cend());
        return json::Value(out);
    }
    if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        json::Value::Array out;
        for (const auto& item : *array) out.push_back(substituteParameters(item, parameters, depth));
        return json::Value(std::move(out));
    }
    if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        json::Value::Object out;
        for (const auto& [key, item] : *object) out[key] = substituteParameters(item, parameters, depth);
        return json::Value(std::move(out));
    }
    return value;
}

// Methods that take a spectrum from the method run just before them.
bool spectrumConsumer(const std::string& method) {
    return method == "sindy" || method == "boundary_eigenfunctions" || method == "eigenvalue_optimization";
}

// A "sequence" stage is a list of post-processing methods chosen by a
// parameter. Each method becomes its own stage (ID_1, ID_2, ...) reading the
// same sources, so each is reused or rerun on its own; a method that takes a
// spectrum is fed the previous method's archive. A method may be named by a
// string, taking its settings from parameters.method_settings.NAME. ${ID.*}
// refers to the last method's stage.
void expandSequences(json::Value& result) {
    const auto& parameters = result.at("parameters");
    json::Value::Array out;
    auto aliases = result.at("stage_aliases").object();  // an expanded job keeps its aliases
    for (const auto& stage : result.at("stages").array()) {
        if (stage.at("sequence").isNull()) { out.push_back(stage); continue; }
        const auto id = stage.at("id").string();
        const auto& items = stage.at("sequence").array();
        if (items.empty()) throw std::invalid_argument("sequence stage '" + id + "' lists no methods");
        json::Value::Array generated;
        std::string previous;
        for (std::size_t k = 0; k < items.size(); ++k) {
            auto method = items[k];
            if (const auto* name = std::get_if<std::string>(&method.storage())) {
                const auto& preset = parameters.at("method_settings").at(*name);
                method = json::Value(json::Value::Object{{"name", json::Value(*name)},
                    {"settings", preset.isNull() ? json::Value(json::Value::Object{})
                                                 : substituteParameters(preset, parameters, 0)}});
            }
            if (!isObject(method) || method.at("name").string().empty())
                throw std::invalid_argument("sequence stage '" + id + "' needs method names or {name, settings} objects");
            if (!previous.empty() && spectrumConsumer(method.at("name").string())
                && method.at("settings").at("spectrum_archive").isNull())
                method.at("settings").at("spectrum_archive") = json::Value("${" + previous + ".archive}");
            const auto stageIdText = id + "_" + std::to_string(k + 1);
            json::Value job(json::Value::Object{{"schema", json::Value(kSchema)}, {"kind", json::Value("postprocess")},
                {"sources", stage.at("sources")}, {"methods", json::Value(json::Value::Array{method})},
                {"label", json::Value(result.at("label").string("experiment") + " / " + stageIdText)}});
            if (!stage.at("archive").isNull()) job.at("archive") = stage.at("archive");
            out.emplace_back(json::Value::Object{{"id", json::Value(stageIdText)}, {"job", std::move(job)}});
            generated.emplace_back(stageIdText);
            previous = stageIdText;
        }
        aliases[id] = json::Value(std::move(generated));
    }
    result.at("stages") = json::Value(std::move(out));
    result.at("stage_aliases") = json::Value(std::move(aliases));
}

json::Value expandDepth(const json::Value& job, std::size_t depth) {
    if (depth > kMaximumExtendsDepth) throw std::invalid_argument("experiment extends chain is too deep");
    if (job.at("kind").string() != "experiment")
        throw std::invalid_argument("expandExperiment expects kind \"experiment\"");
    json::Value result;
    const auto base = job.at("extends").string();
    if (!base.empty()) {
        result = expandDepth(readJson(std::filesystem::u8path(base)), depth + 1);
        result.at("extends_chain").array().insert(result.at("extends_chain").array().begin(),
                                                  json::Value(base));
    } else {
        result.at("stages").array();
        // An already-expanded job keeps the chain it was expanded from.
        result.at("extends_chain") = job.at("extends_chain").isNull()
            ? json::Value(json::Value::Array{}) : job.at("extends_chain");
    }
    for (const auto& [key, value] : job.object())
        if (key != "stages" && key != "overrides" && key != "extends" && key != "extends_chain"
            && key != "parameters")
            result.at(key) = value;
    if (!job.at("parameters").isNull())
        result.at("parameters") = mergeParameters(result.at("parameters"), job.at("parameters"));
    auto& stages = result.at("stages").array();
    for (const auto& stage : job.at("stages").array()) {
        const auto id = stageId(stage);
        auto existing = std::find_if(stages.begin(), stages.end(),
            [&](const json::Value& item) { return item.at("id").string() == id; });
        if (existing != stages.end()) { *existing = stage; continue; }
        const auto before = stage.at("before").string();
        auto position = before.empty() ? stages.end() : std::find_if(stages.begin(), stages.end(),
            [&](const json::Value& item) { return item.at("id").string() == before; });
        if (!before.empty() && position == stages.end())
            throw std::invalid_argument("stage '" + id + "' is to be inserted before unknown stage '" + before + "'");
        auto inserted = stage; inserted.object().erase("before");
        stages.insert(position, std::move(inserted));
    }
    for (const auto& [id, patch] : job.at("overrides").object()) {
        auto target = std::find_if(stages.begin(), stages.end(),
            [&](const json::Value& item) { return item.at("id").string() == id; });
        if (target == stages.end()) throw std::invalid_argument("override names unknown stage: " + id);
        target->at("job") = mergePatch(target->at("job"), patch);
    }
    std::set<std::string> ids;
    for (auto& stage : stages) {
        if (!ids.insert(stageId(stage)).second)
            throw std::invalid_argument("duplicate experiment stage id: " + stage.at("id").string());
        if (!stage.at("sequence").isNull() || !stage.at("sweep").isNull() || !stage.at("select").isNull())
            continue;  // expanded into method, point and selection stages below
        // A whole job may be a parameter ("job":"${params.training}"); it is checked after substitution.
        if (std::holds_alternative<std::string>(stage.at("job").storage())) continue;
        if (!isObject(stage.at("job"))) throw std::invalid_argument("experiment stage needs a job object");
    }
    if (stages.empty()) throw std::invalid_argument("experiment has no stages");
    if (depth == 0) {  // substitute ${params.NAME} once, after the whole chain merged
        for (auto& stage : stages) {
            if (!stage.at("sequence").isNull()) {
                for (const auto* field : {"sequence", "sources", "archive"})
                    stage.at(field) = substituteParameters(stage.at(field), result.at("parameters"), 0);
            } else if (!stage.at("sweep").isNull() || !stage.at("select").isNull()) {
                const auto* field = stage.at("sweep").isNull() ? "select" : "sweep";
                stage.at(field) = substituteParameters(stage.at(field), result.at("parameters"), 0);
            } else {
                auto& stageJob = stage.at("job");
                stageJob = substituteParameters(stageJob, result.at("parameters"), 0);
                if (!isObject(stageJob)) throw std::invalid_argument("experiment stage '" + stage.at("id").string() + "' needs a job object");
                if (stageJob.at("schema").isNull()) stageJob.at("schema") = json::Value(kSchema);
                if (stageJob.at("kind").string() == "experiment")
                    throw std::invalid_argument("experiment stages cannot themselves be experiments; use extends");
                if (stageJob.at("label").isNull())  // label with the final experiment name
                    stageJob.at("label") = json::Value(result.at("label").string("experiment") + " / "
                                                       + stage.at("id").string());
            }
        }
        expandSweepsAndSelects(result);
        expandSequences(result);
        std::set<std::string> all;
        for (const auto& stage : result.at("stages").array())
            if (!all.insert(stage.at("id").string()).second)
                throw std::invalid_argument("duplicate experiment stage id: " + stage.at("id").string());
    }
    return result;
}

}  // namespace

json::Value mergePatch(json::Value target, const json::Value& patch) {
    if (!isObject(patch)) return patch;
    if (!isObject(target)) target = json::Value(json::Value::Object{});
    for (const auto& [key, value] : patch.object()) {
        if (value.isNull()) target.object().erase(key);
        else target.at(key) = mergePatch(target.at(key), value);
    }
    return target;
}

json::Value expandExperiment(const json::Value& job) { return expandDepth(job, 0); }

}  // namespace dynsyslab::workflow
