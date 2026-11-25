#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/workflow/Pipeline.hpp>
#include <dynsyslab/core/OutputPath.hpp>

#include <dynsyslab/workflow/Plan.hpp>

#include <dynsyslab/core/Text.hpp>

#include "WorkflowSupport.hpp"

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

namespace dynsyslab::workflow {
namespace {

constexpr const char* kSchema = "dynsyslab.compute-job/v1";
constexpr std::size_t kMaximumExtendsDepth = 8;

bool isObject(const json::Value& value) {
    return std::holds_alternative<json::Value::Object>(value.storage());
}

json::Value readJson(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::invalid_argument("pipeline base is unavailable: " + path.u8string());
    std::stringstream buffer; buffer << stream.rdbuf();
    return json::parse(buffer.str());
}

std::string stageId(const json::Value& stage) {
    const auto id = stage.at("id").string();
    if (id.empty() || !std::regex_match(id, std::regex("[A-Za-z][A-Za-z0-9_-]*")))
        throw std::invalid_argument("pipeline stage ids must start with a letter and use [A-Za-z0-9_-]");
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
    if (depth > 16) throw std::invalid_argument("pipeline parameters reference each other too deeply");
    static const std::regex pattern(R"(\$\{params\.([A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*)\})");
    // A parameter (or field) that is present with the value null is defined: null is
    // how toggles such as control_limit or target_ramp are switched off. Only a
    // missing key is an error.
    auto lookup = [&](const std::string& path) -> json::Value {
        const auto child = [&](const json::Value& object, const std::string& key) -> json::Value {
            const auto* fields = std::get_if<json::Value::Object>(&object.storage());
            if (!fields || !fields->count(key)) throw std::invalid_argument("undefined pipeline parameter: " + path);
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
                {"label", json::Value(result.at("label").string("pipeline") + " / " + stageIdText)}});
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
    if (depth > kMaximumExtendsDepth) throw std::invalid_argument("pipeline extends chain is too deep");
    if (job.at("kind").string() != "pipeline")
        throw std::invalid_argument("expandPipeline expects kind \"pipeline\"");
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
            throw std::invalid_argument("duplicate pipeline stage id: " + stage.at("id").string());
        if (!stage.at("sequence").isNull()) continue;  // expanded into method stages below
        // A whole job may be a parameter ("job":"${params.training}"); it is checked after substitution.
        if (std::holds_alternative<std::string>(stage.at("job").storage())) continue;
        if (!isObject(stage.at("job"))) throw std::invalid_argument("pipeline stage needs a job object");
    }
    if (stages.empty()) throw std::invalid_argument("pipeline has no stages");
    if (depth == 0) {  // substitute ${params.NAME} once, after the whole chain merged
        for (auto& stage : stages) {
            if (!stage.at("sequence").isNull()) {
                for (const auto* field : {"sequence", "sources", "archive"})
                    stage.at(field) = substituteParameters(stage.at(field), result.at("parameters"), 0);
            } else {
                auto& stageJob = stage.at("job");
                stageJob = substituteParameters(stageJob, result.at("parameters"), 0);
                if (!isObject(stageJob)) throw std::invalid_argument("pipeline stage '" + stage.at("id").string() + "' needs a job object");
                if (stageJob.at("schema").isNull()) stageJob.at("schema") = json::Value(kSchema);
                if (stageJob.at("kind").string() == "pipeline")
                    throw std::invalid_argument("pipeline stages cannot themselves be pipelines; use extends");
                if (stageJob.at("label").isNull())  // label with the final pipeline name
                    stageJob.at("label") = json::Value(result.at("label").string("pipeline") + " / "
                                                       + stage.at("id").string());
            }
        }
        expandSequences(result);
        std::set<std::string> all;
        for (const auto& stage : result.at("stages").array())
            if (!all.insert(stage.at("id").string()).second)
                throw std::invalid_argument("duplicate pipeline stage id: " + stage.at("id").string());
    }
    return result;
}

// ${id.archive} | ${id.archives[N]} | ${id.output:METHOD} | ${id.object:KIND} | ${id.manifest}
// A method may publish several data objects; its first is its primary output, which ${id.archive}
// and ${id.output:METHOD} name. ${id.object:KIND} names the one object of that kind.
const std::regex& referencePattern() {
    static const std::regex pattern(
        R"(\$\{([A-Za-z][A-Za-z0-9_-]*)\.(archive|manifest|archives\[([0-9]+)\]|output:([A-Za-z0-9_]+)|object:([A-Za-z0-9_]+))\})");
    return pattern;
}

struct StageOutput {
    RunResult run;
};

std::string resolveOne(const std::smatch& match, const std::map<std::string, StageOutput>& done) {
    const auto id = match[1].str();
    const auto found = done.find(id);
    if (found == done.end())
        throw std::invalid_argument("reference to stage '" + id + "' that has not run earlier");
    const auto& run = found->second.run;
    const auto field = match[2].str();
    auto text = [](const std::filesystem::path& path) { return path.generic_u8string(); };
    if (field == "manifest") {
        if (run.manifest.empty()) throw std::invalid_argument("stage '" + id + "' published no manifest");
        return text(run.manifest);
    }
    if (field == "archive") {
        // A method that reads observables gets an automatic coordinate change first; that
        // intermediate archive is not the stage's result, so ${id.archive} names the one other.
        if (run.archives.size() > 1 && run.methods.size() == run.archives.size()) {
            std::vector<std::filesystem::path> results;  // each method's primary object
            for (std::size_t k = 0; k < run.archives.size(); ++k)
                if (run.methods[k] != "coordinate_transform" && (k == 0 || run.methods[k] != run.methods[k - 1]))
                    results.push_back(run.archives[k]);
            if (results.size() == 1) return text(results.front());
        }
        if (run.archives.size() != 1)
            throw std::invalid_argument("stage '" + id + "' published " + std::to_string(run.archives.size())
                + " archives; use ${" + id + ".archives[N]} or ${" + id + ".output:METHOD}");
        return text(run.archives.front());
    }
    if (match[3].matched) {
        const auto index = std::stoul(match[3].str());
        if (index >= run.archives.size())
            throw std::invalid_argument("stage '" + id + "' has no archive " + match[3].str());
        return text(run.archives[index]);
    }
    if (match[5].matched) {
        const auto kind = match[5].str();
        std::vector<std::filesystem::path> objects;
        for (const auto& path : run.archives)
            if (dynsyslab::core::dataObjectKind(dynsyslab::core::readNpzIndex(path).metadataJson) == kind)
                objects.push_back(path);
        if (objects.size() != 1)
            throw std::invalid_argument("stage '" + id + "' produced " + std::to_string(objects.size())
                + " objects of kind " + kind);
        return text(objects.front());
    }
    const auto method = match[4].str();
    std::vector<std::filesystem::path> matches;
    for (std::size_t k = 0; k < run.methods.size() && k < run.archives.size(); ++k)
        if (run.methods[k] == method && (k == 0 || run.methods[k - 1] != method))  // its primary object
            matches.push_back(run.archives[k]);
    if (matches.size() != 1)
        throw std::invalid_argument("stage '" + id + "' produced " + std::to_string(matches.size())
            + " archives for method " + method);
    return text(matches.front());
}

json::Value resolve(const json::Value& value, const std::map<std::string, StageOutput>& done) {
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        std::string out; auto begin = text->cbegin(); std::smatch match;
        while (std::regex_search(begin, text->cend(), match, referencePattern())) {
            out.append(begin, match[0].first); out += resolveOne(match, done);
            begin = match[0].second;
        }
        out.append(begin, text->cend());
        return json::Value(out);
    }
    if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        json::Value::Array out; for (const auto& item : *array) out.push_back(resolve(item, done));
        return json::Value(std::move(out));
    }
    if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        json::Value::Object out;
        for (const auto& [key, item] : *object) out[key] = resolve(item, done);
        return json::Value(std::move(out));
    }
    return value;
}

void collectReferences(const json::Value& value, std::set<std::string>& ids) {
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        for (std::sregex_iterator it(text->begin(), text->end(), referencePattern()), end; it != end; ++it)
            ids.insert((*it)[1].str());
    } else for (const auto& item : value.array()) collectReferences(item, ids);
    for (const auto& [key, item] : value.object()) { (void)key; collectReferences(item, ids); }
}

// The stage a sequence id stands for (its last method), or the id itself.
std::vector<std::string> aliasMembers(const json::Value& expanded, const std::string& id) {
    std::vector<std::string> members;
    for (const auto& item : expanded.at("stage_aliases").at(id).array()) members.push_back(item.string());
    return members;
}

// A stage kind "archive" names existing archives instead of computing them,
// e.g. a training trajectory solved earlier: {"kind":"archive","archives":[path]}.
RunResult archiveStage(const json::Value& job) {
    RunResult run{"archive", L"Existing archives", {}, {}, {}};
    for (const auto& item : job.at("archives").array()) {
        const auto path = std::filesystem::u8path(item.string());
        if (!std::filesystem::is_regular_file(dynsyslab::core::extendedLengthPath(path)))
            throw std::invalid_argument("archive stage names a missing file: " + item.string());
        run.archives.push_back(path);
        run.methods.push_back("archive");
    }
    if (run.archives.empty()) throw std::invalid_argument("archive stage needs a nonempty archives list");
    return run;
}

void atomicText(const std::filesystem::path& path, const std::string& contents) {
    (void)dynsyslab::core::prepareOutputFile(path, "pipeline manifest");
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    { std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to create pipeline manifest");
      stream << contents << '\n'; }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish pipeline manifest");
#else
    std::filesystem::rename(temporary, path);
#endif
}

// Stage reuse: an earlier pipeline manifest whose stage job is byte-identical
// (after full resolution) and whose archives all still exist supplies the
// result. Identity cascades: reusing a stage keeps downstream jobs identical.
std::optional<RunResult> reusableStage(const std::filesystem::path& archiveRoot,
                                       const json::Value& resolved) {
    const auto directory = archiveRoot / L"campaigns" / L"pipelines";
    if (!std::filesystem::is_directory(directory)) return std::nullopt;
    // Labels only name archives; they never change numerical work.
    auto identity = [](json::Value job) { job.object().erase("label"); return json::dump(job); };
    const auto wanted = identity(resolved);
    std::vector<std::filesystem::path> manifests;
    for (const auto& entry : std::filesystem::directory_iterator(directory))
        if (entry.path().extension() == L".json") manifests.push_back(entry.path());
    std::sort(manifests.rbegin(), manifests.rend());  // newest first (UTC-stamped names)
    for (const auto& path : manifests) {
        json::Value manifest;
        try { manifest = readJson(path); } catch (const std::exception&) { continue; }
        for (const auto& stage : manifest.at("stages").array()) {  // only completed stages are recorded
            if (identity(stage.at("job")) != wanted) continue;
            RunResult run{stage.at("kind").string(), L"reused from " + path.wstring(), {}, {}, {}};
            bool present = true;
            for (const auto& item : stage.at("archives").array()) {
                run.archives.push_back(std::filesystem::u8path(item.string()));
                present = present && std::filesystem::is_regular_file(dynsyslab::core::extendedLengthPath(run.archives.back()));
            }
            for (const auto& item : stage.at("methods").array()) run.methods.push_back(item.string());
            run.manifest = std::filesystem::u8path(stage.at("manifest").string());
            if (present && !run.archives.empty()) return run;
        }
    }
    return std::nullopt;
}

json::Value pathArray(const std::vector<std::filesystem::path>& paths) {
    json::Value::Array out; for (const auto& path : paths) out.emplace_back(path.generic_u8string());
    return json::Value(std::move(out));
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

json::Value expandPipeline(const json::Value& job) { return expandDepth(job, 0); }

namespace {
json::Value inputValue(const std::string& text) {
    try { return json::parse(text); } catch (const std::exception&) { return json::Value(text); }
}
void setJobInputPath(json::Value& root, const std::string& path, json::Value value) {
    if (path.empty()) throw std::invalid_argument("--set needs a non-empty path");
    json::Value* node = &root; std::size_t begin = 0;
    while (true) {
        const auto dot = path.find('.', begin);
        const auto segment = path.substr(begin, dot == std::string::npos ? std::string::npos : dot - begin);
        if (segment.empty()) throw std::invalid_argument("invalid --set path: " + path);
        json::Value* next;
        if (std::holds_alternative<json::Value::Array>(node->storage())) {
            if (segment.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument("--set path indexes an array with '" + segment + "': " + path);
            const auto index = std::stoul(segment);
            if (index >= node->array().size())
                throw std::invalid_argument("--set path index out of range: " + path);
            next = &node->array()[index];
        } else next = &node->at(segment);
        if (dot == std::string::npos) { *next = std::move(value); return; }
        node = next; begin = dot + 1;
    }
}
}  // namespace

json::Value applyJobInputs(json::Value job, const std::vector<std::pair<std::string, std::string>>& params,
                           const std::vector<std::pair<std::string, std::string>>& sets) {
    if (params.empty() && sets.empty()) return job;
    if (job.at("kind").string() == "job_group") {
        for (const auto& [name, text] : params) job.at("parameters").at(name) = inputValue(text);
        for (const auto& [path, text] : sets) {
            if (path == "__fresh") job.at("fresh") = json::Value(true);
            else if (path.rfind("__", 0) == 0) throw std::invalid_argument("--until and --refresh apply to pipeline jobs");
            else setJobInputPath(job, path, inputValue(text));
        }
        return job;
    }
    json::Value::Object record, paramRecord, setRecord;
    std::string until;  // "__until" carries --until STAGE
    for (const auto& [name, text] : params) paramRecord[name] = inputValue(text);
    for (const auto& [path, text] : sets) {
        if (path == "__until") until = text;
        else if (path == "__fresh") record["fresh"] = json::Value(true);
        else if (path == "__refresh") record["refresh"].array().emplace_back(text);
        else setRecord[path] = inputValue(text);
    }
    record["params"] = json::Value(paramRecord); record["sets"] = json::Value(setRecord);
    if (!until.empty()) {
        if (job.at("kind").string() != "pipeline") throw std::invalid_argument("--until applies only to pipeline jobs");
        record["until"] = json::Value(until);
    }
    if (job.at("kind").string() != "pipeline") {
        if (!params.empty()) throw std::invalid_argument("--param applies only to pipeline jobs; use --set");
        for (const auto& [path, value] : setRecord) setJobInputPath(job, path, value);
        return job;
    }
    // NAME may be a dotted path into the merged parameters (training.time.stop).
    if (!paramRecord.empty()) {
        auto merged = expandPipeline(job).at("parameters");
        for (const auto& [name, value] : paramRecord) setJobInputPath(merged, name, value);
        job.at("parameters") = merged;
    }
    auto expanded = expandPipeline(job);
    auto& stages = expanded.at("stages").array();
    for (const auto& [path, value] : setRecord) {
        const auto dot = path.find('.');
        if (dot == std::string::npos) throw std::invalid_argument("pipeline --set needs STAGE.path: " + path);
        const auto id = path.substr(0, dot);
        auto stage = std::find_if(stages.begin(), stages.end(),
            [&](const json::Value& s) { return s.at("id").string() == id; });
        if (stage == stages.end()) throw std::invalid_argument("--set names unknown stage: " + id);
        setJobInputPath(stage->at("job"), path.substr(dot + 1), value);
    }
    if (!until.empty()) {
        if (const auto members = aliasMembers(expanded, until); !members.empty()) until = members.back();
        auto last = std::find_if(stages.begin(), stages.end(),
            [&](const json::Value& s) { return s.at("id").string() == until; });
        if (last == stages.end()) throw std::invalid_argument("--until names unknown stage: " + until);
        stages.erase(last + 1, stages.end());
    }
    if (record.count("fresh")) expanded.at("reuse_identical_stages") = json::Value(false);
    if (record.count("refresh")) {
        json::Value::Array refresh;  // a sequence id refreshes every method in it
        for (const auto& id : record["refresh"].array()) {
            const auto members = aliasMembers(expanded, id.string());
            if (members.empty()) refresh.push_back(id);
            for (const auto& member : members) refresh.emplace_back(member);
        }
        record["refresh"] = json::Value(std::move(refresh));
        for (const auto& id : record["refresh"].array())
            if (std::none_of(stages.begin(), stages.end(),
                    [&](const json::Value& s) { return s.at("id").string() == id.string(); }))
                throw std::invalid_argument("--refresh names unknown stage: " + id.string());
        expanded.at("refresh_stages") = record["refresh"];
    }
    expanded.at("command_inputs") = json::Value(std::move(record));
    return expanded;
}

json::Value planPipeline(const json::Value& job) {
    const auto expanded = expandPipeline(job);
    json::Value::Array stages; std::set<std::string> earlier;
    for (const auto& stage : expanded.at("stages").array()) {
        const auto id = stage.at("id").string();
        std::set<std::string> references; collectReferences(stage.at("job"), references);
        for (const auto& reference : references)
            if (!earlier.count(reference))
                throw std::invalid_argument("stage '" + id + "' references '" + reference
                                            + "', which is not an earlier stage");
        json::Value::Object entry{{"id", json::Value(id)},
                                  {"kind", json::Value(stage.at("job").at("kind").string())}};
        if (stage.at("job").at("kind").string() == "archive") entry["plan"] = json::Value("existing archives");
        else if (references.empty()) entry["plan"] = planJob(stage.at("job"));
        else {
            json::Value::Array names; for (const auto& r : references) names.emplace_back(r);
            entry["depends_on"] = json::Value(std::move(names));
            entry["plan"] = json::Value("deferred until its inputs are published");
        }
        stages.emplace_back(std::move(entry)); earlier.insert(id);
        for (const auto& [alias, members] : expanded.at("stage_aliases").object())
            if (!members.array().empty() && members.array().back().string() == id) earlier.insert(alias);
    }
    return json::Value(json::Value::Object{
        {"schema", json::Value("dynsyslab.workflow-plan/v1")}, {"kind", json::Value("pipeline")},
        {"status", json::Value("planned")}, {"stage_count", json::Value(double(stages.size()))},
        {"extends_chain", expanded.at("extends_chain")}, {"stages", json::Value(std::move(stages))}});
}

RunResult executePipeline(const json::Value& job, const std::filesystem::path& archiveRoot,
                          const std::filesystem::path& progressPath,
                          const dynsyslab::core::ExecutionContext* execution) {
    const auto expanded = expandPipeline(job);
    (void)planPipeline(job);  // reject forward references before any work runs
    const auto& stages = expanded.at("stages").array();
    const auto label = expanded.at("label").string("pipeline");
    const auto stamp = detail::timestamp();
    std::map<std::string, StageOutput> done;
    RunResult result{"pipeline", L"", {}, {}, {}};
    json::Value::Array record;
    std::wstring summary = L"Native pipeline completed\r\n";
    result.manifest = archiveRoot / L"campaigns" / L"pipelines"
        / dynsyslab::core::wide(detail::safeName(label, "pipeline") + "__" + stamp + ".pipeline.json");
    // Published after every stage, so completed stages of a failed or cancelled
    // pipeline stay recorded (and reusable) and the failure is diagnosable.
    auto publishManifest = [&](const std::string& status, const std::string& error) {
        json::Value::Object manifest{
            {"format", json::Value("dynsyslab.pipeline/v1")}, {"label", json::Value(label)},
            {"created_utc", json::Value(stamp)}, {"status", json::Value(status)},
            {"extends_chain", expanded.at("extends_chain")},
            {"expanded_job", expanded}, {"stages", json::Value(record)}};
        if (!error.empty()) manifest["error"] = json::Value(error);
        atomicText(result.manifest, json::dump(json::Value(std::move(manifest))));
    };
    std::string current;
    try {
    for (std::size_t index = 0; index < stages.size(); ++index) {
        if (execution) execution->checkpoint();
        const auto id = stages[index].at("id").string(); current = id;
        detail::atomicProgress(progressPath, double(index) / double(stages.size()), index,
                               stages.size(), "Running stage " + id);
        const auto resolved = resolve(stages[index].at("job"), done);
        std::optional<RunResult> reused;
        const auto& refresh = expanded.at("refresh_stages").array();
        const bool forced = std::any_of(refresh.begin(), refresh.end(),
            [&](const json::Value& item) { return item.string() == id; });
        const bool existing = resolved.at("kind").string() == "archive";
        if (expanded.at("reuse_identical_stages").boolean(false) && !forced && !existing)
            reused = reusableStage(archiveRoot, resolved);
        auto run = existing ? archiveStage(resolved)
                 : reused ? std::move(*reused) : executeJob(resolved, archiveRoot, {}, execution);
        record.emplace_back(json::Value::Object{
            {"id", json::Value(id)}, {"kind", json::Value(run.kind)}, {"job", resolved},
            {"reused", json::Value(bool(reused))},
            {"archives", pathArray(run.archives)},
            {"methods", [&] { json::Value::Array m; for (const auto& x : run.methods) m.emplace_back(x);
                              return json::Value(std::move(m)); }()},
            {"manifest", json::Value(run.manifest.generic_u8string())}});
        summary += L"\r\n[" + dynsyslab::core::wide(id) + L"] " + run.summary + L"\r\n";
        result.archives.insert(result.archives.end(), run.archives.begin(), run.archives.end());
        for (std::size_t k = 0; k < run.archives.size(); ++k)
            result.methods.push_back(k < run.methods.size() ? id + ":" + run.methods[k] : id);
        done[id] = {std::move(run)};
        for (const auto& [alias, members] : expanded.at("stage_aliases").object())
            if (!members.array().empty() && members.array().back().string() == id) done[alias] = done[id];
        publishManifest(index + 1 == stages.size() ? "complete" : "running", {});
    }
    } catch (const std::exception& error) {
        publishManifest("failed", "stage '" + current + "': " + error.what());
        throw;
    }
    detail::atomicProgress(progressPath, 1.0, stages.size(), stages.size(), "Pipeline complete");
    result.summary = summary + L"\r\nPipeline manifest:\r\n" + result.manifest.wstring();
    return result;
}

}  // namespace dynsyslab::workflow
