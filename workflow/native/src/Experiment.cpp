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

namespace dynsyslab::workflow {
using namespace experiment_detail;
namespace {


// ${id.archive} | ${id.archives[N]} | ${id.output:METHOD} | ${id.object:KIND} | ${id.manifest}
// A method may publish several data objects; its first is its primary output, which ${id.archive}
// and ${id.output:METHOD} name. ${id.object:KIND} names the one object of that kind.
const std::regex& referencePattern() {
    static const std::regex pattern(
        R"(\$\{([A-Za-z][A-Za-z0-9_-]*)\.(archive|manifest|archives\[([0-9]+)\]|output:([A-Za-z0-9_]+)|object:([A-Za-z0-9_]+)|archives)\})");
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
    if (field == "archives")  // only a whole string names the list (resolve)
        throw std::invalid_argument("${" + id + ".archives} must be a whole string; it stands for a list");
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
        // "${id.archives}" as a whole string is the list of the stage's archives.
        if (std::regex_match(*text, match, referencePattern()) && match[2].str() == "archives") {
            const auto found = done.find(match[1].str());
            if (found == done.end())
                throw std::invalid_argument("reference to stage '" + match[1].str() + "' that has not run earlier");
            json::Value::Array list;
            for (const auto& path : found->second.run.archives) list.emplace_back(path.generic_u8string());
            return json::Value(std::move(list));
        }
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
    (void)dynsyslab::core::prepareOutputFile(path, "experiment manifest");
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    { std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to create experiment manifest");
      stream << contents << '\n'; }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish experiment manifest");
#else
    std::filesystem::rename(temporary, path);
#endif
}

// Stage reuse: an earlier experiment manifest whose stage job is byte-identical
// (after full resolution), computed by this same pipeline (its binaries' fingerprint), and whose
// archives all still exist supplies the result. Identity cascades: reusing a stage keeps
// downstream jobs identical. A changed pipeline recomputes, so rerunning an experiment shows it.
std::optional<RunResult> reusableStage(const std::filesystem::path& archiveRoot,
                                       const json::Value& resolved) {
    const auto directory = archiveRoot / L"campaigns" / L"experiments";
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
            if (stage.at("pipeline").string() != detail::pipelineFingerprint()) continue;
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
            else if (path.rfind("__", 0) == 0) throw std::invalid_argument("--until and --refresh apply to experiment jobs");
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
        if (job.at("kind").string() != "experiment") throw std::invalid_argument("--until applies only to experiment jobs");
        record["until"] = json::Value(until);
    }
    if (job.at("kind").string() != "experiment") {
        if (!params.empty()) throw std::invalid_argument("--param applies only to experiment jobs; use --set");
        for (const auto& [path, value] : setRecord) setJobInputPath(job, path, value);
        return job;
    }
    // NAME may be a dotted path into the merged parameters (training.time.stop).
    if (!paramRecord.empty()) {
        auto merged = expandExperiment(job).at("parameters");
        for (const auto& [name, value] : paramRecord) setJobInputPath(merged, name, value);
        job.at("parameters") = merged;
    }
    auto expanded = expandExperiment(job);
    auto& stages = expanded.at("stages").array();
    for (const auto& [path, value] : setRecord) {
        const auto dot = path.find('.');
        if (dot == std::string::npos) throw std::invalid_argument("experiment --set needs STAGE.path: " + path);
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

json::Value planExperiment(const json::Value& job) {
    const auto expanded = expandExperiment(job);
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
        {"schema", json::Value("dynsyslab.workflow-plan/v1")}, {"kind", json::Value("experiment")},
        {"status", json::Value("planned")}, {"stage_count", json::Value(double(stages.size()))},
        {"extends_chain", expanded.at("extends_chain")}, {"stages", json::Value(std::move(stages))}});
}

RunResult executeExperiment(const json::Value& job, const std::filesystem::path& archiveRoot,
                          const std::filesystem::path& progressPath,
                          const dynsyslab::core::ExecutionContext* execution) {
    const auto expanded = expandExperiment(job);
    (void)planExperiment(job);  // reject forward references before any work runs
    const auto& stages = expanded.at("stages").array();
    const auto label = expanded.at("label").string("experiment");
    const auto stamp = detail::timestamp();
    std::map<std::string, StageOutput> done;
    RunResult result{"experiment", L"", {}, {}, {}};
    json::Value::Array record;
    std::wstring summary = L"Native experiment completed\r\n";
    result.manifest = archiveRoot / L"campaigns" / L"experiments"
        / dynsyslab::core::wide(detail::safeName(label, "experiment") + "__" + stamp + ".experiment.json");
    // Published after every stage, so completed stages of a failed or cancelled
    // experiment stay recorded (and reusable) and the failure is diagnosable.
    auto publishManifest = [&](const std::string& status, const std::string& error) {
        json::Value::Object manifest{
            {"schema", json::Value("dynsyslab.experiment/v1")}, {"label", json::Value(label)},
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
        const bool select = resolved.at("kind").string() == "select";
        if (expanded.at("reuse_identical_stages").boolean(false) && !forced && !existing && !select)
            reused = reusableStage(archiveRoot, resolved);
        const auto selection = result.manifest.parent_path()
            / dynsyslab::core::wide(detail::safeName(label, "experiment") + "__" + stamp + "." + id + ".selection.json");
        auto run = existing ? archiveStage(resolved)
                 : select ? runSelect(resolved, selection)
                 : reused ? std::move(*reused) : executeJob(resolved, archiveRoot, {}, execution);
        record.emplace_back(json::Value::Object{
            {"id", json::Value(id)}, {"kind", json::Value(run.kind)}, {"job", resolved},
            {"reused", json::Value(bool(reused))},
            {"pipeline", json::Value(detail::pipelineFingerprint())},  // reused stages passed the same check
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
    detail::atomicProgress(progressPath, 1.0, stages.size(), stages.size(), "Experiment complete");
    result.summary = summary + L"\r\nExperiment manifest:\r\n" + result.manifest.wstring();
    return result;
}

}  // namespace dynsyslab::workflow
