// Runs a flattened job group through executeJob: resume from a ledger, symbolic archive
// references, reference-counted cleanup of intermediates and a free-space floor.
#include <dynsyslab/workflow/JobGroup.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/ArchiveDigest.hpp>
#include <dynsyslab/workflow/DataPack.hpp>
#include <dynsyslab/workflow/Plan.hpp>

#include "WorkflowSupport.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <variant>

namespace dynsyslab::workflow {
namespace {

constexpr double kGiB = 1024.0 * 1024.0 * 1024.0;

struct Finished {
    std::vector<std::filesystem::path> archives;
    std::vector<std::string> methods;
};

std::string selectArchive(const std::string& leaf, const std::string& selector, const Finished& run) {
    const auto text = [](const std::filesystem::path& path) { return path.generic_u8string(); };
    if (selector == "archive") {
        std::vector<std::filesystem::path> results;
        if (run.methods.size() == run.archives.size() && run.archives.size() > 1) {
            for (std::size_t k = 0; k < run.archives.size(); ++k)
                if (run.methods[k] != "coordinate_transform") results.push_back(run.archives[k]);
        } else results = run.archives;
        if (results.size() != 1)
            throw std::invalid_argument("job '" + leaf + "' published " + std::to_string(run.archives.size())
                                        + " archives; name one with archives[N] or output:METHOD");
        return text(results.front());
    }
    std::smatch match;
    static const std::regex indexed(R"(archives\[([0-9]+)\])");
    if (std::regex_match(selector, match, indexed)) {
        const auto index = std::stoul(match[1].str());
        if (index >= run.archives.size())
            throw std::invalid_argument("job '" + leaf + "' has no archive " + match[1].str());
        return text(run.archives[index]);
    }
    if (selector.rfind("output:", 0) == 0) {
        const auto method = selector.substr(7);
        std::vector<std::filesystem::path> found;
        for (std::size_t k = 0; k < run.methods.size() && k < run.archives.size(); ++k)
            if (run.methods[k] == method) found.push_back(run.archives[k]);
        if (found.size() != 1)
            throw std::invalid_argument("job '" + leaf + "' produced " + std::to_string(found.size())
                                        + " archives for method " + method);
        return text(found.front());
    }
    if (selector.rfind("metadata:", 0) == 0) {  // an archive path a stage recorded, e.g. selected_krom
        if (run.archives.size() != 1) throw std::invalid_argument("job '" + leaf + "' must publish one archive for " + selector);
        const auto key = selector.substr(9);
        const auto chosen = json::parse(core::readNpzIndex(run.archives.front()).metadataJson).at(key).string();
        if (chosen.empty()) throw std::invalid_argument("job '" + leaf + "' recorded no " + key);
        return chosen;
    }
    throw std::invalid_argument("job '" + leaf + "': unknown output selector '" + selector + "'");
}

// Every archive a finished job published (all trajectories of a family solve), for ${id.archives}.
json::Value allArchives(const std::string& leaf, const std::map<std::string, Finished>& done) {
    if (!done.count(leaf)) throw std::invalid_argument("reference to job '" + leaf + "', which has not finished");
    json::Value::Array out;
    for (const auto& path : done.at(leaf).archives) out.emplace_back(path.generic_u8string());
    return json::Value(std::move(out));
}

json::Value resolveTokens(const json::Value& value, const std::map<std::string, Finished>& done,
                          json::Value::Object& resolutions) {
    static const std::regex whole(R"(@\{([^|}]+)\|archives\})");
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        std::smatch all;
        if (std::regex_match(*text, all, whole)) {  // ${id.archives}: the list itself
            auto list = allArchives(all[1].str(), done);
            resolutions[*text] = list;
            return list;
        }
        static const std::regex pattern(R"(@\{([^|}]+)\|([^}]+)\})");
        std::string out; auto begin = text->cbegin(); std::smatch match;
        while (std::regex_search(begin, text->cend(), match, pattern)) {
            out.append(begin, match[0].first);
            const auto leaf = match[1].str();
            if (!done.count(leaf)) throw std::invalid_argument("reference to job '" + leaf + "', which has not finished");
            const auto path = selectArchive(leaf, match[2].str(), done.at(leaf));
            resolutions[match[0].str()] = json::Value(path);
            out += path;
            begin = match[0].second;
        }
        out.append(begin, text->cend());
        return json::Value(out);
    }
    if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        json::Value::Array out;
        for (const auto& item : *array) {
            const auto* text = std::get_if<std::string>(&item.storage());
            std::smatch all;
            if (text && std::regex_match(*text, all, whole)) {  // an archive list inside a list (sources): spliced in
                const auto list = resolveTokens(item, done, resolutions);  // held: .array() refers into it
                for (const auto& path : list.array()) out.push_back(path);
                continue;
            }
            out.push_back(resolveTokens(item, done, resolutions));
        }
        return json::Value(std::move(out));
    }
    if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        json::Value::Object out;
        for (const auto& [key, item] : *object) out[key] = resolveTokens(item, done, resolutions);
        return json::Value(std::move(out));
    }
    return value;
}

void writeLedger(const std::filesystem::path& path, const json::Value& ledger) {
    std::filesystem::create_directories(path.parent_path());
    auto partial = path; partial += L".partial";
    { std::ofstream stream(partial, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to write the job ledger: " + path.u8string());
      stream << json::dump(ledger) << '\n'; }
    std::filesystem::rename(partial, path);
}

json::Value readLedger(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) return json::Value(json::Value::Object{});
    std::ostringstream text; text << stream.rdbuf();
    try { return json::parse(text.str()); } catch (const std::exception&) { return json::Value(json::Value::Object{}); }
}

double freeBytes(const std::filesystem::path& path) {
    std::error_code error;
    auto probe = path;
    while (!probe.empty() && !std::filesystem::exists(probe, error) && probe != probe.parent_path()) probe = probe.parent_path();
    const auto space = std::filesystem::space(probe, error);
    return error ? -1.0 : static_cast<double>(space.available);
}

// Largest amount of disk held at once when each intermediate is deleted after its last consumer.
double simulatedPeak(const JobGroupGraph& graph, bool cleanup) {
    std::map<std::string, std::size_t> consumers;
    for (const auto& leaf : graph.leaves) for (const auto& id : leaf.dependencies) ++consumers[id];
    std::map<std::string, double> held;
    double current = 0.0, peak = 0.0;
    for (const auto& leaf : graph.leaves) {
        current += leaf.estimatedBytes; held[leaf.id] = leaf.estimatedBytes;
        peak = std::max(peak, current);
        if (!cleanup) continue;
        for (const auto& id : leaf.dependencies)
            if (--consumers[id] == 0 && !std::any_of(graph.leaves.begin(), graph.leaves.end(),
                    [&](const JobGroupLeaf& l) { return l.id == id && l.keep; })) current -= held[id];
        if (!consumers[leaf.id] && !leaf.keep) current -= held[leaf.id];
    }
    return peak;
}

double minimumFree(const json::Value& job) {
    const auto& runtime = job.at("runtime").at("min_free_gb");
    return (runtime.isNull() ? job.at("cleanup").at("min_free_gb").number(0.0) : runtime.number()) * kGiB;
}

bool cleanupEnabled(const json::Value& job) {
    const auto& runtime = job.at("runtime").at("cleanup");
    return runtime.isNull() ? job.at("cleanup").at("enabled").boolean(false) : runtime.boolean();
}

std::string leafFolder(std::size_t index) {
    char name[16]; std::snprintf(name, sizeof(name), "L%04zu", index + 1);
    return name;
}

}  // namespace

std::filesystem::path jobGroupRoot(const json::Value& job, const std::filesystem::path& archiveRoot) {
    return archiveRoot / L"job-groups" / core::wide(detail::safeName(job.at("label").string("job-group"), "job-group"));
}

json::Value planJobGroup(const json::Value& job) {
    const auto graph = loadJobGroup(job);
    const auto pin = checkDataPackPin(job);
    json::Value::Array jobs;
    double total = 0.0;
    for (const auto& leaf : graph.leaves) {
        json::Value::Object entry{{"id", json::Value(leaf.id)}, {"kind", json::Value(leaf.job.at("kind").string())},
                                  {"file", json::Value(leaf.file.generic_u8string())},
                                  {"estimated_bytes", json::Value(leaf.estimatedBytes)}};
        json::Value::Array depends;
        for (const auto& id : leaf.dependencies) depends.emplace_back(id);
        entry["depends_on"] = json::Value(std::move(depends));
        if (leaf.dependencies.empty()) {
            const auto kind = leaf.job.at("kind").string();
            if (kind != "viewer_figure" && kind != "archive_transform") (void)planJob(leaf.job);
            entry["plan"] = "validated";
        } else entry["plan"] = "deferred until its inputs are published";
        total += leaf.estimatedBytes;
        jobs.emplace_back(std::move(entry));
    }
    return json::Value(json::Value::Object{
        {"schema", json::Value("dynsyslab.workflow-plan/v1")}, {"kind", json::Value("job_group")},
        {"status", json::Value("planned")}, {"job_count", json::Value(double(jobs.size()))},
        {"data_pack", pin}, {"estimated_total_bytes", json::Value(total)},
        {"estimated_peak_bytes", json::Value(simulatedPeak(graph, cleanupEnabled(job)))},
        {"cleanup", json::Value(cleanupEnabled(job))}, {"min_free_gb", json::Value(minimumFree(job) / kGiB)},
        {"jobs", json::Value(std::move(jobs))}});
}

RunResult executeJobGroup(const json::Value& job, const std::filesystem::path& archiveRoot,
                          const std::filesystem::path& progressPath,
                          const core::ExecutionContext* execution) {
    const auto graph = loadJobGroup(job);
    const auto pin = checkDataPackPin(job);  // refuses a wrong release before any work
    const bool cleanup = cleanupEnabled(job);
    const double floor = minimumFree(job);
    const auto root = jobGroupRoot(job, archiveRoot);
    const double peak = simulatedPeak(graph, cleanup);
    if (floor > 0.0) {
        const double available = freeBytes(root);
        if (available >= 0.0 && available - peak < floor)
            throw std::runtime_error("disk preflight: " + std::to_string(available / kGiB) + " GiB free, the run may hold "
                + std::to_string(peak / kGiB) + " GiB, and the floor is " + std::to_string(floor / kGiB) + " GiB");
    }
    const auto ledgerPath = root / L"ledger.json";
    const auto leafRoot = archiveRoot / L"jg" / core::wide(core::sha256Text(root.generic_u8string()).substr(0, 10));
    const bool fresh = job.at("fresh").boolean(false);
    auto previous = fresh ? json::Value(json::Value::Object{}) : readLedger(ledgerPath);
    if (previous.at("schema").string() != "dynsyslab.job-ledger/v1") previous = json::Value(json::Value::Object{});
    const auto& referenceMap = job.at("runtime").at("reference_map");
    const AgreementTolerance tolerance{job.at("runtime").at("relative_tolerance").number(1e-9),
                                       job.at("runtime").at("absolute_tolerance").number(1e-12)};

    std::map<std::string, std::size_t> consumers;
    std::set<std::string> keptLeaves;
    for (const auto& leaf : graph.leaves) {
        for (const auto& id : leaf.dependencies) ++consumers[id];
        if (leaf.keep) keptLeaves.insert(leaf.id);
    }
    std::map<std::string, std::string> identities;
    std::map<std::string, Finished> done;
    std::set<std::string> unavailable;  // failed or blocked
    json::Value::Object entries;
    json::Value ledger(json::Value::Object{
        {"schema", json::Value("dynsyslab.job-ledger/v1")}, {"label", job.at("label")},
        {"job_file", job.at("job_file")}, {"data_pack", pin}, {"status", json::Value("running")}});
    const auto publish = [&](const std::string& status) {
        ledger.at("status") = json::Value(status);
        ledger.at("jobs") = json::Value(entries);
        writeLedger(ledgerPath, ledger);
    };
    std::map<std::string, std::set<std::string>> heldBy;  // archive -> live jobs whose metadata names it
    const auto recordLineage = [&](const std::string& holder, const std::vector<std::filesystem::path>& archives) {
        for (const auto& path : archives) {
            if (path.extension() != L".npz" || !std::filesystem::is_regular_file(core::extendedLengthPath(path))) continue;
            const auto metadata = core::readNpzIndex(path).metadataJson;
            for (const auto& [id, other] : entries) {
                const auto folder = other.at("folder").string();
                if (id != holder && !folder.empty() && metadata.find(folder + "/") != std::string::npos) heldBy[id].insert(holder);
            }
        }
    };
    const auto release = [&](const std::string& id) {  // nothing still to run or still alive needs `id`
        if (!cleanup || keptLeaves.count(id) || !entries.count(id) || consumers[id] != 0) return false;
        auto& entry = entries[id];
        if (entry.at("status").string() != "complete" || entry.at("deleted").boolean(false)) return false;
        for (const auto& holder : heldBy[id])
            if (!entries.count(holder) || !entries[holder].at("deleted").boolean(false)) return false;
        // Windows finishes a delete only after the last handle closes, so a folder can briefly
        // stay non-empty. A folder that stays is recorded and does not fail the science.
        std::error_code error;
        const auto folder = core::extendedLengthPath(std::filesystem::u8path(entry.at("folder").string()));
        for (int attempt = 0; attempt < 25; ++attempt) {
            std::filesystem::remove_all(folder, error);
            if (!error && !std::filesystem::exists(folder)) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
        if (std::filesystem::exists(folder)) { entry.at("cleanup_error") = json::Value(error.message()); return false; }
        entry.at("deleted") = json::Value(true);
        return true;
    };
    const auto sweep = [&] {  // deleting one archive can free the archives it alone was holding
        for (bool changed = true; changed;) {
            changed = false;
            for (const auto& leaf : graph.leaves) changed = release(leaf.id) || changed;
        }
    };

    std::size_t failures = 0;
    std::string firstFailure;
    for (std::size_t index = 0; index < graph.leaves.size(); ++index) {
        if (execution) execution->checkpoint();
        const auto& leaf = graph.leaves[index];
        detail::atomicProgress(progressPath, double(index) / double(graph.leaves.size()), index,
                               graph.leaves.size(), "Running job " + leaf.id);
        std::string identityText = json::dump(leaf.job);
        for (const auto& id : leaf.dependencies) identityText += "|" + identities[id];
        identities[leaf.id] = core::sha256Text(identityText);

        const auto blocker = std::find_if(leaf.dependencies.begin(), leaf.dependencies.end(),
            [&](const std::string& id) { return unavailable.count(id) > 0; });
        if (blocker != leaf.dependencies.end()) {
            unavailable.insert(leaf.id);
            entries[leaf.id] = json::Value(json::Value::Object{{"status", json::Value("blocked")},
                {"file", json::Value(leaf.file.generic_u8string())}, {"blocked_by", json::Value(*blocker)}});
            publish("running");
            continue;
        }
        // Resume: a finished job with the same identity is kept. Its files may already have been
        // cleaned; they are recomputed only if a job still to run needs them.
        const auto& before = previous.at("jobs").at(leaf.id);
        bool reused = before.at("status").string() == "complete"
            && before.at("identity").string() == identities[leaf.id];
        if (reused) {
            Finished finished;
            bool present = true;
            for (const auto& item : before.at("archives").array()) {
                finished.archives.push_back(std::filesystem::u8path(item.string()));
                present = present && std::filesystem::is_regular_file(core::extendedLengthPath(finished.archives.back()));
            }
            for (const auto& item : before.at("methods").array()) finished.methods.push_back(item.string());
            bool neededLater = false;
            for (std::size_t later = index + 1; later < graph.leaves.size() && !neededLater; ++later) {
                const auto& other = graph.leaves[later];
                if (!other.dependencies.count(leaf.id)) continue;
                const auto& state = previous.at("jobs").at(other.id);
                neededLater = state.at("status").string() != "complete";
            }
            if (present || !neededLater) {
                entries[leaf.id] = before;
                entries[leaf.id].at("reused") = json::Value(true);
                recordLineage(leaf.id, finished.archives);
                done[leaf.id] = std::move(finished);
                for (const auto& id : leaf.dependencies) --consumers[id];
                sweep();
                publish("running");
                continue;
            }
        }
        const auto folder = leafRoot / core::wide(leafFolder(index));
        std::error_code ignored;
        std::filesystem::remove_all(core::extendedLengthPath(folder), ignored);  // a partial earlier attempt
        const double available = freeBytes(root);
        if (floor > 0.0 && available >= 0.0 && available - leaf.estimatedBytes < floor) {
            publish("failed");
            throw std::runtime_error("disk floor: job '" + leaf.id + "' would leave less than "
                                     + std::to_string(floor / kGiB) + " GiB free");
        }
        const auto started = std::chrono::steady_clock::now();
        json::Value::Object entry{{"status", json::Value("complete")},
                                  {"file", json::Value(leaf.file.generic_u8string())},
                                  {"identity", json::Value(identities[leaf.id])},
                                  {"folder", json::Value(folder.generic_u8string())}};
        try {
            json::Value::Object resolutions;
            auto resolved = resolveTokens(leaf.job, done, resolutions);
            const auto kind = resolved.at("kind").string();
            if (kind == "viewer_figure") resolved.at("output_directory") = json::Value((root / L"figures").generic_u8string());
            auto run = executeJob(resolved, folder, {}, execution);
            Finished finished{run.archives, run.methods};
            json::Value::Array archives, methods, digests;
            for (const auto& path : run.archives) {
                archives.emplace_back(path.generic_u8string());
                digests.push_back(resultDigest(path));
            }
            for (const auto& method : run.methods) methods.emplace_back(method);
            entry["archives"] = json::Value(std::move(archives));
            entry["methods"] = json::Value(std::move(methods));
            entry["digests"] = json::Value(std::move(digests));
            entry["inputs"] = json::Value(std::move(resolutions));
            const auto& references = referenceMap.at(leaf.id);
            if (!references.isNull()) {
                json::Value::Array compared;
                const auto& list = references.array();
                for (std::size_t k = 0; k < list.size() && k < run.archives.size(); ++k)
                    if (!list[k].isNull())
                        compared.push_back(compareResults(run.archives[k], std::filesystem::u8path(list[k].string()), tolerance));
                // {"derived/some-folder":[reference,...]}: intermediates a stage wrote but did not
                // publish, matched in the order they were written.
                for (const auto& [subfolder, expected] : references.object()) {
                    std::vector<std::filesystem::path> written;
                    const auto directory = core::extendedLengthPath(folder / std::filesystem::u8path(subfolder));
                    if (std::filesystem::is_directory(directory))
                        for (const auto& item : std::filesystem::directory_iterator(directory))
                            if (item.path().extension() == L".npz") written.push_back(item.path());
                    std::sort(written.begin(), written.end(), [](const auto& a, const auto& b) {
                        return std::filesystem::last_write_time(a) < std::filesystem::last_write_time(b); });
                    if (written.size() != expected.array().size())
                        throw std::runtime_error("reference map: " + subfolder + " holds " + std::to_string(written.size()) + " archives");
                    for (std::size_t k = 0; k < written.size(); ++k) {
                        auto item = compareResults(written[k], std::filesystem::u8path(expected.array()[k].string()), tolerance);
                        item.at("reference") = expected.array()[k];
                        compared.push_back(std::move(item));
                    }
                }
                entry["reference_comparison"] = json::Value(std::move(compared));
            }
            entry["seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            done[leaf.id] = std::move(finished);
            entries[leaf.id] = json::Value(std::move(entry));
        } catch (const core::ExecutionCancelled&) {
            publish("cancelled");
            throw;
        } catch (const std::exception& error) {
            ++failures;
            if (firstFailure.empty()) firstFailure = "job '" + leaf.id + "' (" + leaf.file.u8string() + "): " + error.what();
            unavailable.insert(leaf.id);
            entry["status"] = "failed"; entry["error"] = std::string(error.what());
            entries[leaf.id] = json::Value(std::move(entry));
            publish("running");
            continue;
        }
        recordLineage(leaf.id, done[leaf.id].archives);
        for (const auto& id : leaf.dependencies) --consumers[id];
        sweep();
        publish("running");
    }
    if (failures) {
        publish("failed");
        throw std::runtime_error("job group failed: " + firstFailure + " (" + std::to_string(failures) + " failed)");
    }
    std::error_code emptied;  // removes the short leaf folder once cleanup has emptied it
    if (cleanup) std::filesystem::remove(leafRoot, emptied);
    publish("complete");
    detail::atomicProgress(progressPath, 1.0, graph.leaves.size(), graph.leaves.size(), "Job group complete");
    RunResult result{"job_group", L"", {}, ledgerPath, {}};
    for (const auto& leaf : graph.leaves)
        if (leaf.job.at("kind").string() == "viewer_figure" || keptLeaves.count(leaf.id))
            result.archives.insert(result.archives.end(), done[leaf.id].archives.begin(), done[leaf.id].archives.end());
    result.summary = L"Job group completed: " + std::to_wstring(graph.leaves.size()) + L" jobs\r\nLedger:\r\n" + ledgerPath.wstring();
    return result;
}

}  // namespace dynsyslab::workflow
