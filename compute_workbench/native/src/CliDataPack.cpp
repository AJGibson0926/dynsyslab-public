#include "CliDataPack.hpp"

#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/DataPack.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>
#include <dynsyslab/workflow/Pipeline.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Queue.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::cli {
namespace {

namespace wf = dynsyslab::workflow;
namespace json = dynsyslab::core::json;

json::Value readJson(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to open: " + path.u8string());
    std::ostringstream text; text << stream.rdbuf();
    return json::parse(text.str());
}

json::Value rootJob(const std::filesystem::path& pack, const Arguments& args) {
    const auto file = wf::dataPackRootJob(pack);
    auto job = wf::applyJobInputs(wf::bindJobFile(readJson(file), file), args.params, args.sets);
    applyJobGroupRuntime(job, args);
    return job;
}

void print(const json::Value& value) { std::cout << json::dump(value) << '\n'; }

}  // namespace

void applyJobGroupRuntime(json::Value& job, const Arguments& args) {
    if (job.at("kind").string() != "job_group") return;
    auto& runtime = job.at("runtime");
    if (args.allowVersionMismatch) runtime.at("allow_version_mismatch") = json::Value(true);
    if (args.noCleanup) runtime.at("cleanup") = json::Value(false);
    if (args.minFreeGb >= 0.0) runtime.at("min_free_gb") = json::Value(args.minFreeGb);
    if (!args.referenceMap.empty()) runtime.at("reference_map") = readJson(args.referenceMap);
    if (runtime.isNull()) job.object().erase("runtime");
}

int dataPackCommand(const Arguments& args) {
    if (args.positional.size() < 3)
        throw std::invalid_argument("datapack needs an action and a pack folder (see docs/DATA_PACKS.md)");
    const auto action = core::narrow(args.positional[1]);
    const std::filesystem::path target(args.positional[2]);
    if (action == "inspect") { print(wf::inspectDataPack(target)); return 0; }
    if (action == "trace") { print(wf::traceArchiveLineage(target)); return 0; }
    if (action == "rehash") { wf::refreshDataPackHashes(target); print(wf::validateDataPack(target)); return 0; }
    if (action == "validate") {
        const auto report = wf::validateDataPack(target);
        print(report);
        return report.at("valid").boolean() ? 0 : 1;
    }
    if (action == "plan") { print(wf::planJob(rootJob(target, args))); return 0; }
    const auto root = core::prepareOutputDirectory(args.archiveRoot, "--archive-root");
    const auto ledger = [&] { return wf::jobGroupRoot(rootJob(target, args), root) / L"ledger.json"; };
    if (action == "run" || action == "verify") {
        if (action == "run" || !std::filesystem::is_regular_file(ledger()))
            (void)wf::executeJob(rootJob(target, args), root);
        if (action == "run") { print(json::Value(json::Value::Object{{"ledger", json::Value(ledger().generic_u8string())}})); return 0; }
        const auto report = wf::verifyDataPackRun(target, ledger());
        print(report);
        const auto verdict = report.at("verdict").string();
        return verdict == "verified_exact" ? 0 : verdict == "numerically_equivalent" ? 3 : 1;
    }
    if (action == "seal") {
        if (args.output.empty()) throw std::invalid_argument("datapack seal needs --output FOLDER");
        if (!std::filesystem::is_regular_file(ledger())) (void)wf::executeJob(rootJob(target, args), root);
        print(wf::sealDataPack(target, ledger(), args.output));
        return 0;
    }
    if (action == "export") {
        if (args.output.empty()) throw std::invalid_argument("datapack export needs --output FOLDER");
        print(wf::exportDataPack(target, args.output));
        return 0;
    }
    throw std::invalid_argument("unknown datapack action: " + action);
}

}  // namespace dynsyslab::cli
