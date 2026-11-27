#include "CliQueue.hpp"

#include "CliDataPack.hpp"
#include "CliUsage.hpp"

#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Benchmark.hpp>
#include <dynsyslab/workflow/Campaign.hpp>
#include <dynsyslab/workflow/Experiment.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>
#include <dynsyslab/workflow/JobTemplates.hpp>
#include <dynsyslab/workflow/Queue.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::cli {
namespace {

namespace json = dynsyslab::core::json;
using dynsyslab::core::narrow;

// Operands and --options after the command word(s); flags take no value.
struct Tokens {
    std::vector<std::string> operands;
    std::multimap<std::string, std::string> options;
    std::set<std::string> flags;

    Tokens(const Arguments& args, std::size_t first) {
        static const std::set<std::string> flagNames{"--hold", "--any-runtime"};
        for (std::size_t index = first; index < args.positional.size(); ++index) {
            const auto token = narrow(args.positional[index]);
            if (token.rfind("--", 0) != 0) { operands.push_back(token); continue; }
            if (flagNames.count(token)) { flags.insert(token); continue; }
            if (index + 1 >= args.positional.size())
                throw std::invalid_argument(token + " needs a value");
            options.emplace(token, narrow(args.positional[++index]));
        }
    }
    bool has(const std::string& name) const { return options.count(name) != 0; }
    std::string text(const std::string& name, std::string fallback = {}) const {
        const auto found = options.find(name);
        return found == options.end() ? fallback : found->second;
    }
    double number(const std::string& name, double fallback) const {
        return has(name) ? std::stod(text(name)) : fallback;
    }
    std::size_t count(const std::string& name, std::size_t fallback = 0) const {
        const double value = number(name, static_cast<double>(fallback));
        if (value < 0.0) throw std::invalid_argument(name + " must be nonnegative");
        return static_cast<std::size_t>(value);
    }
    std::vector<std::string> all(const std::string& name) const {
        std::vector<std::string> result;
        const auto range = options.equal_range(name);
        for (auto item = range.first; item != range.second; ++item) result.push_back(item->second);
        return result;
    }
};

std::string readText(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to open " + path.u8string());
    std::ostringstream out; out << stream.rdbuf(); return out.str();
}

void print(const std::string& schema, const std::string& operation, json::Value payload) {
    json::Value::Object value;
    value["schema"] = schema;
    value["operation"] = operation;
    value["payload"] = std::move(payload);
    std::cout << json::dump(json::Value(std::move(value))) << '\n';
}

void printQueue(const std::string& operation, json::Value payload) {
    print("dynsyslab.job-queue/v1", operation, std::move(payload));
}

void emit(const json::Value& value, const Arguments& args) {
    if (args.output.empty()) { std::cout << json::dump(value) << '\n'; return; }
    const auto path = dynsyslab::core::prepareOutputFile(args.output, "--output");
    { std::ofstream out(path, std::ios::binary | std::ios::trunc); out << json::dump(value) << '\n'; }
    std::cout << json::dump(json::Value::Object{{"ok", true}, {"path", path.generic_u8string()}}) << '\n';
}

dynsyslab::core::ExecutionLimits limitsFrom(const Tokens& tokens) {
    return {tokens.count("--max-threads"),
            static_cast<std::size_t>(tokens.number("--max-memory-mib", 0.0) * 1024.0 * 1024.0)};
}

// Submission flags are sugar for the job's `queue` block.
void applySubmitFlags(json::Value& job, const Tokens& tokens) {
    auto& queue = job.at("queue");
    if (tokens.has("--retry-limit")) queue.at("retry_limit") = static_cast<double>(tokens.count("--retry-limit"));
    if (tokens.has("--max-concurrent"))
        queue.at("max_concurrent") = static_cast<double>(tokens.count("--max-concurrent"));
    if (tokens.has("--reuse")) queue.at("reuse") = tokens.text("--reuse");
    if (tokens.flags.count("--hold")) queue.at("hold") = true;
    for (const auto& item : tokens.all("--after")) {
        const auto colon = item.find(':');
        json::Value::Object dependency{{"job", item.substr(0, colon)}};
        if (colon != std::string::npos) dependency["condition"] = item.substr(colon + 1);
        queue.at("depends_on").array().emplace_back(std::move(dependency));
    }
}

json::Value runJson(const dynsyslab::workflow::QueueRunResult& result) {
    using State = dynsyslab::workflow::QueueRunResult::State;
    json::Value::Object value;
    if (result.state == State::empty) { value["status"] = "empty"; return json::Value(std::move(value)); }
    value["job"] = dynsyslab::workflow::queueJobJson(*result.job);
    if (result.run) value["result"] = dynsyslab::workflow::queueRunJson(*result.run);
    value["status"] = result.state == State::completed ? "completed"
        : result.state == State::cancelled ? "cancelled"
        : result.state == State::requeued ? "requeued" : "paused";
    if (!result.error.empty()) value["error"] = result.error;
    return json::Value(std::move(value));
}

std::string requiredId(const Tokens& tokens, const char* operation) {
    if (tokens.operands.empty())
        throw std::invalid_argument(std::string("queue ") + operation + " needs a job id");
    return tokens.operands.front();
}

}  // namespace

dynsyslab::core::ExecutionLimits executionLimits(const Arguments& args) {
    return limitsFrom(Tokens(args, 2));
}

int queueCommand(const Arguments& args) {
    if (args.positional.size() < 2) { usage(); return 2; }
    dynsyslab::workflow::JobQueue jobs(args.archiveRoot);
    const auto operation = narrow(args.positional[1]);
    const Tokens tokens(args, 2);
    if (operation == "submit") {
        if (tokens.operands.empty()) throw std::invalid_argument("queue submit needs a job file");
        const std::filesystem::path file = dynsyslab::core::wide(tokens.operands.front());
        auto job = dynsyslab::workflow::applyJobInputs(dynsyslab::workflow::bindJobFile(
            json::parse(readText(file)), file), args.params, args.sets);
        applyJobGroupRuntime(job, args);
        applySubmitFlags(job, tokens);
        dynsyslab::workflow::validateJob(job);
        printQueue("submit", jobs.submitDetailed(job).toJson());
        return 0;
    }
    if (operation == "list" || operation == "tray") {
        json::Value::Array values;
        for (const auto& job : operation == "list" ? jobs.list() : jobs.tray())
            values.push_back(dynsyslab::workflow::queueJobJson(job));
        printQueue(operation, json::Value(std::move(values)));
        return 0;
    }
    if (operation == "status") {
        const auto job = jobs.find(requiredId(tokens, "status"));
        if (!job) return 3;
        printQueue("status", dynsyslab::workflow::queueJobJson(*job));
        return 0;
    }
    if (operation == "cancel" || operation == "pause" || operation == "resume"
        || operation == "release" || operation == "remove") {
        const auto id = requiredId(tokens, operation.c_str());
        if (!jobs.find(id)) return 3;
        if (operation == "remove") {
            jobs.remove(id);
            printQueue("remove", json::Value::Object{{"removed", id}});
            return 0;
        }
        const auto record = operation == "cancel" ? jobs.cancel(id)
            : operation == "pause" ? jobs.pause(id)
            : operation == "resume" ? jobs.resume(id) : jobs.release(id);
        printQueue(operation, dynsyslab::workflow::queueJobJson(record));
        return 0;
    }
    if (operation == "resources") {
        printQueue("resources", dynsyslab::workflow::queueResources(args.archiveRoot));
        return 0;
    }
    if (operation == "watchdog") {
        json::Value::Array findings;
        for (const auto& item : jobs.watchdog(tokens.number("--grace-seconds", 120.0)))
            findings.emplace_back(json::Value::Object{{"job", item.job},
                {"classification", item.classification}, {"worker_alive", item.workerAlive},
                {"heartbeat_age_seconds", item.heartbeatAgeSeconds},
                {"quarantined", item.quarantined}});
        printQueue("watchdog", json::Value(std::move(findings)));
        return 0;
    }
    if (operation == "recover") {
        if (tokens.operands.empty()) throw std::invalid_argument("queue recover needs age seconds");
        printQueue("recover", json::Value::Object{{"recovered_jobs",
            static_cast<double>(jobs.recoverStale(std::stod(tokens.operands.front())))}});
        return 0;
    }
    if (operation == "run-next") {
        try {
            printQueue("run-next", runJson(jobs.runNext("dynsyslab-native-cli", {}, limitsFrom(tokens))));
        } catch (const std::exception& error) {
            printQueue("run-next", json::Value::Object{{"status", "failed"}, {"error", error.what()}});
            return 1;
        }
        return 0;
    }
    if (operation == "serve") {
        dynsyslab::workflow::QueueServeOptions options;
        options.workers = tokens.count("--workers", 1);
        options.admission = {tokens.count("--max-running"), tokens.count("--pool-threads"),
            static_cast<std::size_t>(tokens.number("--pool-memory-mib", 0.0) * 1024.0 * 1024.0)};
        options.limits = limitsFrom(tokens);
        options.pollSeconds = tokens.number("--poll-seconds", 2.0);
        options.maxJobs = tokens.count("--max-jobs");
        options.idleExitSeconds = tokens.number("--idle-exit", -1.0);
        options.watchdogGraceSeconds = tokens.number("--grace-seconds", 120.0);
        options.commandTemplate = tokens.text("--command-template");
        printQueue("serve", jobs.serve(options).toJson());
        return 0;
    }
    throw std::invalid_argument("unknown queue operation");
}

int campaignCommand(const Arguments& args) {
    if (args.positional.size() < 2)
        throw std::invalid_argument("campaign requires list, status, validate, repair, restructure, review or profiles");
    const auto operation = narrow(args.positional[1]);
    const Tokens tokens(args, 2);
    if (operation == "list") { std::cout << json::dump(dynsyslab::workflow::campaignList(args.archiveRoot)) << '\n'; return 0; }
    if (operation == "profiles") { std::cout << json::dump(dynsyslab::workflow::campaignProfiles()) << '\n'; return 0; }
    if (tokens.operands.empty()) throw std::invalid_argument("campaign " + operation + " needs a campaign id");
    if (operation == "status") {
        std::cout << json::dump(dynsyslab::workflow::campaignStatus(args.archiveRoot, tokens.operands.front())) << '\n';
        return 0;
    }
    if (operation == "repair") {
        std::cout << json::dump(dynsyslab::workflow::campaignRepair(
            args.archiveRoot, tokens.operands.front())) << std::endl;
        return 0;
    }
    if (operation == "restructure") {
        std::cout << json::dump(dynsyslab::workflow::campaignRestructure(args.archiveRoot,
            tokens.operands.front(), tokens.count("--points-per-chunk", 0), tokens.text("--reason"))) << '\n';
        return 0;
    }
    if (operation == "review") {
        const auto report = dynsyslab::workflow::campaignReview(args.archiveRoot, tokens.operands.front());
        std::cout << json::dump(report) << '\n';
        return report.at("safe_to_execute").boolean() ? 0 : 4;
    }
    if (operation == "validate") {
        const auto report = dynsyslab::workflow::campaignValidation(args.archiveRoot, tokens.operands.front());
        std::cout << json::dump(report) << '\n';
        return report.at("valid").boolean() ? 0 : 4;
    }
    throw std::invalid_argument("campaign requires list, status, validate, repair, restructure, review or profiles");
}

int benchmarkCommand(const Arguments& args) {
    const auto operation = args.positional.size() > 1 ? narrow(args.positional[1]) : "run";
    const bool explicitRun = operation == "run" || operation == "gate" || operation == "cases";
    const Tokens tokens(args, explicitRun ? 2 : 1);
    if (operation == "cases") {
        json::Value::Array names;
        for (const auto& name : dynsyslab::workflow::benchmarkCaseNames()) names.emplace_back(name);
        std::cout << json::dump(json::Value(std::move(names))) << '\n';
        return 0;
    }
    if (operation == "gate") {
        if (tokens.operands.size() != 2)
            throw std::invalid_argument("benchmark gate needs a reference and a candidate report");
        json::Value::Object required;
        if (tokens.has("--required-speedup")) required["*"] = tokens.number("--required-speedup", 1.0);
        for (const auto& item : tokens.all("--require")) {
            const auto equals = item.find('=');
            if (equals == std::string::npos) throw std::invalid_argument("--require expects CASE=SPEEDUP");
            required[item.substr(0, equals)] = std::stod(item.substr(equals + 1));
        }
        if (required.empty()) required["*"] = 1.0;
        const auto report = dynsyslab::workflow::benchmarkGate(
            json::load(dynsyslab::core::wide(tokens.operands[0])),
            json::load(dynsyslab::core::wide(tokens.operands[1])), json::Value(std::move(required)),
            !tokens.flags.count("--any-runtime"));
        emit(report, args);
        return report.at("passed").boolean() ? 0 : 4;
    }
    if (!explicitRun && operation.rfind("--", 0) != 0)
        throw std::invalid_argument("benchmark requires run, gate or cases");
    dynsyslab::workflow::BenchmarkOptions options;
    options.iterations = tokens.count("--iterations", 5);
    options.warmup = tokens.count("--warmup", 1);
    options.cases = tokens.all("--case");
    emit(dynsyslab::workflow::runBenchmarks(options), args);
    return 0;
}

int templateCommand(const Arguments& args) {
    if (args.positional.size() != 2)
        throw std::invalid_argument("template requires solve, postprocess or estimate");
    emit(dynsyslab::workflow::jobTemplate(narrow(args.positional[1])), args);
    return 0;
}

}  // namespace dynsyslab::cli
