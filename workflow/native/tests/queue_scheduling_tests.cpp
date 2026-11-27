// Queue retry, campaign arrays and registry validation, cycle time, serve, templates,
// benchmark gates and the post-processing memory preflight.
#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Benchmark.hpp>
#include <dynsyslab/workflow/Campaign.hpp>
#include <dynsyslab/workflow/JobTemplates.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Queue.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

namespace json = dynsyslab::workflow::json;
using dynsyslab::core::JobStatus;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value solve() {
    return json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"scheduling-test",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":0.2,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,
        "values":[[-1,0.1]],"ranges":[]},
      "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[-1,0]}},
      "solver":{"method":"RK4"},
      "archive":{"enabled":true,"data_role":"free_training"}
    })");
}

void retry(const std::filesystem::path& root) {
    dynsyslab::workflow::JobQueue queue(root);
    dynsyslab::core::JobStore store(queue.database());
    dynsyslab::core::JobSubmission broken;
    broken.label = "broken"; broken.kind = "solve"; broken.retryLimit = 1;
    broken.manifestJson = R"({"schema":"dynsyslab.compute-job/v1","kind":"no-such-kind"})";
    const auto id = store.submit(broken);
    const auto first = queue.runNext("retry-test");
    require(first.state == dynsyslab::workflow::QueueRunResult::State::requeued, "failure was not requeued");
    bool threw = false;
    try { (void)queue.runNext("retry-test"); } catch (const std::exception&) { threw = true; }
    require(threw && store.find(id)->status == JobStatus::failed && store.find(id)->attempts == 2,
            "retry limit was not enforced");
}

void campaignArray(const std::filesystem::path& root) {
    auto job = solve();
    job.at("label") = "array-campaign";
    job.at("campaign") = json::parse(R"({"enabled":true,"mode":"cartesian","systems_per_job":2,
        "axes":[{"path":"time.stop","start":0.1,"stop":0.3,"points":3}]})");
    job.at("queue") = json::parse(R"({"max_concurrent":1,"retry_limit":1})");
    const auto plan = dynsyslab::workflow::planJob(job);
    require(plan.at("queue_tasks").number() == 2.0, "plan did not report two queue tasks");
    dynsyslab::workflow::JobQueue queue(root);
    const auto submitted = queue.submitDetailed(job);
    require(submitted.jobIds.size() == 2 && !submitted.campaignId.empty(), "array not chunked");
    const auto task = queue.find(submitted.jobIds.front());
    require(task->arrayMaxConcurrent == 1 && task->retryLimit == 1, "array options lost");
    const auto pending = dynsyslab::workflow::campaignValidation(root, submitted.campaignId);
    require(!pending.at("valid").boolean() && pending.at("missing_points").array().size() == 3,
            "validation passed before the campaign ran");
    dynsyslab::workflow::QueueServeOptions options;
    options.idleExitSeconds = 0.0; options.pollSeconds = 0.0;
    const auto summary = queue.serve(options);
    require(summary.completed == 2, "serve did not run both chunks");
    const auto status = dynsyslab::workflow::campaignStatus(root, submitted.campaignId);
    require(status.at("percent").number() == 100.0, "campaign status incomplete");
    const auto report = dynsyslab::workflow::campaignValidation(root, submitted.campaignId);
    require(report.at("valid").boolean(), "campaign validation failed: " + json::dump(report));
    const auto repair = dynsyslab::workflow::campaignRepair(root, submitted.campaignId);
    require(repair.at("resubmitted").array().empty(), "a complete campaign needs no repair");
    require(dynsyslab::workflow::queueResources(root).at("recommended_local_workers").number() >= 1.0,
            "resource report has no worker recommendation");
}

void reuseAndExtension(const std::filesystem::path& root) {
    dynsyslab::workflow::JobQueue queue(root);
    auto job = solve();
    job.at("label") = "reuse";
    const auto id = queue.submit(job);
    require(queue.runNext("reuse").state == dynsyslab::workflow::QueueRunResult::State::completed,
            "reuse source did not run");
    job.at("queue").at("reuse") = "identical";
    job.at("label") = "renamed";  // labels are not scientific content
    require(queue.submitDetailed(job).reusedJob == id, "identical spec was not reused");
    job.at("time").at("stop") = 0.4;
    const auto longer = queue.submitDetailed(job);
    require(longer.reusedJob.empty() && longer.extensionCandidates.array().size() == 1
                && longer.extensionCandidates.array().front().at("relation").string() == "shorter",
            "shorter compatible run was not reported");
    (void)queue.cancel(longer.jobIds.front());
}

void wallTime(const std::filesystem::path& root) {
    dynsyslab::workflow::JobQueue queue(root);
    auto job = solve();
    job.at("label") = "wall-time";
    job.at("time").at("stop") = 1000.0;
    job.at("time").at("output_step") = 0.001;
    job.at("queue").at("wall_time_seconds") = 0.2;
    const auto id = queue.submit(job);
    bool threw = false;
    try { (void)queue.runNext("wall-time"); } catch (const std::exception&) { threw = true; }
    require(threw && queue.find(id)->status == JobStatus::failed
                && queue.find(id)->message.find("wall time") != std::string::npos,
            "the wall-time limit did not stop the job");
}

void cycleTime() {
    auto job = solve();
    job.at("controller") = json::parse(R"({"kind":"waveform","settings":{"frequency_hz":2.0,"amplitude":0.1},
        "target":{"kind":"constant","state":[-1,0]}})");
    job.at("time") = json::parse(R"({"transient_cycles":3,"retained_cycles":2,"samples_per_cycle":10})");
    const auto plan = dynsyslab::workflow::planJob(job);
    require(plan.at("output_samples_per_trajectory").number() == 51.0, "cycle time grid is wrong");
    require(std::abs(plan.at("retained_start").number() - 1.5) < 1e-12, "retained start is wrong");
    job.at("campaign") = json::parse(R"({"enabled":true,"profile":"low",
        "axes":[{"path":"driver.frequency_hz","start":1,"stop":3}]})");
    const auto swept = dynsyslab::workflow::planJob(job);
    require(swept.at("system_count").number() == 31.0 && swept.at("systems_per_job").number() == 10.0,
            "quality profile defaults were not applied");
}

void templatesAndBenchmarks() {
    for (const auto& name : dynsyslab::workflow::jobTemplateNames())
        dynsyslab::workflow::validateJob(dynsyslab::workflow::jobTemplate(name));
    (void)dynsyslab::workflow::planJob(dynsyslab::workflow::jobTemplate("solve"));
    dynsyslab::workflow::BenchmarkOptions options;
    options.iterations = 2; options.warmup = 0; options.cases = {"ode_rk4_duffing"};
    const auto report = dynsyslab::workflow::runBenchmarks(options);
    require(report.at("samples").array().size() == 1, "benchmark case missing");
    const auto gate = dynsyslab::workflow::benchmarkGate(report, report,
        json::Value::Object{{"*", 0.5}});
    require(gate.at("passed").boolean(), "a report must pass a 0.5x gate against itself");
    const auto strict = dynsyslab::workflow::benchmarkGate(report, report,
        json::Value::Object{{"ode_rk4_duffing", 2.0}});
    require(!strict.at("passed").boolean(), "a report cannot be 2x faster than itself");
}

void memoryPreflight(const std::filesystem::path& root) {
    auto source = solve();
    source.at("time").at("stop") = 5.0;
    const auto run = dynsyslab::workflow::executeJob(source, root);
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess",
        "label":"preflight","sources":[],"methods":[{"name":"hdmd","enabled":true,
        "settings":{"observable":0,"matrix_size":100,"sampling_window":300}}],
        "archive":{"enabled":true}})");
    job.at("sources").array().emplace_back(run.archives.front().generic_u8string());
    const auto plan = dynsyslab::workflow::planJob(job);
    const auto& geometry = plan.at("geometry");
    require(geometry.at("sources").array().front().at("samples").number() == 501.0, "source not inspected");
    require(geometry.at("estimated_peak_bytes").number() > 0.0, "no working-set estimate");
    job.at("execution") = json::parse(R"({"max_memory_mib":0.01})");
    bool refused = false;
    try { (void)dynsyslab::workflow::executeJob(job, root); }
    catch (const dynsyslab::core::ExecutionLimitExceeded&) { refused = true; }
    require(refused, "the memory preflight did not refuse the job");
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-queue-scheduling-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    try {
        retry(root); campaignArray(root); reuseAndExtension(root); wallTime(root); cycleTime();
        templatesAndBenchmarks(); memoryPreflight(root);
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_queue_scheduling_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_queue_scheduling_tests: FAIL: " << error.what() << '\n';
        std::filesystem::remove_all(root, ignored);
        return 1;
    }
}
