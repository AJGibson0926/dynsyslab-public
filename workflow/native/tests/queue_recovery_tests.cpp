// Queue recovery: history policy, the deadline policy (checkpoint, requeue, continue), pause
// at a checkpoint, extension of a shorter run, verified input staging, campaign restructure
// and coverage review.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/workflow/Queue.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

namespace json = dynsyslab::workflow::json;
using dynsyslab::core::JobStatus;
using State = dynsyslab::workflow::QueueRunResult::State;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value solve(const std::string& label, double stop) {
    auto job = json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":1,"output_step":0.1,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"values":[[-1,0.1]]},
      "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[-1,0]}},
      "solver":{"method":"RK4"},"archive":{"enabled":true,"data_role":"free_training"}
    })");
    job.at("label") = label;
    job.at("time").at("stop") = stop;
    return job;
}

std::vector<double> values(const std::string& archive) {
    return dynsyslab::core::readNpz(std::filesystem::u8path(archive)).at("values").doubles();
}

json::Value resultOf(const dynsyslab::workflow::JobQueue& queue, const std::string& id) {
    return json::parse(queue.find(id)->message);
}

void historyPolicy(const std::filesystem::path& root) {
    dynsyslab::workflow::JobQueue queue(root);
    auto job = solve("policy", 1.0);
    job.at("history_policy") = json::parse(R"({"enabled":true,"estimated_seconds":400,
        "if_estimated_seconds_over":300,"checkpoint_interval_seconds":5,"pause_at_checkpoint":true})");
    auto history = json::parse(queue.find(queue.submit(job))->manifestJson).at("history");
    require(history.at("checkpoint_wall_seconds").number() == 5.0 && history.at("pause_at_checkpoint").boolean()
            && history.at("policy_resolution").at("reason").string() == "policy_matched", "policy applied");
    job.at("label") = "policy-short";
    job.at("history_policy").at("estimated_seconds") = 10.0;
    history = json::parse(queue.find(queue.submit(job))->manifestJson).at("history");
    require(history.at("checkpoint_wall_seconds").isNull()
            && history.at("policy_resolution").at("reason").string() == "below_estimate_threshold", "policy skipped");
    // campaign.checkpoint_profile "recoverable" goes through the same policy.
    auto campaign = solve("recoverable", 0.2);
    campaign.at("campaign") = json::parse(R"({"enabled":true,"checkpoint_profile":"recoverable",
        "axes":[{"path":"time.stop","values":[0.2]}]})");
    campaign.at("queue") = json::parse(R"({"estimated_seconds":1000})");
    const auto ids = queue.submitDetailed(campaign).jobIds;
    history = json::parse(queue.find(ids.front())->manifestJson).at("history");
    require(history.at("checkpoint_wall_seconds").number() == 600.0, "recoverable profile checkpoints");
    std::cout << "history policy\n";
}

void deadline(const std::filesystem::path& root, const std::vector<double>& reference) {
    dynsyslab::workflow::JobQueue queue(root);
    auto job = solve("deadline", 1.0);
    // The checkpoint deadline falls 0.1 ms after the start: every attempt checkpoints at its
    // first sample, pauses and is requeued; the next attempt continues from that checkpoint.
    job.at("queue") = json::parse(R"({"wall_time_seconds":60,"checkpoint_before_deadline_seconds":59.9999,
        "max_requeues":40})");
    const auto id = queue.submit(job);
    std::size_t requeues = 0;
    for (int attempt = 0; attempt < 40; ++attempt) {
        const auto run = queue.runNext("deadline-test");
        if (run.state == State::completed) break;
        require(run.state == State::requeued, "deadline checkpoint requeues the job");
        ++requeues;
    }
    require(queue.find(id)->status == JobStatus::completed && requeues > 0, "deadline job completed after requeues");
    const auto result = resultOf(queue, id);
    const auto& archives = result.at("archives").array();
    require(result.at("continuation").at("rounds").number() == static_cast<double>(requeues)
            && archives.size() == requeues + 1, "one archive per round");
    const auto last = values(archives.back().string());
    require(std::abs(last[last.size() - 2] - reference[reference.size() - 2]) < 1.0e-13
            && std::abs(last.back() - reference.back()) < 1.0e-13, "continued solve reproduces the uninterrupted one");
    std::cout << "deadline: " << requeues << " requeues\n";
}

void pauseAtCheckpoint(const std::filesystem::path& root, const std::vector<double>& reference) {
    dynsyslab::workflow::JobQueue queue(root);
    auto job = solve("pause-checkpoint", 1.0);
    job.at("history") = json::parse(R"({"pause_at_checkpoint":true,"checkpoint_interval":0.3})");
    const auto id = queue.submit(job);
    // A pause file present at the first checkpoint stands in for `queue pause` mid-solve.
    auto paused = job;
    paused.at("history").at("pause_file") = (root / "user-pause").generic_string();
    { std::ofstream touch(root / "user-pause"); }
    const auto pausedId = queue.submit(paused);
    (void)queue.runNext("pause-test");  // the first job runs through
    const auto first = queue.runNext("pause-test");
    require(first.state == State::paused && queue.find(pausedId)->status == JobStatus::paused, "paused at a checkpoint");
    std::filesystem::remove(root / "user-pause");
    (void)queue.resume(pausedId);
    require(queue.runNext("pause-test").state == State::completed, "resumed from the checkpoint");
    const auto archives = resultOf(queue, pausedId).at("archives").array();
    require(archives.size() == 2, "paused prefix and continuation");
    const auto last = values(archives.back().string());
    require(std::abs(last.back() - reference.back()) < 1.0e-13, "pause and resume reproduce the solve");
    require(queue.find(id)->status == JobStatus::completed, "unpaused job completed");
    std::cout << "pause at checkpoint\n";
}

void extension(const std::filesystem::path& root, const std::vector<double>& reference) {
    dynsyslab::workflow::JobQueue queue(root);
    const auto shortId = queue.submit(solve("extend-me", 0.5));
    require(queue.runNext("extend-test").state == State::completed, "short run");
    auto longer = solve("extend-me", 1.0);
    longer.at("queue") = json::parse(R"({"reuse":"extend"})");
    const auto submitted = queue.submitDetailed(longer);
    require(submitted.extendedFrom == shortId, "extension found the shorter run");
    const auto manifest = json::parse(queue.find(submitted.jobIds.front())->manifestJson);
    require(manifest.at("initial_conditions").at("mode").string() == "continue_from"
            && manifest.at("time").at("start").number() == 0.5, "continuation job");
    require(queue.runNext("extend-test").state == State::completed, "extension ran");
    const auto archive = resultOf(queue, submitted.jobIds.front()).at("archives").array().front().string();
    const auto shorter = values(resultOf(queue, shortId).at("archives").array().front().string());
    const auto extended = values(archive);
    require(extended[0] == shorter[shorter.size() - 2] && extended[1] == shorter.back(), "starts at the shorter run's end");
    require(std::abs(extended.back() - reference.back()) < 1.0e-13, "extension reproduces the full run");
    std::cout << "extension\n";
}

void staging(const std::filesystem::path& root) {
    dynsyslab::workflow::JobQueue queue(root);
    const auto input = root / "input.txt";
    std::filesystem::create_directories(root);
    { std::ofstream out(input); out << "staged input"; }
    auto job = solve("staging", 0.2);
    job.at("queue") = json::Value::Object{{"staging", json::Value::Object{
        {"inputs", json::Value::Array{json::Value(input.generic_string())}}}}};
    const auto id = queue.submit(job);
    require(queue.runNext("staging-test").state == State::completed, "staged job ran");
    const auto manifest = resultOf(queue, id).at("staging");
    require(manifest.at("inputs").array().front().at("sha256").string() == dynsyslab::core::sha256File(input),
            "staging digest recorded");
    require(!std::filesystem::exists(std::filesystem::u8path(manifest.at("staging_directory").string())),
            "scratch cleaned");
    std::cout << "staging\n";
}

void restructure(const std::filesystem::path& root) {
    dynsyslab::workflow::JobQueue queue(root);
    auto job = solve("restructure-campaign", 0.2);
    job.at("campaign") = json::parse(R"({"enabled":true,"systems_per_job":4,
        "axes":[{"path":"time.stop","values":[0.1,0.2,0.3,0.4]}]})");
    const auto submitted = queue.submitDetailed(job);
    require(submitted.jobIds.size() == 1, "one chunk");
    const auto report = dynsyslab::workflow::campaignRestructure(root, submitted.campaignId, 2, "test");
    const auto& changes = report.at("restructured").array();
    require(changes.size() == 1 && changes.front().at("children").array().size() == 2
            && queue.find(submitted.jobIds.front())->status == JobStatus::cancelled, "chunk split");
    const auto review = dynsyslab::workflow::campaignReview(root, submitted.campaignId);
    require(review.at("safe_to_execute").boolean() && review.at("pending_chunks").number() == 2.0
            && review.at("needs_repair").array().empty(), "review after restructure");
    while (queue.runNext("restructure-test").state != State::empty) {}
    const auto validation = dynsyslab::workflow::campaignValidation(root, submitted.campaignId);
    require(validation.at("valid").boolean(), "restructured campaign validates");
    std::cout << "restructure and review\n";
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-queue-recovery-test";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    try {
        dynsyslab::workflow::JobQueue queue(root / "reference");
        const auto id = queue.submit(solve("reference", 1.0));
        (void)queue.runNext("reference");
        const auto reference = values(json::parse(queue.find(id)->message).at("archives").array().front().string());
        historyPolicy(root / "policy");
        deadline(root / "deadline", reference);
        pauseAtCheckpoint(root / "pause", reference);
        extension(root / "extension", reference);
        staging(root / "staging");
        restructure(root / "restructure");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::filesystem::remove_all(root, ignored);
    return 0;
}
