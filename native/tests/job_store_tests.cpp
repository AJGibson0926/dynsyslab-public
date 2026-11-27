#include <dynsyslab/core/CampaignRegistry.hpp>
#include <dynsyslab/core/JobStore.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
using dynsyslab::core::DependencyCondition;
using dynsyslab::core::JobStatus;
using dynsyslab::core::JobStore;
using dynsyslab::core::JobSubmission;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

JobSubmission job(const std::string& label) {
    JobSubmission result;
    result.label = label; result.kind = "solve"; result.manifestJson = "{\"label\":\"" + label + "\"}";
    return result;
}

void basics(JobStore& store) {
    const auto id = store.submit("test", "solve", "{\"kind\":\"solve\"}");
    require(store.find(id)->status == JobStatus::queued, "queued job");
    const auto claimed = store.claimNext("test-worker");
    require(claimed && claimed->id == id && claimed->attempts == 1, "atomic claim");
    store.updateProgress(id, 0.5, "halfway");
    require(store.find(id)->progress == 0.5, "progress update");
    store.complete(id, "done");
    require(store.find(id)->status == JobStatus::completed, "completed job");
    require(!store.claimNext("idle-worker"), "empty queue");
}

void retries(JobStore& store) {
    auto item = job("retry"); item.retryLimit = 1;
    const auto id = store.submit(item);
    require(store.claimNext("w")->id == id, "retry claim 1");
    require(store.fail(id, "boom"), "first failure is requeued");
    require(store.find(id)->status == JobStatus::queued, "requeued status");
    require(store.claimNext("w")->attempts == 2, "second attempt");
    require(!store.fail(id, "boom"), "second failure is terminal");
    require(store.find(id)->status == JobStatus::failed, "terminal failure");
}

void dependencies(JobStore& store) {
    const auto first = store.submit(job("first"));
    auto afterOk = job("after-ok"); afterOk.dependencies = {{first, DependencyCondition::afterok}};
    auto afterNotOk = job("after-not-ok");
    afterNotOk.dependencies = {{first, DependencyCondition::afternotok}};
    auto afterAny = job("after-any"); afterAny.dependencies = {{first, DependencyCondition::afterany}};
    const auto ok = store.submit(afterOk), notOk = store.submit(afterNotOk), any = store.submit(afterAny);
    require(store.claimNext("w")->id == first, "dependency source runs first");
    require(!store.claimNext("w"), "dependents wait for a running source");
    store.fail(first, "failed on purpose");
    const auto next = store.claimNext("w");
    require(next && next->id == notOk, "afternotok runs after a failure");
    require(store.find(ok)->status == JobStatus::cancelled, "afterok is cancelled after a failure");
    store.complete(notOk, "done");
    require(store.claimNext("w")->id == any, "afterany runs after any end");
    store.complete(any, "done");
    using dynsyslab::core::dependencyState;
    using dynsyslab::core::DependencyState;
    require(dependencyState(DependencyCondition::afterok, std::nullopt) == DependencyState::impossible,
            "unknown dependency is impossible");
    require(dependencyState(DependencyCondition::afternotok, JobStatus::completed)
                == DependencyState::impossible, "afternotok after success is impossible");
}

void arrays(JobStore& store) {
    std::vector<std::string> source, sink;
    for (std::size_t index = 0; index < 3; ++index) {
        auto task = job("array-" + std::to_string(index));
        task.arrayId = "array-a"; task.arrayIndex = index; task.arrayMaxConcurrent = 1;
        source.push_back(store.submit(task));
    }
    for (std::size_t index = 0; index < 2; ++index) {
        auto task = job("corr-" + std::to_string(index));
        task.arrayId = "array-b"; task.arrayIndex = index;
        task.dependencies = {{"array-a", DependencyCondition::aftercorr}};
        sink.push_back(store.submit(task));
    }
    const auto a0 = store.claimNext("w");
    require(a0 && a0->id == source[0], "array task 0 first");
    require(!store.claimNext("w"), "array cap and aftercorr hold the rest");
    store.complete(source[0], "done");
    const auto a1 = store.claimNext("w");
    require(a1 && a1->id == source[1], "array task 1 after task 0 ends");
    const auto b0 = store.claimNext("w");
    require(b0 && b0->id == sink[0], "aftercorr task 0 follows array task 0");
    store.complete(sink[0], "done"); store.fail(source[1], "bad");
    require(store.claimNext("w")->id == source[2], "array task 2");
    store.complete(source[2], "done");
    require(!store.claimNext("w"), "nothing left to claim");
    require(store.find(sink[1])->status == JobStatus::cancelled, "aftercorr cancelled with its partner");
    require(store.arrayTasks("array-a").size() == 3, "array task listing");
}

void claimsAndAdmission(JobStore& store) {
    auto first = job("spec-a"); first.specKey = "same";
    auto second = job("spec-b"); second.specKey = "same";
    auto heavy = job("heavy"); heavy.resources.memoryBytes = 600;
    const auto a = store.submit(first), b = store.submit(second), h = store.submit(heavy);
    const dynsyslab::core::AdmissionPolicy policy{0, 0, 1000};
    auto big = job("big"); big.resources.memoryBytes = 500;
    require(store.claimNext("w", policy)->id == a, "spec owner claimed");
    const auto next = store.claimNext("w", policy);
    require(next && next->id == h, "same spec waits; the next job is admitted");
    store.complete(a, "done");
    require(store.claimNext("w", policy)->id == b, "spec free after its owner ends");
    const auto blocked = store.submit(big);
    require(!store.claimNext("w", policy), "memory pool blocks a job that does not fit");
    store.complete(h, "done"); store.complete(b, "done");
    require(store.claimNext("w", policy)->id == blocked, "job admitted once memory frees");
    store.complete(blocked, "done");
    require(store.findBySpec("same")->id == b, "newest completed job of a spec");
    require(dynsyslab::core::admissionDenial({4, 0}, {}, {0, 2, 0}).empty(),
            "a lone job is never blocked by the pool");
}

void pauseTrayWatchdog(JobStore& store) {
    auto held = job("tray"); held.reserved = true;
    const auto tray = store.submit(held);
    require(!store.claimNext("w"), "reserved jobs are not claimed");
    store.release(tray);
    require(store.pause(tray) == JobStatus::paused, "queued job pauses at once");
    require(!store.claimNext("w"), "paused jobs are not claimed");
    store.resume(tray);
    const auto running = store.claimNext("w", {}, 4000000000ull);  // a pid that is not alive
    require(running && running->id == tray, "resumed job claimed");
    require(store.pause(tray) == JobStatus::running && store.find(tray)->pauseRequested,
            "running job gets a pause request");
    const auto findings = store.watchPausedJobs(120.0);
    require(findings.size() == 1 && findings.front().quarantined
                && findings.front().classification == "exited_without_acknowledgement",
            "watchdog quarantines an exited worker");
    require(store.find(tray)->status == JobStatus::quarantined, "quarantined status");
    store.resume(tray);
    require(store.claimNext("w")->id == tray, "reviewed job runs again");
    store.pause(tray);
    store.acknowledgePause(tray, "paused at a safe boundary");
    require(store.find(tray)->status == JobStatus::paused, "acknowledged pause");
    auto removable = job("removable"); removable.reserved = true;
    const auto gone = store.submit(removable);
    store.remove(gone);
    require(!store.find(gone), "reserved job removed");
}

void registry(const std::filesystem::path& database) {
    dynsyslab::core::CampaignRegistry campaigns(database);
    campaigns.registerPlan({"c1", "campaign", "sweep", 4, "{\"first\":true}", {}});
    campaigns.registerPlan({"c1", "campaign", "sweep", 4, "{\"second\":true}", {}});
    require(campaigns.plan("c1")->planJson == "{\"first\":true}", "original plan is kept");
    campaigns.recordChunk({"c1", 0, 2, "job-a", {}, {}});
    campaigns.recordChunk({"c1", 0, 2, {}, "manifest.json", {}});
    const auto chunks = campaigns.chunks("c1");
    require(chunks.size() == 1 && chunks.front().jobId == "job-a"
                && chunks.front().manifestPath == "manifest.json", "chunk upsert keeps both fields");
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-native-job-store-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    try {
        {
            JobStore store(root / "state.sqlite");
            basics(store); retries(store); dependencies(store); arrays(store);
            claimsAndAdmission(store); pauseTrayWatchdog(store);
        }
        registry(root / "state.sqlite");
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_job_store_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_job_store_tests: FAIL: " << error.what() << '\n';
        std::filesystem::remove_all(root, ignored); return 1;
    }
}
