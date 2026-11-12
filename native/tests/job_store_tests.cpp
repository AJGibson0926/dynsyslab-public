#include <dynsyslab/core/JobStore.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-native-job-store-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    try {
        dynsyslab::core::JobStore store(root / "state.sqlite");
        const auto id = store.submit("test", "solve", "{\"kind\":\"solve\"}");
        require(store.find(id)->status == dynsyslab::core::JobStatus::queued, "queued job");
        const auto claimed = store.claimNext("test-worker");
        require(claimed && claimed->id == id && claimed->attempts == 1, "atomic claim");
        store.updateProgress(id, 0.5, "halfway");
        require(store.find(id)->progress == 0.5, "progress update");
        store.complete(id, "done");
        require(store.find(id)->status == dynsyslab::core::JobStatus::completed, "completed job");
        require(!store.claimNext("idle-worker"), "empty queue");
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_job_store_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_job_store_tests: FAIL: " << error.what() << '\n';
        std::filesystem::remove_all(root, ignored); return 1;
    }
}
