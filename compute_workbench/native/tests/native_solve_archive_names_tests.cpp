// Solve archive publication: repeated runs in the same second never share (or replace) an archive,
// manifest or batch identity, and archive roots deeper than MAX_PATH publish complete results.
#include <dynsyslab/workflow/Solve.hpp>

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LongPath.hpp>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace {
namespace fs = std::filesystem;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }

const char* kJob = R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"same-second",
  "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
  "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,"values":[[0.1,0]],"ranges":[]},
  "solver":{"method":"RK4","step_size":0.005},"archive":{"enabled":true}})";
}  // namespace

int main() {
    const auto root = fs::temp_directory_path() / ("dynsyslab-solve-names-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        const auto job = dynsyslab::core::json::parse(kJob);
        std::set<fs::path> archives, manifests;
        for (int run = 0; run < 3; ++run) {  // well inside one second
            const auto result = dynsyslab::workflow::executeSolve(job, root / "archives", root / "progress.json");
            require(result.archives.size() == 1, "one archive per run");
            archives.insert(result.archives.front());
            manifests.insert(result.hierarchyManifest);
        }
        require(archives.size() == 3 && manifests.size() == 3, "same-second runs get distinct archives and manifests");
        for (const auto& path : archives) require(fs::is_regular_file(path), "every archive kept");

        auto deep = root / "archives";
        while (deep.wstring().size() < 300) deep /= "a-long-directory-name-for-max-path";
        const auto result = dynsyslab::workflow::executeSolve(job, deep, root / "progress.json");
        require(result.archives.size() == 1, "deep root archive");
        require(result.hierarchyManifest.wstring().size() > 260, "manifest path is past MAX_PATH");
        require(fs::is_regular_file(dynsyslab::core::extendedLengthPath(result.archives.front()))
                && fs::is_regular_file(dynsyslab::core::extendedLengthPath(result.hierarchyManifest)),
                "deep root publishes archive and manifest");

        std::error_code ignored; fs::remove_all(dynsyslab::core::extendedLengthPath(root), ignored);
        std::cout << "dynsyslab_native_solve_archive_names_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored; fs::remove_all(dynsyslab::core::extendedLengthPath(root), ignored);
        std::cerr << "dynsyslab_native_solve_archive_names_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
