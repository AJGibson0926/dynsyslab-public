#include "dynsyslab_viewer/NativePayload.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    const auto path = std::filesystem::temp_directory_path()
        / "dynsyslab-koopman-validation-payload.npz";
    try {
        dynsyslab::core::NpzRecord archive;
        const std::vector<double> first{
            10.0, 0.01, 0.03, 0.02, 0.01, 0.04,
            10.0, 0.02, 0.05, 0.03, 0.02, 0.06};
        const std::vector<double> second{
            10.0, 0.04, 0.08, 0.05, 0.03, 0.09,
            10.0, 0.05, 0.10, 0.06, 0.04, 0.11};
        archive.arrays["one_step_metrics"] =
            dynsyslab::core::NpyArray::fromDoubles({2, 6}, first);
        archive.arrays["anchored_long_term_metrics"] =
            dynsyslab::core::NpyArray::fromDoubles({2, 6}, second);
        archive.arrays["generator_metrics"] =
            dynsyslab::core::NpyArray::fromDoubles({2, 6}, second);
        archive.arrays["horizons"] =
            dynsyslab::core::NpyArray::fromDoubles({3}, {0.1, 0.5, 1.0});
        archive.arrays["horizon_relative_l2"] =
            dynsyslab::core::NpyArray::fromDoubles(
                {2, 3}, {0.01, 0.05, 0.2, 0.02, 0.08, 0.3});
        archive.arrays["metadata_json"] =
            dynsyslab::core::NpyArray::fromString(
                R"({"schema":"dynsyslab.derived/v1","label":"KROM validation","method":"koopman-validation","dataset_role":"test","outside_training_fraction":0.25})");
        dynsyslab::core::writeNpzAtomic(path, archive);

        const auto payload = dsl::buildNativePayload(path, "koopman_validation");
        require(payload.figures.size() == 4,
                "validation dashboard figure count is wrong");
        require(payload.figures[0].traces.size() == 3,
                "validation metric comparisons are missing");
        require(payload.figures[3].traces[0].heatmap
                    && payload.figures[3].traces[0].zShape
                        == std::vector<int>({2, 3}),
                "finite-horizon heatmap is malformed");
        require(payload.warnings.size() == 1,
                "extrapolation warning was not displayed");
        archive.arrays["generator_available"]=dynsyslab::core::NpyArray::fromDoubles({2},{1,0});
        dynsyslab::core::writeNpzAtomic(path,archive);
        const auto discrete=dsl::buildNativePayload(path,"koopman_validation");
        require(discrete.figures[0].traces[2].x==std::vector<double>{0.},
                "unavailable generator metric presented as zero error");
        require(discrete.warnings.size()==2,"missing generator availability warning");
        std::error_code ignored;
        std::filesystem::remove(path, ignored);
        std::cout << "dynsyslab_koopman_validation_payload_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored;
        std::filesystem::remove(path, ignored);
        std::cerr << "dynsyslab_koopman_validation_payload_tests: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
