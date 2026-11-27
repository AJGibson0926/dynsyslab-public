// End-to-end RPE chain on the one krom object: free training solve -> identity coordinates ->
// H-DMD -> SINDy eigenfunctions -> krom_extract (one oscillatory pair) -> KLQR tracking a
// ramped multi-sine radius target. It runs the tracked experiment
// examples/native/experiments/rpe-viscous-klqr-multisine/base-r.json, shortened, and checks the
// frozen qualitative outcome: a published krom with an oscillatory pair and closed-loop
// tracking of the radius.
// H-DMD: Arbabi & Mezic 2017, SIAM J. Appl. Dyn. Syst. 16:2096-2126; SINDy eigenfunctions: Kaiser,
// Kutz & Brunton 2021, Mach. Learn.: Sci. Technol. 2:035023; KLQR: Kaiser, Kutz & Brunton 2021;
// Rayleigh-Plesset: Plesset & Prosperetti 1977, Annu. Rev. Fluid Mech. 9:145-185.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Experiment.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / ("dynsyslab-rpe-klqr-pipeline-" + std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::error_code ignored;
    try {
        std::ifstream stream(std::filesystem::path(DYNSYSLAB_RPE_EXPERIMENT), std::ios::binary);
        require(bool(stream), "tracked experiment is missing");
        std::ostringstream text; text << stream.rdbuf();
        auto job = dynsyslab::workflow::bindJobFile(json::parse(text.str()), DYNSYSLAB_RPE_EXPERIMENT);
        job = dynsyslab::workflow::applyJobInputs(std::move(job), {
            {"method_settings", R"({"hdmd":{"observable":0,"matrix_size":200,"sampling_window":3000,"sample_stride":1,
                "rank":0,"svd_tolerance":1e-10,"mode_sort":"terminal-amplitude-energy"},
                "sindy":{"gamma_method":"system_rhs","polynomial_degree":4,"threshold":0.05,"max_iterations":20,
                "normalize":true,"reuse_conjugate_fits":true,"compute_training_values":false,
                "archive_storage":"coefficient-first","reconstruction":"none"}})"}},
            {{"train.time.stop", "8"}, {"control.time.stop", "4"}});
        const auto result = dynsyslab::workflow::executeJob(job, root);
        std::filesystem::path control, krom;
        for (const auto& path : result.archives) {
            const auto meta = json::parse(core::readNpzIndex(path).metadataJson);
            if (meta.at("method").string() == "krom-extract") krom = path;
            if (!meta.at("controller").isNull() || path.u8string().find("control") != std::string::npos) control = path;
        }
        require(!krom.empty(), "the pipeline published no extracted krom");
        const auto model = core::loadKoopmanModel(krom);
        bool oscillatory = false;
        for (const auto& rate : model.continuousEigenvalues()) oscillatory = oscillatory || std::abs(rate.imag()) > 0.1;
        require(model.liftedDimension() == 2 && oscillatory, "the krom is not one oscillatory pair");
        require(!control.empty(), "the pipeline published no controlled run");
        const auto run = core::readNpz(control);
        const auto time = run.at("time").doubles();
        const auto values = run.at("values").doubles();
        const auto target = run.at("control_target_values").doubles();
        const std::size_t n = run.at("values").shape.back();
        double squares = 0; std::size_t count = 0;
        for (std::size_t k = 0; k < time.size(); ++k)
            if (time[k] >= 2.0) { squares += std::pow(values[k * n] - target[k * n], 2); ++count; }
        const double rmse = std::sqrt(squares / double(count));
        require(std::isfinite(rmse) && rmse < 0.25, "closed-loop radius tracking RMSE " + std::to_string(rmse));
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_rpe_klqr_pipeline_tests: PASS (tracking RMSE " << rmse << ")\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_rpe_klqr_pipeline_tests: FAIL: " << error.what() << '\n';
        std::filesystem::remove_all(root, ignored);
        return 1;
    }
}
