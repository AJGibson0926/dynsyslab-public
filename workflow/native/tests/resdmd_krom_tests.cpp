// ResDMD end to end: a ResDMD job (Colbrook & Townsend, Comm. Pure Appl. Math. 77, 2024) publishes a
// krom that loads through loadKoopmanModel like H-DMD's and SINDy's; observable_reconstruction fits a
// readout C on top and input_map fits an input map B from a forced record, each as a child krom that
// references its parent and stores only what it adds.
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
namespace w = dynsyslab::workflow;
namespace fs = std::filesystem;
namespace json = c::json;

void check(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

// A lightly damped linear oscillator x'' + 0.1 x' + x = u (Duffing x'' = -delta x' + alpha x - beta x^3 + u, alpha = -1, beta = 0): lambda = -0.05 +- i sqrt(0.9975).
constexpr const char* kSystem = R"({"name":"duffing","parameters":{"alpha":-1,"beta":0,"delta":0.1}})";

fs::path solve(const fs::path& root, bool forced) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"resdmd-source",
      "time":{"start":0,"stop":20,"output_step":0.05},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[1.0,0.0]]},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12},"archive":{"enabled":true}})");
    job.object()["system"] = json::parse(kSystem);
    if (forced) {
        job.at("solver").object()["control_evaluation"] = "sampled";
        job.at("solver").object()["control_sample_time"] = 0.05;
        job.object()["controller"] = json::parse(
            R"({"kind":"random_pulse","settings":{"model_amplitude":0.2,"model_pulse_duration":0.5,"seed":3}})");
    }
    return w::executeJob(job, root).archives.front();
}

w::PostprocessResult postprocess(const fs::path& source, const std::string& methods, const fs::path& root) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"resdmd-krom","methods":)"
                           + methods + "}");
    job.at("sources") = json::Value::Array{source.generic_u8string()};
    return w::executePostprocess(job, root, {});
}

bool isChildOf(const fs::path& child, const fs::path& parent) {
    const auto meta = json::parse(c::readNpz(child).metadataJson());
    const auto& block = meta.at("krom").at("parent");
    return !block.isNull() && fs::equivalent(fs::u8path(block.at("path").string()), parent);
}

}  // namespace

int main() {
    fs::path root;
    try {
        root = fs::temp_directory_path() / ("dynsyslab-resdmd-krom-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(root);
        const auto free = solve(root, false);
        const std::complex<double> expected(-0.05, std::sqrt(0.9975));

        // 1. ResDMD on a shared Hankel matrix (examples/native/spectral-koopman-job.json) publishes a krom.
        const auto spectral = postprocess(free, R"([{"name":"hankel_matrix","settings":{"observables":[0,1],
            "window":{"delays":10,"columns":200}}},{"name":"resdmd","settings":{"use_previous_hankel":true,"rank":2,
            "dictionary":{"kind":"linear","constant":false},"residual_grid":[[-0.05,1.0],[0.0,0.0]]}}])", root);
        const auto resdmd = spectral.archives.back();
        check(c::dataObjectKind(c::readNpz(resdmd).metadataJson()) == c::object_kind::krom && c::isKrom(c::readNpz(resdmd)),
              "ResDMD publishes a krom");
        const auto model = c::loadKoopmanModel(resdmd);
        check(model.liftedDimension() == 2 && model.samples().samples() > 0, "the ResDMD krom has its spectrum and samples");
        // ResDMD's own output survives publication (Colbrook & Townsend, Comm. Pure Appl. Math. 2024):
        // each eigenpair carries its one-step residual (near zero here: the linear dictionary spans this
        // linear system's eigenfunctions) and the spectrum the pseudospectral minimum at each grid point.
        for (std::size_t j = 0; j < 2; ++j) {
            const auto pair = model.eigenpair(j);
            check(pair.residual && pair.residual->value < 1e-6 && pair.residual->definition == "one-step"
                  && pair.residual->scope == "empirical-training" && pair.eigenfunction, "the krom's eigenpairs carry ResDMD's residuals");
        }
        const auto& pseudo = model.spectrum().pseudospectrum();
        check(pseudo && pseudo->grid.size() == 2 && pseudo->minimum.size() == 2 && pseudo->definition == "one-step",
              "the krom keeps the pseudospectrum grid");
        check(pseudo->minimum[1] > 0.5, "z = 0 lies far from this spectrum (|mu| near 1)");
        {   // Stored once: on the spectrum, not in a separate block or the embedded operator contract.
            const auto archive = c::readNpz(resdmd);
            check(archive.contains("eigenpair_residual") && archive.contains("pseudospectrum_grid")
                  && !archive.contains("spectral_data_json"), "the residuals live in the krom's spectrum");
            check(c::readSpectralKromArchive(archive).spectral.residuals.eigenpair.size() == 2,
                  "the spectral fit restores its residuals from the krom");
        }
        const auto rates = model.continuousEigenvalues();
        check(std::abs(rates[0] - expected) < 2e-2 || std::abs(rates[0] - std::conj(expected)) < 2e-2,
              "ResDMD finds the oscillator's eigenvalues");
        check(model.representation().name == "interpolation" && model.representationNames().size() == 2,
              "a delay-window expansion and a state interpolation; the state one is evaluated");
        const double x0[2] = {0.5, -0.2};
        check(std::isfinite(model.evaluate(x0)[0].real()), "the ResDMD krom evaluates at a state");
        bool noMap = false;
        try { (void)model.controlMap(x0, "constant"); }
        catch (const std::out_of_range& error) { noMap = std::string(error.what()).find("no input map B") != std::string::npos; }
        check(noMap, "a krom without B says so when control needs one");

        // 2. H-DMD and SINDy kroms load through the same loader.
        const auto hdmd = postprocess(free, R"([{"name":"hdmd","settings":{"observable":0,"matrix_size":200,
            "sampling_window":100,"rank":2,"window_alignment":"end"}}])", root).archives.front();
        const auto sindy = postprocess(free, R"([{"name":"sindy","settings":{"polynomial_degree":1,"threshold":0,
            "max_iterations":1,"gamma_method":"system_rhs","eigenvalues":"from_data"}}])", root).archives.front();
        for (const auto& path : {hdmd, sindy}) {
            check(c::dataObjectKind(c::readNpz(path).metadataJson()) == c::object_kind::krom, "H-DMD and SINDy publish kroms");
            const auto other = c::loadKoopmanModel(path);
            check(other.liftedDimension() >= 2 && std::isfinite(other.evaluate(x0)[0].real()), "loads like the ResDMD krom");
        }
        {   // SINDy attaches its generator residual to each pair (Kaiser, Kutz & Brunton 2021); H-DMD its one-step residual.
            const auto fitted = c::loadKoopmanModel(sindy);
            for (std::size_t j = 0; j < fitted.liftedDimension(); ++j) {
                const auto& residual = fitted.spectrum().residual(j);
                check(residual && residual->definition == "generator" && std::isfinite(residual->value) && residual->value < 1e-3,
                      "SINDy publishes a generator residual per pair");
            }
            check(!c::readNpz(sindy).contains("training_residual_rms"), "SINDy keeps no parallel residual array");
            const auto delay = c::loadKoopmanModel(hdmd);
            std::size_t visited = 0;
            for (const auto& pair : delay.pairs()) {
                const auto& residual = pair.residual;
                check(residual && residual->definition == "one-step" && residual->scope == "empirical-training"
                      && std::isfinite(residual->value) && residual->value >= 0 && pair.eigenfunction,
                      "H-DMD publishes a one-step residual per pair");
                check(pair.eigenvalue.index == visited && pair.eigenvalue.discrete == delay.spectrum().eigenvalue(visited).discrete,
                      "pairs() visits the pairs in mode order");
                ++visited;
            }
            check(visited == delay.liftedDimension() && delay.pairs().size() == visited, "pairs() visits every pair");
            check(!c::readNpz(hdmd).contains("eigenfunction_koopman_residual"), "H-DMD keeps no parallel residual array");
            check(!delay.spectrum().pseudospectrum() && delay.spectrum().continuous().empty(), "H-DMD has no pseudospectrum");
        }
        {   // control_target lifts a recorded trajectory over its own clock when no time is given.
            const auto lifted = c::readNpz(postprocess(free, R"([{"name":"control_target","settings":{"krom":")"
                + sindy.generic_u8string() + R"(","state_target":{"kind":"trajectory","source":")" + free.generic_u8string()
                + R"("}}}])", root).archives.front());
            const auto clock = c::readNpz(free).at("time").doubles(), times = lifted.at("time").doubles();
            check(std::abs(times.front() - clock.front()) < 1e-12 && std::abs(times.back() - clock.back()) < 2e-3
                  && lifted.at("values").shape.front() == times.size(), "a trajectory target spans its own clock");
        }
        {   // krom_extract keeps the selected pairs with their residuals (a child: they come from the parent).
            const auto extracted = postprocess(resdmd, R"([{"name":"krom_extract","settings":{"krom":")"
                + resdmd.generic_u8string() + R"(","mode_selection":{"count":1,"rule":"oscillatory_pairs"}}}])", root).archives.front();
            check(isChildOf(extracted, resdmd), "the extracted krom references the ResDMD krom");
            const auto reduced = c::loadKoopmanModel(extracted);
            check(reduced.liftedDimension() == 2, "one oscillatory pair");
            for (std::size_t j = 0; j < 2; ++j)
                check(reduced.spectrum().residual(j) && reduced.spectrum().residual(j)->definition == "one-step",
                      "krom_extract keeps the selected pairs' residuals");
            check(reduced.spectrum().pseudospectrum().has_value(), "and the spectrum's pseudospectrum");
        }

        // 3. C on top: observable_reconstruction publishes a child krom with the readout.
        const auto withC = postprocess(free, R"([{"name":"observable_reconstruction","settings":{"krom":")"
            + resdmd.generic_u8string() + R"(","observable":"position"}}])", root).archives.front();
        check(isChildOf(withC, resdmd), "the readout krom references the ResDMD krom");
        const auto childArchive = c::readNpz(withC);
        check(!childArchive.contains("eigenfunction_values") && !childArchive.contains("interpolation_points")
              && childArchive.at("reconstruction_names").string() == "position", "it stores only its readout");
        const auto readout = c::loadKoopmanModel(withC);
        check(readout.hasReconstruction("position") && readout.spectrum().continuousRates() == rates,
              "C over the ResDMD eigenfunctions");
        {   // At a training state (the interpolation is local) C psi reproduces the position.
            const auto source = c::readNpz(free);
            const auto values = (source.contains("values") ? source.at("values") : source.at("states")).doubles();
            const double* state = values.data() + 2 * 150;
            check(std::abs(readout.reconstruct("position", state)[0].real() - state[0]) < 5e-2 * (1 + std::abs(state[0])),
                  "the readout reproduces the position");
        }

        // 4. B on top: input_map on a forced record publishes a child of the readout krom.
        const auto forced = solve(root / "forced", true);
        const auto withB = postprocess(forced, R"([{"name":"input_map","settings":{"krom":")"
            + withC.generic_u8string() + R"(","representation":"constant"}}])", root).archives.front();
        check(isChildOf(withB, withC), "the input-map krom references the readout krom");
        const auto controlled = c::loadKoopmanModel(withB);
        check(controlled.hasControlMap("constant") && controlled.hasReconstruction("position")
              && controlled.inputDimension() == 1, "B and C on the ResDMD krom");
        const auto b = controlled.controlMap(x0, "constant");
        check(b.size() == 2 && std::isfinite(b[0].real()) && std::abs(b[0]) > 0, "a fitted input map");
        check(!c::readNpz(withB).contains("reconstruction_names"), "the input-map krom stores no readout of its own");

        std::error_code ignored;
        fs::remove_all(root, ignored);
        std::cout << "ResDMD krom tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ResDMD krom tests failed: " << error.what() << '\n';
        std::error_code ignored;
        if (!root.empty()) fs::remove_all(root, ignored);
        return 1;
    }
}
