// SINDy takes what it is given: H-DMD's eigenfunctions (it adds a calibrated polynomial to the same
// eigenfunctions), or training data alone (it finds the eigenvalues too). A damped linear oscillator,
// x'' + 0.2 x' + x = 0 (Duffing with alpha = -1, beta = 0 in this suite's sign convention), has
// eigenvalues -0.1 ± i sqrt(0.99) and linear eigenfunctions, so every representation can match it;
// a hardening Duffing oscillator (beta = 1) checks the default fit where they cannot.
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace dynsyslab;
namespace json = core::json;
namespace {
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

std::filesystem::path eigenfunctionsWithBase(const workflow::PostprocessResult& result) {
    for (const auto& path : result.archives) {
        const auto meta = json::parse(core::readNpz(path).metadataJson());
        if (!meta.at("krom").at("parent").isNull()) return path;
    }
    throw std::runtime_error("SINDy published no child of H-DMD's krom");
}

workflow::PostprocessResult postprocess(const std::filesystem::path& source, const std::string& methods,
                                        const std::filesystem::path& root) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"sindy-inputs","methods":)"
                           + methods + "}");
    job.at("sources") = json::Value::Array{source.generic_u8string()};
    return workflow::executePostprocess(job, root, {});
}
}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-sindy-inputs-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto solve = workflow::executeJob(json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve",
            "system":{"name":"duffing","parameters":{"alpha":-1,"beta":0,"delta":0.2}},
            "time":{"start":0,"stop":20,"output_step":0.01},
            "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[1.0,0.0]]},
            "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12},"archive":{"enabled":true}})"), root);
        const auto trajectory = solve.archives.front();
        const std::complex<double> lambda(-0.1, std::sqrt(0.99));
        const std::string sindy = R"({"name":"sindy","settings":{"polynomial_degree":2,"threshold":0,"max_iterations":1,
            "gamma_method":"system_rhs","reconstruction":"none")";

        // H-DMD -> SINDy: one set of eigenfunctions with three representations.
        const auto chain = postprocess(trajectory, R"([{"name":"hdmd","settings":{"observable":0,"matrix_size":300,
            "sampling_window":200,"rank":2,"window_alignment":"end"}},)" + sindy + "}}]", root);
        const auto combinedPath = eigenfunctionsWithBase(chain);
        const auto combined = core::readNpz(combinedPath);
        const auto set = core::readKoopmanModel(combined);
        check(set.representationNames() == std::vector<std::string>{"polynomial", "interpolation", "time-delay"},
              "the polynomial joins H-DMD's interpolation and time-delay representations, as the default");
        check(!combined.contains("window_map") && !combined.contains("eigenfunction_values"),
              "H-DMD's window map and samples are kept by reference, not copied");
        for (const double mismatch : combined.at("agreement_mismatch").doubles())
            check(mismatch < 0.05, "the calibrated polynomial agrees with H-DMD's values");
        const auto rates = set.spectrum().continuousRates();
        check(std::any_of(rates.begin(), rates.end(), [&](auto r) { return std::abs(r - lambda) < 1e-2; }),
              "the eigenvalues are H-DMD's");

        // K-ROMs on these eigenfunctions use the representation asked for (default: the polynomial).
        const auto krom = [&](const std::string& representation) {
            return postprocess(combinedPath, R"([{"name":"krom_extract","settings":{"selected_modes":[0,1])"
                + (representation.empty() ? std::string() : ",\"eigenfunction_representation\":\"" + representation + "\"")
                + "}}]", root).archives.front();
        };
        const auto named = [&](const std::filesystem::path& path) {
            return json::parse(core::readNpz(path).metadataJson()).at("eigenfunction_representation_name").string();
        };
        check(named(krom("")) == "polynomial" && named(krom("interpolation")) == "interpolation",
              "K-ROMs record the representation they were built on");
        {   // A time-delay K-ROM is a K-ROM; it evaluates on windows only, so a state is refused.
            const auto delay = krom("time-delay");
            check(named(delay) == "time-delay", "the time-delay K-ROM records its representation");
            bool refused = false;
            const double state[2] = {1.0, 0.0};
            try { (void)core::loadKoopmanModel(delay).evaluate(state); } catch (const std::exception&) { refused = true; }
            check(refused, "a time-delay K-ROM refuses a state (no state map)");
        }

        // A nonlinear oscillator whose H-DMD values cover only a late window near rest (as in a
        // decaying RPE bubble): SINDy's default fit is the generator null space (Kaiser, Kutz &
        // Brunton 2021), so its polynomials stay eigenfunctions wherever the references reach, and a
        // mode that disagrees with its references keeps unit RMS over the training rows, not their scale.
        {
            const auto hard = workflow::executeJob(json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve",
                "system":{"name":"duffing","parameters":{"alpha":-1,"beta":1,"delta":0.2}},
                "time":{"start":0,"stop":30,"output_step":0.01},
                "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[2.0,0.0]]},
                "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12},"archive":{"enabled":true}})"), root).archives.front();
            const auto fitted = eigenfunctionsWithBase(postprocess(hard, R"([{"name":"hdmd","settings":{"observable":0,
                "matrix_size":300,"sampling_window":1500,"rank":0,"svd_tolerance":1e-10,"window_alignment":"end"}},
                {"name":"sindy","settings":{"polynomial_degree":5,"threshold":0.05,"max_iterations":20,
                "gamma_method":"system_rhs","reconstruction":"none"}}])", root));
            const auto archive = core::readNpz(fitted);
            const auto meta = json::parse(archive.metadataJson());
            check(meta.at("sindy_input").at("fit").string() == "nullspace", "SINDy's default fit is the generator null space");
            const auto model = core::loadKoopmanModel(fitted);
            const auto states = core::readNpz(hard).at("values").doubles();
            const std::size_t rows = states.size() / 2, modes = model.liftedDimension();
            std::vector<double> squared(modes, 0.0);
            for (std::size_t k = 0; k < rows; ++k) {
                const auto values = model.evaluate(states.data() + 2 * k);
                for (std::size_t j = 0; j < modes; ++j) squared[j] += std::norm(values[j]);
            }
            const auto rates = model.spectrum().continuousRates();
            std::vector<double> residual;  // the pairs' generator residuals, of the published eigenfunctions
            for (std::size_t j = 0; j < modes; ++j) residual.push_back(model.spectrum().residual(j).value().value);
            std::size_t dominant = 0;  // the mode nearest the linearization's -0.1 + i sqrt(0.99)
            for (std::size_t j = 1; j < modes; ++j)
                if (std::abs(rates[j] - lambda) < std::abs(rates[dominant] - lambda)) dominant = j;
            // Relative generator residual of the published eigenfunction (both residual and RMS after
            // calibration; the former training_residual_rms was the uncalibrated fit's, which made this
            // ratio look 1/|c| smaller). Measured: nullspace 0.37, anchored 0.85.
            check(residual[dominant] / std::sqrt(squared[dominant] / double(rows)) < 0.5,
                  "the oscillatory eigenfunction satisfies the generator equation better than the anchored fit (relative residual below 0.5)");
            const auto& flagged = meta.at("agreement").at("flagged_modes").array();
            check(!flagged.empty(), "some H-DMD modes disagree with every polynomial");
            for (const auto& mode : flagged) {
                const std::size_t j = std::size_t(mode.number());
                check(meta.at("agreement").at("calibration").array().at(j).string() == "normalized"
                      && std::abs(std::sqrt(squared[j] / double(rows)) - 1.0) < 1e-9,
                      "a mode that disagrees with its references is normalized to unit RMS");
            }
        }

        // SINDy from training data alone: it finds the eigenvalues as well.
        const auto alone = postprocess(trajectory, "[" + sindy + R"(,"eigenvalues":"from_data"}}])", root);
        const auto found = core::readKoopmanModel(core::readNpz(alone.archives.front())).spectrum().continuousRates();
        check(std::any_of(found.begin(), found.end(), [&](auto r) { return std::abs(r - lambda) < 1e-6; })
              && std::any_of(found.begin(), found.end(), [&](auto r) { return std::abs(r - std::conj(lambda)) < 1e-6; }),
              "eigenvalues from data: -0.1 ± i sqrt(0.99)");
        std::filesystem::remove_all(root);
        std::cout << "SINDy input tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SINDy input tests failed: " << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
