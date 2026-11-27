// Phase D identification in the job pipeline: a checksum-verified IC table, pooled H-DMD and pooled
// SINDy over several trajectories (each publishing a krom), SINDy on a state-space domain with a
// radial dictionary, the H-DMD -> SINDy readout transfer gate, the observable spectrum comparison,
// delay taps at physical-time lags, HAVOK derivative order, weak-form vector-field SINDy, the SVD
// driver and analysis in archive coordinates.
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
namespace w = dynsyslab::workflow;
namespace fs = std::filesystem;
namespace json = c::json;

void check(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

w::PostprocessResult post(const fs::path& source, const std::string& methods, const fs::path& root) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"identification","methods":)" + methods + "}");
    job.at("sources") = json::Value::Array{source.generic_u8string()};
    return w::executePostprocess(job, root, {});
}

json::Value meta(const fs::path& path) { return json::parse(c::readNpz(path).metadataJson()); }

}  // namespace

int main() {
    fs::path root;
    try {
        root = fs::temp_directory_path() / ("dynsyslab-identification-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(root);
        // 1. Initial conditions from a checksum-verified table (states as columns, on the unit circle).
        const auto table = root / "ics.csv";
        { std::ofstream out(table); out << "1,0,0.6\n0,1,0.8\n"; }
        auto solve = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve",
            "system":{"name":"duffing","parameters":{"alpha":-1,"beta":0,"delta":0.2}},
            "time":{"start":0,"stop":12,"output_step":0.02},
            "initial_conditions":{"mode":"file","count":3,"clusters":1,"layout":"columns","norm":1.0},
            "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12},"archive":{"enabled":true}})");
        solve.at("initial_conditions").object()["path"] = table.generic_u8string();
        solve.at("initial_conditions").object()["sha256"] = c::sha256File(table);
        const auto runs = w::executeJob(solve, root).archives;
        check(runs.size() == 3, "three trajectories from the IC table");
        {
            auto stale = solve;
            stale.at("initial_conditions").object()["sha256"] = std::string(64, '0');
            bool refused = false;
            try { (void)w::executeJob(stale, root / "stale"); } catch (const std::exception&) { refused = true; }
            check(refused, "an IC table that fails its checksum is refused");
        }
        const std::string extra = "[\"" + runs[1].generic_u8string() + "\",\"" + runs[2].generic_u8string() + "\"]";
        const std::complex<double> lambda(-0.1, std::sqrt(0.99));

        // 2. Pooled H-DMD over the three trajectories: one krom, no window across a boundary.
        const auto pooled = post(runs[0], R"([{"name":"hdmd","settings":{"observable":0,"matrix_size":200,"sampling_window":100,
            "rank":2,"window_alignment":"end","additional_sources":)" + extra + "}}]", root).archives.front();
        check(c::dataObjectKind(c::readNpz(pooled).metadataJson()) == c::object_kind::krom, "pooled H-DMD publishes a krom");
        const auto hdmd = c::loadKoopmanModel(pooled);
        check(meta(pooled).at("pooled").at("trajectories").number() == 3, "the krom records its three trajectories");
        check(std::abs(hdmd.continuousEigenvalues()[0] - lambda) < 1e-3 || std::abs(hdmd.continuousEigenvalues()[0] - std::conj(lambda)) < 1e-3,
              "pooled H-DMD finds the oscillator's eigenvalues");
        check(hdmd.samples().samples() == 3 * 200 + 1, "every window of every trajectory is a sample");

        // 3. SINDy over the same trajectories, anchored to the pooled H-DMD krom, with its readouts transferred.
        const auto sindy = post(runs[0], R"([{"name":"sindy","settings":{"krom":")" + pooled.generic_u8string()
            + R"(","polynomial_degree":1,"include_constant":false,"threshold":0,"max_iterations":1,"gamma_method":"system_rhs",
            "readout_transfer":true,"rescale":[1.0,2.0],"additional_sources":)" + extra + "}}]", root).archives.front();
        const auto sindyMeta = meta(sindy);
        check(sindyMeta.at("sindy_training").at("input").string() == "trajectories"
              && sindyMeta.at("sindy_training").at("trajectory_samples").array().size() == 3, "SINDy pools the trajectories");
        check(sindyMeta.at("readout_transfer").at("passed").boolean(), "the H-DMD readouts pass the transfer gate");
        const auto transferred = c::loadKoopmanModel(sindy);
        check(!hdmd.reconstructions().empty() && transferred.hasReconstruction(hdmd.reconstructions().begin()->first)
              && transferred.representation().name == "polynomial",
              "the SINDy krom keeps the transferred readout and evaluates its polynomial");

        // 4. SINDy on a state-space domain with a radial dictionary (rates from the vector field).
        const auto domain = post(runs[0], R"([{"name":"sindy","settings":{"eigenvalues":[[-0.1,0.99498743710662],[-0.1,-0.99498743710662]],
            "domain":{"lower":[-1,-1],"upper":[1,1],"density":9},"gamma_method":"system_rhs","threshold":0,"max_iterations":1,
            "dictionary":{"kind":"gaussian","width":0.8,"include_constant":true,"centers":[[-1,-1],[-1,1],[1,-1],[1,1],[0,0],[0.5,0],[0,0.5],[-0.5,0],[0,-0.5]]},
            "reconstruction":"none"}}])", root).archives.front();
        const auto onDomain = c::loadKoopmanModel(domain);
        check(meta(domain).at("sindy_training").at("input").string() == "domain"
              && onDomain.expansion().dictionary().kinds().back() == "gaussian", "domain SINDy fits in the radial dictionary");

        // 5. The spectra of two observables of one trajectory: one hdmd stage per observable with identical
        //    settings, then observable_spectrum_comparison of their kroms.
        const auto fitOf = [&](int column) {
            return post(runs[0], R"([{"name":"hdmd","settings":{"observable":)" + std::to_string(column)
                + R"(,"matrix_size":200,"sampling_window":100,"rank":2}}])", root).archives.front();
        };
        auto compare = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"spectra",
            "methods":[{"name":"observable_spectrum_comparison","settings":{"names":["x","v"]}}]})");
        compare.at("sources") = json::Value::Array{fitOf(0).generic_u8string(), fitOf(1).generic_u8string()};
        const auto report = c::readNpz(w::executePostprocess(compare, root, {}).archives.front());
        check(report.at("pairwise_discrete_nearest_rms").doubles()[1] < 1e-3, "two observables of one oscillator give nearby spectra");

        // 6. Delay taps at physical-time lags, then H-DMD on them; HAVOK with fourth-order derivatives.
        const auto taps = post(runs[0], R"([{"name":"coordinate_transform","settings":{"kind":"delay_taps","columns":[0],"lags":[0,0.05,0.13]}}])", root).archives.front();
        const auto tapped = c::readNpz(taps);
        check(tapped.at("values").shape[1] == 3 && meta(taps).at("delay_taps").at("lags").array().size() == 3, "three taps of one channel");
        const auto havok = post(runs[0], R"([{"name":"havok","settings":{"delays":30,"rank":4,"observable":0,"derivative_order":4}}])", root);
        check(meta(havok.archives.front()).at("derivative_order").number() == 4, "HAVOK records its derivative order");

        // 7. The SVD driver and archive coordinates are per-stage settings.
        const auto jacobi = post(runs[0], R"([{"name":"dmd","settings":{"rank":2,"svd_driver":"gesvd","analysis_coordinates":"archive"}}])", root);
        const auto jacobiModel = c::loadKoopmanModel(jacobi.archives.front());
        check(std::abs(jacobiModel.continuousEigenvalues()[0].real() + 0.1) < 1e-3, "DMD with the Jacobi SVD driver");
        check(meta(jacobi.archives.front()).at("analysis_coordinates").string() == "archive", "archive coordinates are recorded");
        check(c::readNpz(jacobi.archives.front()).contains("quality_factor"), "the spectrum carries its rate, frequency and Q summaries");

        // 8. POD modes as functions of the state: at a sample state the interpolated amplitudes
        //    reconstruct it exactly (full rank).
        {
            const auto first = c::readNpz(runs[0]).at("values").doubles();
            const auto x0 = json::dump(json::Value::Array{first[0], first[1]});
            const auto pod = c::readNpz(post(runs[0], R"([{"name":"pod","settings":{"rank":2,"amplitude_queries":[)" + x0
                + "]}}]", root).archives.front());
            const auto rebuilt = pod.at("pod_reconstruction_queries").doubles();
            check(pod.at("pod_amplitude_queries").shape == std::vector<std::size_t>{1, 2}
                  && std::abs(rebuilt[0] - first[0]) < 1e-9
                  && std::abs(rebuilt[1] - first[1]) < 1e-9, "POD amplitudes reconstruct a sample state");
        }

        // 9. Delay pairs (midpoint, secant) at lags 0 and 5 samples; their delay consistency is small,
        //    since the pairs are second-order estimates of (value, rate) at shifted times.
        {
            const auto pairs = post(runs[0], R"([{"name":"delay_pairs","settings":{"lags":[0,5],"half_width":1}}])", root).archives.front();
            const auto pairArchive = c::readNpz(pairs);
            check(pairArchive.at("values").shape[1] == 4 && meta(pairs).at("delay_pairs").at("latency").number() > 0, "delay pairs carry their latency");
            const double dt = pairArchive.at("time").doubles()[1] - pairArchive.at("time").doubles()[0];
            const auto report = post(pairs, R"([{"name":"delay_consistency","settings":{"lags":[0,)" + json::dump(json::Value(5 * dt))
                + "]}}]", root).archives.front();
            check(c::readNpz(report).at("value_relative_error").doubles()[0] < 1e-3, "secant pairs are delay consistent");
        }

        std::error_code ignored;
        fs::remove_all(root, ignored);
        std::cout << "identification workflow tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "identification workflow tests failed: " << error.what() << '\n';
        std::error_code ignored;
        if (!root.empty()) fs::remove_all(root, ignored);
        return 1;
    }
}
