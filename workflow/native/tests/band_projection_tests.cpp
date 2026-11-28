// A system's known spectrum against a spectral fit, from standard pipeline jobs only (one job group):
// a solve of a family of conservative pendulum librations, band_projection of the records onto the
// catalogue band libration_phase (weighted Birkhoff averages in its exact phase), a pooled mpEDMD fit
// with canonical-area weights from the catalogue invariant and band, Rigged DMD on the fit, the
// reference band added beside the Rigged band of the same observable, and packet_validation on held-out
// orbits. Checks: the projection's Fourier coefficients rebuild the angle at the start of each record,
// the reference density integrates to the observable's mass, and both bands sit on one spectrum.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>

#include <algorithm>
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
constexpr double pi = 3.14159265358979323846;
void check(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

std::string solve(const std::string& values) {
    return R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"librations",
      "system":{"name":"nonlinear_pendulum","parameters":{"damping_n_m_s":0}},
      "time":{"start":0,"stop":80,"output_step":0.05},
      "initial_conditions":{"mode":"explicit","count":)" + std::to_string(std::count(values.begin(), values.end(), '[') - 1) + R"(,"clusters":1,"values":)" + values + R"(},
      "solver":{"method":"RK45","rtol":1e-11,"atol":1e-13},"archive":{"enabled":true}})";
}
std::string post(const std::string& label, const std::string& sources, const std::string& method, const std::string& settings) {
    return R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":")" + label + R"(","archive":{"enabled":true},"sources":)"
           + sources + R"(,"methods":[{"name":")" + method + R"(","settings":)" + settings + "}]}";
}

}  // namespace

int main(int, char** argv) {
    fs::path root;
    try {
        root = fs::temp_directory_path() / ("dynsyslab-band-projection-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(root);
        const std::string kernel = R"("sample_interval":0.05,"bandwidth":0.1,"order":2)";
        const std::string weights = R"({"canonical_area":{"invariant":"energy","band":"libration_phase"}})";
        std::string group = std::string(R"({"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"known-spectrum-vs-fit","jobs":[)")
            + R"({"id":"train","job":)" + solve("[[0.6,0],[0.8,0],[1.0,0],[1.2,0],[1.4,0]]") + "},"
            + R"({"id":"held","job":)" + solve("[[0.7,0],[1.1,0]]") + "},"
            + R"({"id":"reference","keep":true,"inputs":{"train":"${train.archives}"},"job":)" + post("reference", R"(["${inputs.train}"])", "band_projection",
                R"({"band":"libration_phase","observables":["angle"],"harmonics":31,"weights":)" + weights
                + R"(,"invariant":"energy",)" + kernel + R"(,"angles":{"count":256},"packet_centers":[0.3]})") + "},"
            + R"({"id":"fit","inputs":{"train":"${train.archives}"},"job":)"
            + post("fit", R"(["${inputs.train}"])", "mpedmd",
                R"({"observables":[0],"window":{"delays":8},"assume_measure_preserving":true,"pool":true,
                    "residual_grid":{"real":[-1.2,1.2,5],"imag":[-1.2,1.2,5]},
                    "dictionary":{"kind":"linear","constant":false},"source_weights":)" + weights + "}") + "},"
            + R"({"id":"rigged","inputs":{"fit":"${fit.archive}"},"job":)" + post("rigged", R"(["${inputs.fit}"])", "rigged_measure",
                R"({"measures":[{"name":"angle","band":"rigged","reconstruction":{"name":"coordinates","output":0},
                    "angles":{"count":256},"bandwidth":0.1,"order":2},
                   {"name":"angle","band":"packets","reconstruction":{"name":"coordinates","output":0},"angles":[0.3],"bandwidth":0.1,"order":2}]})") + "},"
            + R"({"id":"both","keep":true,"inputs":{"train":"${train.archives}","krom":"${rigged.archive}"},"job":)" + post("both", R"(["${inputs.train}"])", "band_projection",
                R"({"band":"libration_phase","observables":["angle"],"harmonics":31,"weights":)" + weights + "," + kernel
                + R"(,"angles":{"count":256},"krom":"${inputs.krom}","measures":[{"observable":"angle","as":"coordinates","output":0}]})") + "},"
            + R"({"id":"packets","keep":true,"inputs":{"held":"${held.archives}","krom":"${rigged.archive}"},"job":)" + post("packets", R"(["${inputs.held}"])", "packet_validation",
                R"({"krom":"${inputs.krom}","observables":[0],"window":{"delays":8},"band":"libration_phase","harmonics":31,
                    "packets":[{"observable":"angle","measure":{"band":"packets","observable":"coordinates","output":0}}]})") + "}";
        // The Viewer draws the spectrum from a job (when the Viewer CLI is built beside this test):
        // the krom's atoms and both bands' densities, and its eigenvalues over the pseudospectrum.
        const bool viewer = fs::is_regular_file(fs::absolute(fs::u8path(argv[0])).parent_path() / "DynSysLabViewerCli.exe");
        const std::string figure = std::string(R"(,{"id":"measures","keep":true,"inputs":{"krom":"${both.archives[1]}"},"job":{"kind":"viewer_figure",)")
            + R"("source":"${inputs.krom}","output":"spectral-measures.png","figure":{"view":"spectral_measures","width":900,"height":700}}})";
        group += (viewer ? figure : std::string()) + "]}";
        const auto result = w::executeJobGroup(json::parse(group), root);
        fs::path reference, both, packets;
        for (const auto& path : result.archives) {
            const auto name = path.filename().u8string();
            if (name.find("band-projection") != std::string::npos && path.parent_path().filename() == "report" && reference.empty()) reference = path;
            if (name.find("band-projection") != std::string::npos && path.parent_path().filename() == "krom") both = path;
            if (name.find("packet-validation") != std::string::npos) packets = path;
        }
        check(!reference.empty() && !both.empty() && !packets.empty(), "every stage published");
        if (viewer) check(std::any_of(result.archives.begin(), result.archives.end(), [](const fs::path& p) {
                              return p.filename() == "spectral-measures.png" && fs::file_size(p) > 10000; }),
                          "the Viewer rendered the spectral measures");
        // The coefficients rebuild the angle at phase zero of orbit 0 (starting at rest: phase 0, angle 0.6).
        const auto r = c::readNpz(reference);
        const auto a = r.at("fourier_coefficients").complexDoubles();
        const auto h = r.at("fourier_harmonics").doubles().size();
        // the libration phase is pi/2 at the turning point theta = 0.6 (alpha = pi/2 there)
        std::complex<double> rebuilt{};
        for (std::size_t k = 0; k < h; ++k) rebuilt += a[k] * std::exp(std::complex<double>(0, (double(k) - double(h / 2)) * pi / 2));
        check(std::abs(rebuilt - 0.6) < 1e-6, "the band projection rebuilds the angle at the record's start: " + std::to_string(rebuilt.real()));
        double integral = 0, mass = 0;
        for (const double d : r.at("reference_density__angle").doubles()) integral += d * 2 * pi / 256;
        const auto weight = r.at("orbit_weight").doubles();
        for (std::size_t j = 0; j < weight.size(); ++j)
            for (std::size_t k = 0; k < h; ++k) mass += weight[j] * std::norm(a[j * h + k]);
        check(std::abs(integral - mass) < 1e-6 * mass, "the reference density integrates to the observable's mass");
        // Both bands on one spectrum, for one observable.
        const auto model = c::loadKoopmanModel(both);
        check(model.spectrum().band("fourier-reference") && model.spectrum().band("rigged")
              && model.spectrum().decomposition("coordinates").measures.size() == 3, "reference and Rigged bands on one spectrum");
        check(model.spectrum().pseudospectrum() && model.spectrum().pseudospectrum()->grid.size() == 25, "the residual lattice");
        const auto p = c::readNpz(packets);
        check(p.at("packet_rms__angle__0").doubles().size() == 2 && p.at("reference_packet_rms__angle__0").doubles().size() == 2,
              "packets validated on the held-out orbits");
        std::cout << "band projection tests passed\n";
    } catch (const std::exception& error) {
        std::cerr << "band projection tests failed: " << error.what() << '\n';
        std::error_code ignored;
        if (!root.empty()) fs::remove_all(root, ignored);
        return 1;
    }
    std::error_code ignored;
    fs::remove_all(root, ignored);
    return 0;
}
