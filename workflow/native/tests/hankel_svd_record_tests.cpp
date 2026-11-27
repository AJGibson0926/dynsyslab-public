// H-SVD takes time-series records and emits a Hankel basis plus a coordinate time series that
// H-DMD consumes directly; the state domain of that H-DMD resolves to the original solve by time.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/HankelBasis.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value job(const std::filesystem::path& source, const char* method, json::Value::Object settings) {
    return json::Value::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", method},
        {"sources", json::Value::Array{source.generic_u8string()}},
        {"methods", json::Value::Array{json::Value::Object{{"name", method}, {"enabled", true}, {"settings", settings}}}},
        {"archive", json::Value::Object{{"enabled", true}}}};
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-hsvd-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto sourcePath = root / "trajectory.npz";
        constexpr std::size_t samples = 80;
        std::vector<double> time(samples), states;
        for (std::size_t k = 0; k < samples; ++k) {
            const double t = 0.05 * double(k);
            time[k] = t;
            states.push_back(std::cos(1.7 * t) + 0.2 * std::cos(3.1 * t));
            states.push_back(-1.7 * std::sin(1.7 * t) - 0.62 * std::sin(3.1 * t));
        }
        core::NpzRecord source;
        source.arrays["time"] = core::NpyArray::fromDoubles({samples}, time);
        source.arrays["values"] = core::NpyArray::fromDoubles({samples, 2}, states);
        source.arrays["metadata_json"] = core::NpyArray::fromString(
            R"({"schema":"dynsyslab.trajectory/v1","label":"source","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{},"integration_form":"dimensionless"})");
        core::writeNpzAtomic(sourcePath, source);

        // H-SVD on position only, causal windows, whitened coordinates.
        const auto svdRun = dynsyslab::workflow::executeJob(job(sourcePath, "hankel_svd", {
            {"observables", json::Value::Array{"state_0"}}, {"window", json::Value::Object{{"delays", 12.0}}},
            {"window_alignment", "end"}, {"rank", 4.0}, {"coordinate_scaling", "whitened"}}), root / "svd");
        const auto svd = core::readNpz(svdRun.archives.front());
        const auto meta = json::parse(svd.metadataJson());
        require(meta.at("hankel_basis").at("schema").string() == core::kHankelBasisSchema, "basis block written");
        require(svd.at("values").shape == std::vector<std::size_t>({samples - 11, 4}), "one coordinate row per window");
        require(std::abs(svd.at("time").doubles().front() - time[11]) < 1e-12, "causal windows are labelled by their newest sample");

        // The stored coordinates are the basis applied to the record's windows.
        const auto basis = core::hankelBasisFromArchive(meta.at("hankel_basis"), svd.at("basis_modes").doubles(),
                                                        svd.at("singular_values").doubles(), svd.at("mean").doubles());
        std::vector<double> positions(samples);
        for (std::size_t k = 0; k < samples; ++k) positions[k] = states[2 * k];
        const core::TimeSeries record(time, states, 2);
        const auto q = basis.coordinates(record, {30}, core::HankelCoordinateScaling::whitened);
        const auto stored = svd.at("values").doubles();
        for (std::size_t k = 0; k < 4; ++k)
            require(std::abs(q[k] - stored[(30 - 11) * 4 + k]) < 1e-9 * (1 + std::abs(q[k])), "basis reproduces stored coordinates");

        // H-DMD consumes the coordinate record as a time series; its state domain is the original solve.
        const auto dmdRun = dynsyslab::workflow::executeJob(job(svdRun.archives.front(), "hdmd", {
            {"observables", json::Value::Array{"hsvd[0]", "hsvd[1]"}}, {"matrix_size", 20.0}, {"sampling_window", 4.0},
            {"rank", 4.0}, {"svd_tolerance", 1.0e-12}}), root / "dmd");
        const auto dmd = core::readNpz(dmdRun.archives.front());
        const auto set = core::readKoopmanModel(dmd);
        require(core::domainName(set.representation("interpolation").function->coordinates()) == "state", "the full state is the default domain");
        const auto& points = set.interpolation().points();
        const auto t0 = set.samples().time.front();
        const auto row = std::size_t(std::lround(t0 / 0.05));
        require(std::abs(points[0] - states[2 * row]) < 1e-12, "base points come from the solve, found by time");
        {   // Asked for, the interpolation is on H-DMD's own input channels instead.
            const auto onInput = dynsyslab::workflow::executeJob(job(svdRun.archives.front(), "hdmd", {
                {"observables", json::Value::Array{"hsvd[0]", "hsvd[1]"}}, {"matrix_size", 20.0}, {"sampling_window", 4.0},
                {"rank", 4.0}, {"svd_tolerance", 1.0e-12}, {"interpolation_domain", "input"}}), root / "dmd-input");
            const auto inputSet = core::readKoopmanModel(core::readNpz(onInput.archives.front()));
            require(inputSet.recognizeDomain({"hsvd[0]", "hsvd[1]"}) == "hsvd[0],hsvd[1]", "H-SVD coordinates name the input domain");
        }
        const auto errorsRun = dynsyslab::workflow::executeJob(job(dmdRun.archives.front(), "eigenfunction_errors", {
            {"interpolation_methods", json::Value::Array{json::Value::Object{{"name", "inverse_distance"}},
                                                         json::Value::Object{{"name", "local_linear"}, {"neighbors", 12.0}}}},
            {"include_expansion", true}}), root / "errors");
        const auto errors = core::readNpz(errorsRun.archives.front());
        require(errors.at("value_error").shape == std::vector<std::size_t>({3, set.liftedDimension()}), "one error row per method plus the expansion");
        require(json::parse(errors.metadataJson()).at("rows").array().size() == 3, "rows are named");
        std::filesystem::remove_all(root);
        std::cout << "hankel svd record tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
