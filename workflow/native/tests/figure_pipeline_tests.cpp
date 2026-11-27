// figure_dataset operations, bifurcation_points, viewer_figure curve sources and the phase-cover
// clock check.
#include "../src/ViewerFigureJob.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace fs = std::filesystem;
namespace core = dynsyslab::core;
namespace json = dynsyslab::core::json;
namespace workflow = dynsyslab::workflow;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }

fs::path run(const fs::path& root, const std::string& name, double amplitude, double dt) {
    std::vector<double> time, values;
    for (int i = 0; i < 400; ++i) {
        time.push_back(i * dt);
        values.push_back(amplitude * std::cos(time.back()));
        values.push_back(-amplitude * std::sin(time.back()));
    }
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({time.size()}, time);
    record.arrays["values"] = core::NpyArray::fromDoubles({time.size(), 2}, values);
    record.arrays["metadata_json"] = core::NpyArray::fromString(
        R"({"schema":"dynsyslab.solution/v1","system":"duffing","system_parameters":{"driver":{"amplitude":)"
        + std::to_string(amplitude) + "}}}");
    const auto path = root / (name + ".npz");
    core::writeNpzAtomic(path, record);
    return path;
}

json::Value postprocess(const std::string& method, const json::Value& settings, const std::vector<fs::path>& sources) {
    json::Value::Array list;
    for (const auto& source : sources) list.emplace_back(source.generic_u8string());
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.compute-job/v1")},
        {"kind", json::Value("postprocess")}, {"label", json::Value(method + "-test")}, {"sources", json::Value(list)},
        {"methods", json::Value(json::Value::Array{json::Value(json::Value::Object{{"name", json::Value(method)},
            {"settings", settings}})})}, {"archive", json::Value(json::Value::Object{{"enabled", json::Value(true)}})}});
}

}  // namespace

int main() {
    const auto root = fs::temp_directory_path() / ("dynsyslab-figure-pipeline-"
        + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        fs::create_directories(root);
        const auto a = run(root, "a", 1.0, 0.05), b = run(root, "b", 2.0, 0.025);
        const auto in = json::parse(R"({"run":")" + a.generic_u8string() + R"(","other":")" + b.generic_u8string() + "\"}");
        auto settings = json::parse(R"({"arrays":[
            {"name":"radius","from":"run","array":"values","columns":[0,1],"operation":"norm"},
            {"name":"gap","from":"run","array":"values","column":0,"operation":"difference",
             "minus":{"from":"other","array":"values","column":0,"onto":{"from":"run","array":"time"}}},
            {"name":"slope","from":"run","array":"values","column":0,"operation":"diff","per_time":true}]})");
        settings.at("inputs") = in;
        const auto table = core::readNpz(workflow::executePostprocess(postprocess("figure_dataset", settings, {a}), root).archives.front());
        const auto radius = table.at("radius").doubles(), gap = table.at("gap").doubles(), slope = table.at("slope").doubles();
        require(radius.size() == 400 && std::abs(radius[17] - 1.0) < 1e-12, "norm across columns");
        require(std::abs(gap[100] - (std::cos(5.0) - 2 * std::cos(5.0))) < 1e-3, "difference onto one clock");
        require(slope.size() == 399 && std::abs(slope[200] + std::sin(10.025)) < 1e-3, "diff per time is the derivative");

        bool rejected = false;
        const auto points = core::readNpz(workflow::executePostprocess(postprocess("bifurcation_points",
            json::parse(R"({"parameter":"system_parameters.driver.amplitude","tail_fraction":0.5})"), {a, b}), root).archives.front());
        const auto parameter = points.at("parameter").doubles(), maximum = points.at("maximum").doubles();
        require(!parameter.empty() && parameter.size() == maximum.size(), "one parameter per maximum");
        for (std::size_t i = 0; i < parameter.size(); ++i)
            require(std::abs(maximum[i] - parameter[i]) < 1e-3, "the maxima of A cos t are A");

        auto figure = json::parse(R"({"kind":"viewer_figure","output":"x.png","figure":{"curves":[
            {"series":"values","label":"a"},{"series":"values","label":"b","source":")" + b.generic_u8string() + R"("}]}})");
        figure.at("source") = json::Value(a.generic_u8string());
        require(workflow::detail::hasCurveSources(figure), "per-curve sources are detected");
        const auto scratch = root / "curves.npz";
        const auto resolved = workflow::detail::assembleCurveSources(figure, scratch);
        const auto combined = core::readNpz(scratch);
        require(resolved.at("source").string() == scratch.u8string(), "the figure reads the scratch archive");
        require(combined.at("curve1_values").shape == std::vector<std::size_t>({400, 2}), "curves share the figure clock");
        require(std::abs(combined.at("curve1_values").doubles()[2 * 100] - 2 * std::cos(5.0)) < 1e-3, "resampled onto the clock");
        require(std::isnan(combined.at("curve1_values").doubles()[2 * 399]), "NaN outside the curve's own span");
        auto page = resolved;  // a model-scale page: the ordinate fits curve 0 (amplitude 1) only
        page.at("figure").at("fit_y_to") = json::parse("[0]");
        const auto fitted = workflow::detail::fitOrdinate(page).at("figure").at("y_range").array();
        require(fitted.size() == 2 && fitted[1].number() < 1.2 && fitted[0].number() > -1.2, "fit_y_to ignores the larger curve");

        const auto spectrum = [&](const std::string& name, double shift, int matrixSize) {
            core::NpzRecord record;
            record.arrays["continuous_eigenvalues"] = core::NpyArray::fromComplex({2}, {{-0.1, 1.0}, {-0.1 + shift, -1.0}});
            record.arrays["discrete_eigenvalues"] = core::NpyArray::fromComplex({2}, {{0.9, 0.1}, {0.9 + shift, -0.1}});
            record.arrays["metadata_json"] = core::NpyArray::fromString(R"({"method":"hdmd","observable":")" + name
                + R"(","matrix_size":)" + std::to_string(matrixSize) + R"(,"sampling_window":10,"sample_stride":1,"delta_t":0.1})");
            const auto path = root / (name + ".spectrum.npz");
            core::writeNpzAtomic(path, record);
            return path;
        };
        const auto x = spectrum("x", 0.0, 40), v = spectrum("v", 0.3, 40);
        const auto compared = core::readNpz(workflow::executePostprocess(postprocess("observable_spectrum_comparison",
            json::Value(json::Value::Object{}), {x, v}), root).archives.front());
        const auto rms = compared.at("pairwise_discrete_nearest_rms").doubles();
        require(rms.size() == 4 && rms[0] == 0.0 && rms[1] == rms[2], "symmetric distance matrix");
        require(std::abs(rms[1] - std::sqrt(0.5 * (0.04 / 2 + 0.09 / 2))) < 1e-12, "nearest-neighbour RMS of one moved eigenvalue");
        rejected = false;
        try { (void)workflow::executePostprocess(postprocess("observable_spectrum_comparison", json::Value(json::Value::Object{}),
                                                             {x, spectrum("w", 0.0, 60)}), root); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "spectra fitted with different settings are not compared");

        {   // control_score estimation scores on a shadow run (no reference signal needed)
            core::NpzRecord shadow;
            std::vector<double> t, x, te, xh, innovation, y, applied, commanded;
            for (int i = 0; i <= 1000; ++i) { t.push_back(0.01 * i); x.push_back(std::sin(t.back())); x.push_back(std::cos(t.back())); }
            for (int i = 0; i <= 1000; ++i) {
                te.push_back(0.01 * i);
                xh.push_back(std::sin(te.back()) + 0.1); xh.push_back(std::cos(te.back()));  // bias 0.1 on component 0
                innovation.push_back(i % 2 ? 1.0 : -1.0);  // lag-1 autocorrelation -1
                y.push_back(std::sin(te.back()));
                applied.push_back(1.0 + std::sin(te.back())); commanded.push_back(1.1 + std::sin(te.back()));
            }
            shadow.arrays["time"] = core::NpyArray::fromDoubles({t.size()}, t);
            shadow.arrays["values"] = core::NpyArray::fromDoubles({t.size(), 2}, x);
            shadow.arrays["estimator_time"] = core::NpyArray::fromDoubles({te.size()}, te);
            shadow.arrays["state_estimate"] = core::NpyArray::fromDoubles({te.size(), 2}, xh);
            shadow.arrays["innovation"] = core::NpyArray::fromDoubles({te.size(), 1}, innovation);
            shadow.arrays["measurement_values"] = core::NpyArray::fromDoubles({te.size(), 1}, y);
            shadow.arrays["applied_input"] = core::NpyArray::fromDoubles({te.size(), 1}, applied);
            shadow.arrays["commanded_input"] = core::NpyArray::fromDoubles({te.size(), 1}, commanded);
            shadow.arrays["metadata_json"] = core::NpyArray::fromString(R"({"schema":"dynsyslab.solution/v1"})");
            const auto shadowPath = root / "shadow.npz";
            core::writeNpzAtomic(shadowPath, shadow);
            const auto scores = core::readNpz(workflow::executePostprocess(postprocess("control_score", json::parse(
                R"({"windows":{"late":[4,10]},"estimation":{"components":[0,1],"innovation_lags":[0.01],"shadow_force":true}})"),
                {shadowPath}), root).archives.front());
            require(std::abs(scores.at("estimate_bias_0").doubles()[0] - 0.1) < 1e-9 && std::abs(scores.at("estimate_error_rms_0").doubles()[0] - 0.1) < 1e-9,
                    "estimate bias and rms error per component");
            require(std::abs(scores.at("estimate_error_rms_1").doubles()[0]) < 1e-9, "an exact component has no error");
            require(std::abs(scores.at("innovation_autocorrelation_0.01").doubles()[0] + 1.0) < 1e-2, "alternating innovations are anti-correlated");
            require(std::abs(scores.at("shadow_force_error_rms").doubles()[0] - 0.1) < 1e-9
                    && std::abs(scores.at("shadow_force_correlation").doubles()[0] - 1.0) < 1e-9, "shadow force error and correlation");
        }

        auto query = json::parse(R"({"schema":"dynsyslab.analytic-query/v1","operation":"phase_cover",
            "system":{"name":"isotonic_oscillator"},"initial_states":[[1.4,0]],"prime":2,"samples":64})");
        const auto cover = workflow::evaluateAnalyticQuery(query);
        require(cover.at("schema").string() == "dynsyslab.phase-cover-geometry/v1" && cover.at("orbits").array().size() == 1,
                "a one-level cover");
        query.at("initial_states") = json::parse("[[1.4,0],[2.0,0]]");  // an isochronous flow has no 2-cover
        rejected = false;
        try { (void)workflow::evaluateAnalyticQuery(query); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "a failed clock check is an error");

        // The case-study experiments that replaced the shadow and tuning scripts plan: their sweeps read
        // object-valued points (${sweep.setup.krom}) and their stages expand.
        const auto study = fs::path(__FILE__).parent_path().parent_path().parent_path().parent_path()
            / "examples/native/case-studies/rpe-echo-mode-selection-20261116";
        for (const auto* file : {"closed-loop/shadow-experiment.json", "tuning/tuning-experiment.json"}) {
            const auto path = study / file;
            const auto plan = json::dump(workflow::planJob(workflow::bindJobFile(json::load(path.wstring()), path)));
            require(plan.find("${sweep.") == std::string::npos, "every sweep reference resolves");
            require(plan.find("score") != std::string::npos && plan.find("_1_run") != std::string::npos, "the sweeps expand into run and score stages");
        }
        fs::remove_all(root);
        std::cout << "dynsyslab_figure_pipeline_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored; fs::remove_all(root, ignored);
        std::cerr << "dynsyslab_figure_pipeline_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
