// The ftle and ssm post-processing methods through executeJob: FTLE ridges of the double-well
// Duffing oscillator lie on its homoclinic loop H = v^2/2 - x^2/2 + x^4/4 = 0; the SSM of the
// hardening Duffing oscillator has the exact conservative backbone pi sqrt(1 + b a^2) / (2 K(m)),
// from the vector field and from trajectories (SSMLearn style).
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/EllipticIntegrals.hpp>
#include <dynsyslab/core/GridFieldArchive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
using O = json::Value::Object;
using A = json::Value::Array;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string duffingMeta(double alpha, double beta, double delta) {
    return R"({"schema":"dynsyslab.trajectory/v1","label":"source","system":"duffing","system_preset":"double_well_reference",)"
           R"("resolved_system_parameters":{"alpha":)" + std::to_string(alpha) + R"(,"beta":)" + std::to_string(beta)
         + R"(,"delta":)" + std::to_string(delta) + R"(,"input_gain":1},"integration_form":"dimensionless"})";
}

// An RK4 trajectory of x'' = -delta x' + alpha x - beta x^3 written as a source archive.
std::filesystem::path trajectory(const std::filesystem::path& path, double alpha, double beta, double delta,
                                 double x, double v, double dt, std::size_t samples) {
    std::vector<double> time, values;
    const auto f = [&](double a, double b, double& da, double& db) { da = b; db = -delta * b + alpha * a - beta * a * a * a; };
    for (std::size_t k = 0; k < samples; ++k) {
        time.push_back(double(k) * dt);
        values.push_back(x); values.push_back(v);
        for (int s = 0; s < 4; ++s) {
            const double h = dt / 4;
            double a1, b1, a2, b2, a3, b3, a4, b4;
            f(x, v, a1, b1); f(x + h / 2 * a1, v + h / 2 * b1, a2, b2);
            f(x + h / 2 * a2, v + h / 2 * b2, a3, b3); f(x + h * a3, v + h * b3, a4, b4);
            x += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4); v += h / 6 * (b1 + 2 * b2 + 2 * b3 + b4);
        }
    }
    core::NpzRecord source;
    source.arrays["time"] = core::NpyArray::fromDoubles({samples}, time);
    source.arrays["values"] = core::NpyArray::fromDoubles({samples, 2}, values);
    source.arrays["metadata_json"] = core::NpyArray::fromString(duffingMeta(alpha, beta, delta));
    core::writeNpzAtomic(path, source);
    return path;
}

dynsyslab::workflow::RunResult run(const std::vector<std::filesystem::path>& sources, const std::filesystem::path& root,
                                   const std::string& method, O settings) {
    A list;
    for (const auto& s : sources) list.push_back(s.generic_u8string());
    const json::Value job = O{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", method},
        {"sources", list}, {"methods", A{O{{"name", method}, {"enabled", true}, {"settings", std::move(settings)}}}},
        {"archive", O{{"enabled", true}}}};
    return dynsyslab::workflow::executeJob(job, root);
}

std::size_t output(const core::GridField& field, const std::string& name) {
    const auto at = std::find(field.outputNames.begin(), field.outputNames.end(), name);
    require(at != field.outputNames.end(), "field output " + name);
    return std::size_t(at - field.outputNames.begin());
}

void ftleTest(const std::filesystem::path& root) {
    const auto source = trajectory(root / "double-well.npz", 1, 1, 0, 0.5, 0, 0.05, 40);
    const auto result = run({source}, root / "ftle", "ftle", {{"domain", O{{"lower", A{-1.6, -1.0}}, {"upper", A{1.6, 1.0}},
        {"points", A{97.0, 61.0}}}}, {"integration_time", 6.0}, {"maximum_step", 0.01}});
    const auto archive = core::readNpz(result.archives.front());
    const auto series = core::readGridFieldArchive(archive);
    const auto& field = series.frames.front();
    require(series.metadata.at("schema").string() == "dynsyslab.grid-field/v1", "ftle publishes a grid field");
    require(archive.contains("lcs_points") && archive.contains("cauchy_green_eigenvalues_forward"), "LCS curves archived");
    // Ridges of the strongest stretching, and the strongest LCS of each direction, lie on the loop.
    const auto H = [](const double* q) { return q[1] * q[1] / 2 - q[0] * q[0] / 2 + std::pow(q[0], 4) / 4; };
    for (const char* direction : {"forward", "backward"}) {
        const auto sigma = output(field, std::string("ftle_") + direction), ridge = output(field, std::string("ridge_") + direction);
        double top = 0;
        for (std::size_t p = 0; p < field.grid.size(); ++p) top = std::max(top, field.value(p, sigma).real());
        std::size_t ridges = 0, onLoop = 0;
        for (std::size_t p = 0; p < field.grid.size(); ++p) {
            if (field.value(p, ridge).real() != 1 || field.value(p, sigma).real() < 0.7 * top) continue;
            double q[2];
            field.grid.point(p, q);
            if (std::abs(q[0]) < 0.3) continue;  // near the saddle the loop's two branches meet
            ++ridges;
            onLoop += std::abs(H(q)) < 0.03;
        }
        require(ridges > 20 && onLoop * 10 >= ridges * 8, std::string("FTLE ridges on the homoclinic loop, ") + direction);
    }
    const auto points = archive.at("lcs_points").doubles(), offsets = archive.at("lcs_offsets").doubles();
    const auto sign = archive.at("lcs_direction").doubles();
    for (const double want : {1.0, -1.0}) {
        const auto k = std::size_t(std::find(sign.begin(), sign.end(), want) - sign.begin());  // strongest first
        require(k < sign.size(), "an LCS in each direction");
        double worst = 0;
        for (auto s = std::size_t(offsets[k]); s < std::size_t(offsets[k + 1]); ++s) worst = std::max(worst, std::abs(H(&points[2 * s])));
        require(worst < 0.02 && offsets[k + 1] - offsets[k] > 50, "the strongest LCS is the homoclinic loop");
    }
    output(field, "lcs_repelling");
    output(field, "lcs_attracting");
}

double exactFrequency(double b, double a) {  // m = b a^2 / (2 (1 + b a^2)); Kovacic and Brennan 2011, section 2.2
    const double s = b * a * a;
    return 3.141592653589793 * std::sqrt(1 + s) / (2 * core::ellipticK(s / (2 * (1 + s))));
}

void checkBackbone(const dynsyslab::workflow::RunResult& result, double tolerance, const std::string& what) {
    require(result.archives.size() == 2, what + ": the model and its backbone");
    const auto model = core::readNpz(result.archives[0]);
    require(model.contains("dynamics_model_json") && model.contains("ssm_parameterization"), what + ": dynamical model arrays");
    const auto meta = json::parse(model.metadataJson());
    require(meta.at("schema").string() == "dynsyslab.spectral-submanifold/v1", what + ": schema");
    require(meta.at("data_object").at("kind").string() == "dynamical_model", what + ": a dynamical_model object");
    const auto backbone = core::readGridFieldArchive(core::readNpz(result.archives[1])).frames.front();
    const auto w = output(backbone, "frequency"), a = output(backbone, "amplitude_position");
    bool checked = false;
    for (std::size_t k = 1; k < backbone.grid.size(); ++k) {
        const double a0 = backbone.value(k - 1, a).real(), a1 = backbone.value(k, a).real();
        if (a0 <= 0.2 && 0.2 <= a1) {
            const double omega = backbone.value(k - 1, w).real()
                               + (backbone.value(k, w).real() - backbone.value(k - 1, w).real()) * (0.2 - a0) / (a1 - a0);
            require(std::abs(omega - exactFrequency(0.5, 0.2)) < tolerance, what + ": backbone frequency at a = 0.2");
            checked = true;
        }
    }
    require(checked, what + ": backbone reaches a = 0.2");
}

void ssmTest(const std::filesystem::path& root) {
    // From the system's own drift (hardening Duffing: alpha = -1, beta = 0.5, light damping).
    const auto source = trajectory(root / "hardening.npz", -1, 0.5, 0.004, 0.3, 0, 0.05, 200);
    checkBackbone(run({source}, root / "ssm", "ssm", {{"order", 7.0}, {"backbone", O{{"radius", 0.2}, {"points", 201.0}}}}),
                  3e-5, "parameterization SSM");
    // From trajectories: the reduced dynamics regressed on the data.
    std::vector<std::filesystem::path> sources;
    for (int k = 0; k < 3; ++k)
        sources.push_back(trajectory(root / ("decay-" + std::to_string(k) + ".npz"), -1, 0.5, 0.02,
                                     0.35 * std::cos(1.1 * k), 0.35 * std::sin(1.1 * k), 0.05, 3000));
    checkBackbone(run(sources, root / "ssm-data", "ssm", {{"mode", "data"}, {"dynamics_order", 5.0}, {"equilibrium", A{0.0, 0.0}},
                      {"backbone", O{{"points", 201.0}}}}), 1e-3, "data-driven SSM");
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-geometry-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        ftleTest(root);
        ssmTest(root);
        std::filesystem::remove_all(root);
        std::cout << "geometry post-processing tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
