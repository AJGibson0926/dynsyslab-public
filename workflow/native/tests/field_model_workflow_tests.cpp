// The field_model post-processing method through executeJob: a polynomial field of the damped
// single-well Duffing oscillator x'' = -0.1 x' - x - x^3 pooled over two trajectories, scored on a
// held-out trajectory by rate error and free rollouts, calibrated by fixed-support shooting; and a
// compact radial field of the same records.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FieldModels.hpp>
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

json::Value metadataOf(const std::filesystem::path& path) {
    return json::parse(core::readNpz(path).metadataJson());
}

void fieldModelTest(const std::filesystem::path& root) {
    const std::vector<std::filesystem::path> training = {trajectory(root / "a.npz", -1, 1, 0.1, 1.2, 0, 0.02, 600),
                                                         trajectory(root / "b.npz", -1, 1, 0.1, 0, 1.0, 0.02, 600)};
    const auto heldOut = trajectory(root / "c.npz", -1, 1, 0.1, -0.8, 0.5, 0.02, 600);
    const auto result = run(training, root / "polynomial", "field_model", {{"kind", "polynomial"},
        {"library", O{{"family", "polynomial"}, {"degree", 3.0}, {"threshold", 0.01}, {"ridge", 1e-10}}},
        {"validation_sources", A{heldOut.generic_u8string()}}, {"rollout", O{{"horizon", 2.0}, {"starts", 3.0}}},
        {"calibrate", O{{"duration", 1.0}, {"maximum_evaluations", 10.0}}},
        {"controls", O{{"restart_time", 4.0}, {"compare_from", 6.0}, {"visits", A{json::Value(A{1.0, 2.0})}}}}});
    const auto meta = metadataOf(result.archives.front());
    const auto& validation = meta.at("validation");
    require(validation.at("rate_relative_l2").number(1) < 1e-3, "held-out rate error of the polynomial field");
    require(validation.at("rollout_completion").number() == 1 && validation.at("rollout_relative_rms").number(1) < 1e-2,
            "held-out rollouts of the polynomial field");
    require(meta.at("calibration").at("shots").number() == 4, "two shooting starts per training record");
    const auto& control = meta.at("controls").at(heldOut.generic_u8string());
    require(control.at("step_change_over_record_rms").number(1) < 1e-6, "halving the step barely changes the rollout");
    require(!control.at("restart_on_common").isNull() && control.at("restart_on_common").array().at(0).number(1) < 1e-2, "a restart stays on the held-out record: " + json::dump(control));
    require(control.at("visits").array().at(0).at("distance").array().size() == 2, "visit separation per coordinate");
    const auto definition = json::parse(core::readNpz(result.archives.front()).at("field_model_json").string());
    const double x[2] = {0.5, -0.2};
    const auto rate = core::fieldFromDefinition(definition).evaluateReal(x);
    require(std::abs(rate[1] - (0.02 - 0.5 - 0.125)) < 1e-3, "the archived definition evaluates the field");

    const auto compact = run(training, root / "compact", "field_model", {{"kind", "compact_radial"},
        {"options", O{{"stride", 10.0}, {"neighbors", 6.0}, {"affine_tail", true}}},
        {"validation_sources", A{heldOut.generic_u8string()}}});
    require(metadataOf(compact.archives.front()).at("validation").at("rate_relative_l2").number(1) < 0.1, "compact radial field on held-out rates");
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-field-model-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        fieldModelTest(root);
        std::filesystem::remove_all(root);
        std::cout << "field model workflow tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
