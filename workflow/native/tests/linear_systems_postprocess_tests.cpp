// `linear_systems` workflow stage on a raw A/B/C/D archive: archive contents and metadata.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include "../src/BalancedTruncationSettings.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path() / ("dynsyslab-linear-systems-" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        core::NpzRecord matrices;  // stable lags -1, -5, -30 plus an unstable mode +0.2
        matrices.arrays["A"] = core::NpyArray::fromDoubles({4, 4}, {-1, 0, 0, 0, 0, -5, 0, 0, 0, 0, -30, 0, 0, 0, 0, 0.2});
        matrices.arrays["B"] = core::NpyArray::fromDoubles({4, 1}, {1, 1, 1, 1});
        matrices.arrays["C"] = core::NpyArray::fromDoubles({1, 4}, {1, 0.5, 0.1, 1});
        core::writeNpzAtomic(root / "matrices.npz", matrices);
        const json::Value job = json::Value::Object{
            {"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", "linear systems"},
            {"sources", json::Value::Array{(root / "matrices.npz").generic_u8string()}},
            {"methods", json::Value::Array{json::Value::Object{{"name", "linear_systems"}, {"enabled", true},
                {"settings", json::Value::Object{{"reduction_order", 2.0}, {"frequency_points", 50.0},
                    {"lqr_state_weight", 1.0}, {"lqr_input_weight", 0.1},
                    {"lqe_process_weight", 1.0}, {"lqe_measurement_weight", 0.01},
                    {"ltr_weights", json::Value::Array{0.0, 1e4}}, {"loop_frequency_points", 400.0}}}}}}};
        const auto result = dynsyslab::workflow::executeJob(job, root / "out");
        const auto record = core::readNpz(result.archives.front());
        for (const char* name : {"system_A", "system_poles", "hankel_singular_values", "hankel_error_bound",
                                 "reduced_A", "reduced_poles", "controllability_gramian", "observability_gramian",
                                 "controllability_gramian_eigenvalues", "frequencies", "full_singular_values",
                                 "reduced_singular_values", "full_phase", "transmission_zeros", "modal_damping",
                                 "lqr_gain", "reduced_lqr_lifted_poles", "lqe_gain", "reduced_lqe_lifted_poles"})
            require(record.contains(name), name);
        require(record.at("hankel_singular_values").shape[0] == 3, "stable part HSVs");
        require(record.at("reduced_A").shape[0] == 2, "reduced order includes the unstable mode");
        const auto meta = json::parse(record.metadataJson());
        require(meta.at("schema").string() == "dynsyslab.linear-systems/v1", "schema");
        require(meta.at("reduction").at("unstable_order").number() == 1, "unstable order");
        const double error = meta.at("reduction_error_hinf").number();
        require(error <= meta.at("reduction").at("error_bound").number() * (1 + 1e-6)
                && error >= meta.at("reduction").at("lower_bound").number() * (1 - 1e-6), "error within balanced bounds");
        require(meta.at("reduced_lqr_lifted_stable").boolean(), "reduced LQR lifted gain stabilizes the full model");
        // Coordinated output feedback (LQG): separation gives 2n poles; reduced compensator closes the full plant.
        require(record.at("lqg_closed_loop_poles").shape[0] == 8, "LQG closed loop has 2n poles");
        require(record.at("reduced_lqg_full_closed_loop_poles").shape[0] == 6, "reduced compensator on the full plant");
        require(record.at("reduced_lqg_compensator_A").shape[0] == 2, "reduced compensator order");
        require(meta.at("output_feedback").at("closed_loop_stable").boolean(), "full LQG closed loop stable");
        require(!meta.at("output_feedback").at("reduced_closed_loop_stable").isNull(), "reduced LQG stability reported");
        // Classical frequency domain: channel transfer function, Bode data and loop robustness.
        for (const char* name : {"transfer_function_poles", "transfer_function_denominator", "transfer_function_zeros__y0_u0",
                                 "transfer_function_numerator__y0_u0", "bode_frequencies", "bode_magnitude_db", "bode_phase_deg",
                                 "loop__lqr_input__response", "loop__kalman_output__return_difference",
                                 "loop__lqg_input__frequencies", "loop__lqg_output__response", "loop__ltr_1_input__response"})
            require(record.contains(name), name);
        require(record.at("transfer_function_denominator").shape[0] == 5, "denominator of a 4-state channel");
        require(record.at("bode_magnitude_db").shape[0] == 50, "Bode on the frequency grid");
        const auto& loops = meta.at("loop_robustness");
        require(meta.at("loop_names").array().size() == 6, "LQR, Kalman, LQG input/output and two LTR loops");
        require(loops.at("lqr_input").at("minimum_return_difference").number() >= 1 - 1e-6, "LQR guarantee recorded");
        require(loops.at("kalman_output").at("minimum_return_difference").number() >= 1 - 1e-6, "Kalman guarantee recorded");
        require(loops.at("lqg_input").at("closed_loop_stable").boolean(), "LQG input loop stable");
        require(loops.at("ltr_1_input").at("ltr_weight").number() == 1e4, "LTR weight recorded");
        require(meta.at("transfer_functions").array().size() == 1, "one SISO channel");
        // KLQE/KLQR stage option parser.
        using dynsyslab::workflow::parseBalancedTruncation;
        require(!parseBalancedTruncation(json::Value{}) && !parseBalancedTruncation(json::Value(false)), "absent option");
        require(parseBalancedTruncation(json::Value(true))->order == 0, "true selects the minimal order");
        require(parseBalancedTruncation(json::Value(3.0))->order == 3, "number selects the order");
        const auto parsed = parseBalancedTruncation(json::Value(json::Value::Object{
            {"error_tolerance", 0.05}, {"method", "residualization"}, {"stability_margin", 1e-6}}));
        require(parsed->order == 0 && parsed->errorTolerance == 0.05 && parsed->stabilityMargin == 1e-6
                && parsed->method == core::BalancedReductionMethod::residualization, "object form");
        bool rejected = false;
        try { (void)parseBalancedTruncation(json::Value(2.5)); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "fractional order rejected");
        std::filesystem::remove_all(root);
        std::cout << "dynsyslab_linear_systems_postprocess_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        if (!root.empty()) { std::error_code ignored; std::filesystem::remove_all(root, ignored); }
        std::cerr << "dynsyslab_linear_systems_postprocess_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
