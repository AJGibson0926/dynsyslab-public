// Phase E pipeline: classical estimators, multi-input and SDC control, PCHIP/expression/lifted
// targets, the training-domain policy, mode rules, validation windows, the reconstruction gate,
// symmetric-pulse input maps, K-ROM tools, control-map path audits and training plans.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>
#include <dynsyslab/workflow/KromAuditQuery.hpp>
#include <dynsyslab/workflow/KromTrainingPlan.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;
using J = json::Value;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

J solve(double stop, J controller = J::Object{{"kind", "none"}}) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"phase-e",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":1,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,"values":[[-0.8,0.1]],"ranges":[]},
      "solver":{"method":"RK4"},"archive":{"enabled":true,"data_role":"free_training"}})");
    job.at("time").at("stop") = stop;
    job.at("controller") = controller;
    return job;
}

std::filesystem::path run(const J& job, const std::filesystem::path& root) {
    const auto result = dynsyslab::workflow::executeJob(job, root);
    require(!result.archives.empty(), "job published no archive");
    return result.archives.front();
}

std::filesystem::path postprocess(const std::string& method, const J& settings, const std::vector<std::filesystem::path>& sources,
                                  const std::filesystem::path& root) {
    J::Array list;
    for (const auto& path : sources) list.emplace_back(path.generic_u8string());
    const J job = J::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", "phase-e-" + method},
        {"sources", J(std::move(list))},
        {"methods", J::Array{J::Object{{"name", method}, {"enabled", true}, {"settings", settings}}}},
        {"archive", J::Object{{"enabled", true}}}};
    return run(job, root);
}

std::filesystem::path kromFixture(const std::filesystem::path& root) {
    auto fit = core::polynomialEigenfunctions(2, {1, 0, 0, 1}, {{1, 0}, {0, 0}, {0, 0}, {1, 0}},
        core::Spectrum::fromContinuousRates({{-1.0, 0.0}, {-2.0, 0.0}}), {0.0, 0.0}, {1.0, 1.0});
    core::SystemConfiguration system;
    system.identifier = "duffing"; system.preset = "double_well_reference";
    core::KoopmanModel model(std::move(fit), std::move(system));
    model.setControlMap("constant", core::ConstantControlMap{{{0.0, 0.0}, {1.0, 0.0}}, 2, 1, "fixture"});
    model.setReconstruction("state", {{{1, 0}, {0, 0}, {0, 0}, {1, 0}}, 2, 2, 0, 0, 0});
    std::map<std::string, core::NpyArray> extra;
    std::vector<double> states;
    for (int k = 0; k < 40; ++k) states.insert(states.end(), {-1.2 + 0.03 * k, -0.4 + 0.02 * k});
    extra["training_states"] = core::NpyArray::fromDoubles({40, 2}, states);
    const auto path = root / "fixture-krom.npz";
    core::saveKoopmanModel(path, model, R"({"schema":"dynsyslab.derived/v1","label":"fixture","delta_t":0.01})",
                           std::move(extra));
    return path;
}

void estimatorsAndControl(const std::filesystem::path& root) {
    auto job = solve(2.0);
    job.at("estimator") = json::parse(R"({"kind":"lqe","C":[1,0],"linearization_state":[-1,0],
        "process_covariance":1,"measurement_covariance":0.001,"sample_time":0.01,"initial_estimate":[-1,0]})");
    const auto archive = core::readNpz(run(job, root));
    require(archive.contains("state_estimate") && archive.contains("observer_poles"), "LQE traces were not archived");
    job.at("estimator") = json::parse(R"({"kind":"luenberger","C":[1,0],"L":[2,1],"linearization_state":[-1,0],"sample_time":0.01})");
    (void)run(job, root);
    job.at("estimator") = json::parse(R"({"kind":"static_map","C":[1,0],"M":[1,0],"sample_time":0.01})");
    (void)run(job, root);
    // Two inputs through an explicit B with additive actuation.
    const auto multi = core::readNpz(run(solve(5.0, json::parse(R"({"kind":"lqr","settings":{"actuation":"additive_state",
        "A":[0,1,2,0],"B":[1,0,0,1],"Q":[10,10],"R":[1,1],"control_limit":[5,5]},
        "target":{"kind":"constant","state":[-1,0]}})")), root));
    require(multi.at("control_values").shape.back() == 2, "multi-input control did not archive two inputs");
    const auto values = multi.at("values").doubles();
    require(std::abs(values[values.size() - 2] + 1.0) < 0.05, "multi-input LQR did not reach the target");
    // SDRE with a user factorisation (here constant, the Jacobian at the target).
    (void)run(solve(0.5, json::parse(R"({"kind":"sdre","settings":{"sdc_A":["0","1","2","0"],"Q":[10,10],"R":1},
        "target":{"kind":"constant","state":[-1,0]}})")), root);
}

void targets(const std::filesystem::path& root) {
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({4}, {0, 1, 2, 3});
    record.arrays["values"] = core::NpyArray::fromDoubles({4, 1}, {0, 1, 1, 2});
    const auto path = root / "target.npz";
    core::writeNpzAtomic(path, record);
    const auto pchip = dynsyslab::workflow::resolveControlTarget(J::Object{{"kind", "trajectory"},
        {"source", path.generic_u8string()}, {"interpolation", "pchip"}}, 1);
    require(pchip(1.0)[0] == 1.0 && std::abs(pchip(1.5)[0] - 1.0) < 1e-12, "PCHIP keeps a flat segment flat");
    const auto expression = dynsyslab::workflow::resolveControlTarget(J::Object{{"kind", "expression"},
        {"components", J::Array{"a*sin(t)", "a*cos(t)"}}, {"parameters", J::Object{{"a", 2.0}}}}, 2);
    require(std::abs(expression(0.5)[0] - 2 * std::sin(0.5)) < 1e-12, "expression target");
    const auto envelope = dynsyslab::workflow::controlTargetEnvelope(J::Object{{"kind", "rpe_canonical"},
        {"name", "rpe_multisine_offset_4_mixed"}}, 2, 0, 10);
    require(envelope.method == "analytic" && envelope.lower[0] < 4 && envelope.upper[0] > 4, "analytic envelope");
}

void kromPipeline(const std::filesystem::path& root) {
    const auto krom = kromFixture(root);
    const auto trajectory = run(solve(1.0), root);
    // KLQR with the domain policy (warn) and a lifted target.
    const auto klqr = core::readNpz(run(solve(0.3, J::Object{{"kind", "klqr"},
        {"settings", J::Object{{"krom", krom.generic_u8string()}, {"Q", 10.0}, {"R", 1.0}}},
        {"target", J::Object{{"kind", "constant"}, {"state", J::Array{-1.0, 0.0}}}}}), root));
    require(klqr.contains("domain_report_counts"), "training-domain report was not archived");
    (void)run(solve(0.3, J::Object{{"kind", "klqr"},
        {"settings", J::Object{{"krom", krom.generic_u8string()}, {"Q", 10.0}, {"R", 1.0}, {"target_domain_policy", "ignore"}}},
        {"target", J::Object{{"kind", "constant"}, {"space", "lifted"}, {"state", J::Array{-1.0, 0.0}}}}}), root);
    bool refused = false;
    try {
        (void)run(solve(0.3, J::Object{{"kind", "klqr"}, {"settings", J::Object{{"krom", krom.generic_u8string()},
            {"Q", 10.0}, {"R", 1.0}, {"target_domain_policy", "error"}}},
            {"target", J::Object{{"kind", "constant"}, {"state", J::Array{5.0, 0.0}}}}}), root);
    } catch (const std::exception&) { refused = true; }
    require(refused, "target_domain_policy error did not refuse a target outside the domain");
    // Multi-input KLQR and KMPC: two inputs through an explicit additive B, exact map grad psi B.
    for (const char* kind : {"klqr", "kmpc"}) {
        const auto multi = core::readNpz(run(solve(3.0, J::Object{{"kind", kind},
            {"settings", J::Object{{"krom", krom.generic_u8string()}, {"actuation", "additive_state"},
                {"B", J::Array{1.0, 0.0, 0.0, 1.0}}, {"Q", 100.0}, {"R", J::Array{1.0, 1.0}}, {"sample_time", 0.01}}},
            {"target", J::Object{{"kind", "constant"}, {"state", J::Array{-1.0, 0.0}}}}}), root));
        require(multi.at("control_values").shape.back() == 2, std::string(kind) + " did not archive two inputs");
        require(multi.contains("domain_report_counts"), std::string(kind) + " multi-input domain report");
        const auto values = multi.at("values").doubles();
        require(std::abs(values[values.size() - 2] + 1.0) < 0.05, std::string(kind) + " multi-input did not reach the target");
    }
    // Mode rules.
    for (const auto& rule : {J(J::Object{{"rule", "nearest_frequency"}, {"target_rate", J::Array{-1.0, 0.0}}}),
                             J(J::Object{{"rule", "validation_score"}, {"validation_sources", J::Array{trajectory.generic_u8string()}}}),
                             J(J::Object{{"rule", "state_observability"}, {"max_modes", 2.0}})}) {
        const auto extracted = core::readNpz(postprocess("krom_extract", J::Object{{"mode_selection", rule}}, {krom}, root));
        require(!json::parse(extracted.metadataJson()).at("mode_selection_report").isNull(), "mode rule report missing");
    }
    // Validation windows.
    const auto validation = core::readNpz(postprocess("koopman_validation", J::Object{{"krom", krom.generic_u8string()},
        {"windows", J::Array{J::Array{0.0, 0.5}, J::Array{0.5, J()}}}, {"binned_defect", J::Object{{"axis", 0.0}, {"bins", 4.0}}},
        {"history", true}}, {trajectory}, root));
    require(validation.contains("window_anchored_relative_l2") && validation.contains("psi_history"), "validation windows missing");
    const auto cells = core::readNpz(postprocess("koopman_validation", J::Object{{"krom", krom.generic_u8string()},
        {"segments", J::Object{{"starts", J::Array{0.0, 0.4, 0.8}}, {"horizon", 0.05}}},
        {"binned_defect", J::Object{{"axes", J::Array{0.0, 1.0}}, {"bins", J::Array{3.0, 3.0}}, {"horizon", 0.05}}}},
        {trajectory}, root));
    require(cells.contains("segment_horizon_relative_l2") && cells.contains("cell_log10_defect")
            && cells.contains("cell_defect_occupied"), "segment errors or cell defect missing");
    // Reconstruction gate: psi = x reconstructs the state exactly.
    const auto gated = core::readNpz(postprocess("observable_reconstruction", J::Object{{"krom", krom.generic_u8string()},
        {"observable", "state"}, {"quality", J::Object{{"validation_fraction", 0.25}}}}, {trajectory}, root));
    require(json::parse(gated.metadataJson()).at("validated").boolean(), "exact readout failed the gate");
    refused = false;
    try { (void)postprocess("observable_reconstruction", J::Object{{"krom", krom.generic_u8string()}, {"observable", "state"},
              {"quality", J::Object{{"maximum_condition_number", 0.5}}}}, {trajectory}, root); }
    catch (const std::exception&) { refused = true; }
    require(refused, "a failing reconstruction was published without allow_quality_failure");
    // Symmetric-pulse input map from paired one-sample runs.
    auto pulse = solve(0.01, J::Object{{"kind", "symmetric_pulse"}, {"settings", J::Object{{"magnitude", 0.5}}}});
    pulse.at("initial_conditions") = json::parse(R"({"mode":"explicit","count":2,"clusters":1,"seed":0,
        "values":[[-0.8,0.1],[-0.8,0.1]],"ranges":[]})");
    pulse.at("archive").at("data_role") = "forced_training";
    const auto pulses = dynsyslab::workflow::executeJob(pulse, root);
    require(pulses.archives.size() == 2, "a pulse pair needs two archives");
    (void)postprocess("input_map", J::Object{{"krom", krom.generic_u8string()}, {"representation", "constant"},
        {"fit", "symmetric_pulse"}, {"map_name", "pulse"}}, pulses.archives, root);
    // K-ROM tools.
    (void)postprocess("control_error_compare", J::Object{{"candidates", J::Array{J::Object{{"krom", krom.generic_u8string()},
        {"map", "constant"}}}}, {"grid", J::Object{{"axis0", J::Array{-1.0, 1.0}}, {"axis1", J::Array{-1.0, 1.0}},
        {"columns", 5.0}, {"rows", 5.0}}}}, {trajectory}, root);
    const auto rollout = core::readNpz(postprocess("krom_rollout", J::Object{{"krom", krom.generic_u8string()}}, {trajectory}, root));
    require(rollout.contains("values"), "krom_rollout did not decode the state");
    (void)postprocess("pod_galerkin", J::Object{{"rank", 2.0}}, {trajectory}, root);
    (void)postprocess("krom_refit", J::Object{{"krom", krom.generic_u8string()}, {"prune_threshold", 1e-6}}, {trajectory}, root);
    {   // gauge_reference to the K-ROM itself, without a refit: every gauge factor is 1.
        const auto aligned = postprocess("krom_refit", J::Object{{"krom", krom.generic_u8string()}, {"refit", false},
            {"gauge_reference", krom.generic_u8string()}}, {trajectory}, root);
        const auto meta = json::parse(dynsyslab::core::readNpz(aligned).metadataJson());
        bool unit = false;
        for (const auto* node : {&meta, &meta.at("krom_tool")})
            for (const auto& factor : node->at("gauge_factors").array()) {
                unit = true;
                require(std::abs(factor.array()[0].number() - 1) < 1e-9 && std::abs(factor.array()[1].number()) < 1e-9, "gauge_reference self-alignment");
            }
        require(unit, "gauge factors reported");
    }
    {   // orient + coefficient_summary: ranked terms per mode, largest first.
        const auto oriented = postprocess("krom_refit", J::Object{{"krom", krom.generic_u8string()}, {"refit", false},
            {"orient", true}, {"coefficient_summary", 1e-3}}, {trajectory}, root);
        const auto meta = json::parse(dynsyslab::core::readNpz(oriented).metadataJson());
        const auto& node = meta.at("coefficient_summary").isNull() ? meta.at("krom_tool") : meta;
        const auto& modes = node.at("coefficient_summary").array();
        require(!modes.empty() && !modes.front().array().empty(), "coefficient summary reported");
        const auto& rows = modes.front().array();
        for (std::size_t k = 1; k < rows.size(); ++k)
            require(rows[k].at("magnitude").number() <= rows[k - 1].at("magnitude").number(), "terms ranked by magnitude");
    }
    (void)postprocess("reconstruction_transfer", J::Object{{"krom", krom.generic_u8string()},
        {"source_krom", krom.generic_u8string()}}, {trajectory}, root);
    // Path diagnostics of the stored map.
    const auto audit = dynsyslab::workflow::evaluateKromAuditQuery(J::Object{{"schema", "dynsyslab.krom-audit-query/v1"},
        {"action", "control_map_path"}, {"krom", krom.generic_u8string()}, {"input_map", "constant"},
        {"trajectory", trajectory.generic_u8string()}});
    require(audit.at("stable_fraction").number() >= 0.0 && !audit.at("relative_error").array().empty(), "control_map_path audit");
    // Training plan from a target envelope and a domain.
    const auto plan = dynsyslab::workflow::designKromTrainingPlan(J::Object{{"krom", krom.generic_u8string()},
        {"requirements", J::Array{J::Object{{"name", "box"}, {"domain", J::Object{{"minimum", J::Array{-1.5, -0.6}},
            {"maximum", J::Array{0.0, 0.6}}}}},
            J::Object{{"name", "orbit"}, {"target", J::Object{{"kind", "expression"},
                {"components", J::Array{"-1+0.1*sin(t)", "0.1*cos(t)"}}}}, {"probes", 8.0}}}},
        {"points_per_axis", J::Array{3.0, 3.0}}, {"time", J::Object{{"stop", 0.1}, {"output_step", 0.01}}}});
    require(plan.at("probes").array().size() == 8 && !plan.at("plan_fingerprint_sha256").string().empty(), "training plan");
    (void)dynsyslab::workflow::planJob(plan.at("solve_job"));
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / ("dynsyslab-phase-e-" + std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::error_code ignored;
    try {
        std::filesystem::create_directories(root);
        estimatorsAndControl(root); targets(root); kromPipeline(root);
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_control_estimation_workflow_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_control_estimation_workflow_tests: FAIL: " << error.what() << '\n';
        std::filesystem::remove_all(root, ignored);
        return 1;
    }
}
