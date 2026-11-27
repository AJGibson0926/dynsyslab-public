// Solve-job features of the phase-B native port: designed inputs (tabulated, command tape,
// telegraph-plus-multisine, symmetric pulse pairs), IC modes (stratified, target tube,
// concatenate, repeat, clip), guards with rejection, post-solve resampling, checkpoints with
// pause and resume, continuation, dimensional integration, preset drivers, Marmottant shell
// laws, the new ODE methods and the adaptive method of lines.

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/Solve.hpp>

#include "../src/SolveContinuation.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace json = dynsyslab::core::json;
namespace wf = dynsyslab::workflow;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value job(const std::string& system, const std::string& time, const std::string& controller,
                const std::string& initial, const std::string& solver = R"({"method":"RK4"})",
                const std::string& extra = "") {
    return json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"features",
        "system":)" + system + R"(,"time":)" + time + R"(,"controller":)" + controller
        + R"(,"initial_conditions":)" + initial + R"(,"solver":)" + solver
        + R"(,"campaign":{"enabled":false},"preflight":{"enabled":false})" + extra + "}");
}

const std::string duffing = R"({"name":"duffing","preset":"double_well_reference","parameters":{}})";
const std::string shortTime = R"({"start":0,"stop":0.5,"output_step":0.1})";
const std::string none = R"({"kind":"none","settings":{}})";

json::Value metadata(const std::filesystem::path& path) {
    return json::parse(dynsyslab::core::readNpz(path).at("metadata_json").string());
}

std::vector<double> array(const std::filesystem::path& path, const char* name) {
    return dynsyslab::core::readNpz(path).at(name).doubles();
}

void testDesignedInputs(const std::filesystem::path& root) {
    auto run = wf::executeJob(job(duffing, shortTime,
        R"({"kind":"tabulated","settings":{"times":[0,0.5],"values":[0,1]}})",
        R"({"mode":"default","count":1})"), root);
    auto control = array(run.archives[0], "control_values");
    require(std::abs(control[2] - 0.4) < 1.0e-12, "tabulated input is piecewise linear");

    // MATLAB Ubig layout: 5 transitions x 2 trajectories, source scale 0.5.
    run = wf::executeJob(job(duffing, shortTime,
        R"({"kind":"command_tape","settings":{"interval":0.1,"layout":"transition_major","scale":0.5,
            "values":[[1,-1],[2,-2],[3,-3],[4,-4],[5,-5]]}})",
        R"({"mode":"default","count":2})"), root);
    control = array(run.archives[1], "control_values");
    require(control[0] == -0.5 && control[2] == -1.5, "command tape is interval-held, scaled, transposed");
    const auto design = metadata(run.archives[0]).at("input_design");
    require(design.at("kind").string() == "command_tape" && design.at("tape_sha256").string().size() == 64,
            "command tape digest");

    run = wf::executeJob(job(duffing, R"({"start":0,"stop":2,"output_step":0.01})",
        R"({"kind":"telegraph_multisine","settings":{"seed":742391,"interval":0.01,"transition_count":200,
            "amplitudes":[0.45,0.6],"hold_steps":[10,20],"multisine_fraction":0.28,
            "angular_frequencies":[1.5,3,4.8,7]}})",
        R"({"mode":"default","count":3})"), root);
    for (std::size_t member = 0; member < 3; ++member) {
        control = array(run.archives[member], "control_values");
        const auto meta = metadata(run.archives[member]).at("input_design");
        const double amplitude = meta.at("amplitude").number();
        double mean = 0.0, peak = 0.0;
        for (std::size_t k = 0; k < 200; ++k) { mean += control[k]; peak = std::max(peak, std::abs(control[k])); }
        require(std::abs(mean / 200.0) < 1.0e-12, "telegraph-multisine tape is zero mean");
        require(peak <= amplitude + 1.0e-12, "telegraph-multisine tape is peak clipped");
        require(meta.at("hold_steps").number() == (member < 2 ? 10.0 : 20.0), "per-member hold steps");
    }
}

void testPairsAndInitialConditions(const std::filesystem::path& root) {
    const auto rpe = R"({"name":"rpe","preset":"air_water_20c","parameters":{}})";
    const auto initial = R"({"mode":"concatenate","count":12,"repeat":2,"clip":[[0.72,1.28],null],"parts":[
        {"mode":"target_tube","count":4,"seed":3,"offset":1,"amplitudes":[0.09,0.12],"angular_frequencies":[4,3],
         "phases":[0,3.141592653589793],"tau_span":[0,20],"jitter":[0.025,0.08]},
        {"mode":"uniform","count":2,"seed":4,"ranges":[[0.76,1.24],[-0.9,0.9]]}]})";
    const auto spec = job(rpe, R"({"start":0,"stop":0.002,"output_step":0.001})",
        R"({"kind":"symmetric_pulse","settings":{"magnitude":1}})", initial);
    const auto set = wf::resolveOdeInitialConditions(spec);
    require(set.trajectoryCount == 12, "pair IC count");
    for (std::size_t pair = 0; pair < 6; ++pair)
        require(set.values[4 * pair] == set.values[4 * pair + 2] && set.values[4 * pair + 1] == set.values[4 * pair + 3],
                "repeat gives both pulse polarities the same anchor");
    for (std::size_t row = 0; row < 12; ++row)
        require(set.values[2 * row] >= 0.72 && set.values[2 * row] <= 1.28, "IC clip");
    const auto run = wf::executeJob(spec, root);
    const auto minus = metadata(run.archives[4]).at("input_design"), plus = metadata(run.archives[5]).at("input_design");
    require(minus.at("pair_id").number() == 2.0 && plus.at("pair_id").number() == 2.0, "pair_id");
    require(minus.at("polarity").number() == -1.0 && plus.at("polarity").number() == 1.0, "polarity");

    const auto stratified = job(duffing, shortTime, none,
        R"({"mode":"stratified","count":10,"seed":7,"ranges":[[0,1],[-1,1]]})");
    const auto strata = wf::resolveOdeInitialConditions(stratified);
    std::set<int> used;
    for (std::size_t row = 0; row < 10; ++row) used.insert(static_cast<int>(std::floor(strata.values[2 * row] * 10.0)));
    require(used.size() == 10, "stratified IC uses every stratum once");
}

void testRejectAndResample(const std::filesystem::path& root) {
    auto spec = job(duffing, shortTime, none,
        R"({"mode":"explicit","count":3,"values":[[0.1,0],[0.5,0],[0.2,0]]})",
        R"({"method":"DOP853","on_failure":"reject","state_bounds":[[-0.3,0.3],null]})");
    auto run = wf::executeJob(spec, root);
    require(run.archives.size() == 2, "reject drops the trajectory that leaves the bounds");
    std::filesystem::path manifest;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root))
        if (entry.path().filename().string().find(".outcome.json") != std::string::npos) manifest = entry.path();
    require(!manifest.empty(), "solve outcome record");
    std::ifstream stream(manifest);
    const auto outcome = json::parse(std::string((std::istreambuf_iterator<char>(stream)), {}));
    require(outcome.at("accepted_trajectory_count").number() == 2.0
            && outcome.at("rejected").array().front().at("trajectory_index").number() == 1.0, "rejected member recorded");

    spec = job(duffing, R"({"start":0,"stop":0.5,"output_step":0.1,"sample_times":[0,0.05,0.12,0.5]})", none,
               R"({"mode":"default","count":1})", R"({"method":"Radau"})", R"(,"interpolation":{"method":"pchip"})");
    run = wf::executeJob(spec, root);
    const auto time = array(run.archives[0], "time");
    require(time.size() == 4 && time[2] == 0.12, "resampled output times");
    require(array(run.archives[0], "raw_time").size() == 6, "raw adaptive series kept");
    const auto meta = metadata(run.archives[0]);
    require(meta.at("interpolation").at("method").string() == "pchip"
            && meta.at("integration_method").string() == "Radau", "resampling and method metadata");
}

// coordinate_transform kind "resample": a non-uniform solve output onto a uniform PCHIP clock.
void testCoordinateResample(const std::filesystem::path& root) {
    const auto solved = wf::executeJob(job(duffing, R"({"start":0,"stop":1,"output_times":[0,0.1,0.15,0.4,0.7,1]})", none,
        R"({"mode":"explicit","count":1,"values":[[0.1,0]]})", R"({"method":"DOP853"})"), root);
    const auto post = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"uniform",
        "sources":[")" + solved.archives[0].generic_string() + R"("],
        "methods":[{"name":"coordinate_transform","enabled":true,
                    "settings":{"kind":"resample","time_step":0.25,"method":"pchip","store_values":true}}],
        "archive":{"enabled":true}})");
    const auto published = wf::executeJob(post, root).archives.front();
    const auto time = array(published, "time");
    require(time.size() == 5 && std::abs(time[1] - 0.25) < 1.0e-15, "resampled uniform clock");
    const auto meta = metadata(published);
    require(meta.at("resampling").at("method").string() == "pchip" && meta.at("uniform_time_steps").boolean(),
            "resampling provenance");
    const auto values = array(published, "values"), raw = array(solved.archives[0], "values");
    require(std::abs(values[0] - raw[0]) < 1.0e-15 && std::abs(values[8] - raw[10]) < 1.0e-15, "resampling keeps the end samples");
}

// Delay and lattice solves take history.checkpoint_interval; partial archives are removed at the end.
void testDelayLatticeCheckpoints(const std::filesystem::path& root) {
    auto delay = wf::executeJob(json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"mg",
        "system":{"name":"mackey_glass","family":"dde","preset":"classic_chaos","parameters":{"history_value":1.2}},
        "time":{"start":0,"stop":20,"output_step":0.1},"solver":{"step_size":0.1},"history":{"checkpoint_interval":5},
        "initial_conditions":{"mode":"explicit","count":1,"values":[[1.2]]},"controller":{"kind":"none"},
        "campaign":{"enabled":false},"preflight":{"enabled":false}})"), root);
    require(array(delay.archives[0], "time").size() == 201, "checkpointed delay solve completes");
    auto lattice = wf::executeJob(json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"ising",
        "system":{"name":"ising","family":"lattice","preset":"critical","parameters":{"size":4,"temperature":2.269}},
        "time":{"start":0,"stop":12,"output_step":1},"solver":{"method":"metropolis"},"history":{"checkpoint_interval":5},
        "initial_conditions":{"mode":"default","count":1,"seed":4},"controller":{"kind":"none"},
        "campaign":{"enabled":false},"preflight":{"enabled":false}})"), root);
    require(array(lattice.archives[0], "time").size() == 12, "checkpointed lattice solve completes");
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root))
        require(entry.path().filename().string().find(".partial.npz") == std::string::npos, "partial archives removed");
}

void testPauseResumeContinue(const std::filesystem::path& root) {
    std::filesystem::create_directories(root);
    const auto pause = root / "pause-now";
    { std::ofstream touch(pause); touch << "pause"; }
    auto spec = job(duffing, R"({"start":0,"stop":1,"output_step":0.1})", none,
        R"({"mode":"explicit","count":2,"values":[[0.1,0],[0.2,0]]})", R"({"method":"RK4"})",
        R"(,"history":{"checkpoint_interval":0.3,"pause_file":")" + pause.generic_string() + R"("})");
    const auto run = wf::executeJob(spec, root);
    require(run.archives.size() == 1, "pause archives only the partial trajectory");
    const auto paused = metadata(run.archives[0]);
    require(paused.at("integration_paused").at("completed_samples").number() == 4.0, "paused at the first checkpoint");
    std::filesystem::path outcomePath;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root))
        if (entry.path().filename().string().find(".outcome.json") != std::string::npos) outcomePath = entry.path();
    std::ifstream outcomeStream(outcomePath);
    const auto outcome = json::parse(std::string((std::istreambuf_iterator<char>(outcomeStream)), {}));
    require(outcome.at("paused").boolean() && outcome.at("resume_jobs").array().size() == 2
            && outcome.at("not_started").array().size() == 1, "pause outcome lists the resume jobs");
    std::filesystem::remove(pause);
    // Resume: continue the paused archive and solve the member that had not started.
    const auto extension = wf::detail::extensionJob(run.archives[0], 1.0);
    const auto continued = wf::executeJob(extension, root);
    const auto child = metadata(continued.archives[0]);
    require(std::abs(child.at("continuation").at("parent_time").number() - 0.3) < 1.0e-12
            && child.at("continuation").at("parent_archive_sha256").string().size() == 64, "child archive links its parent");
    const auto parentValues = array(run.archives[0], "values");
    const auto childValues = array(continued.archives[0], "values");
    require(childValues[0] == parentValues[6] && childValues[1] == parentValues[7], "continuation starts at the parent's last state");
    // A free run over [0, 1] reaches the same final state (fixed RK4 steps are restart-exact).
    auto full = job(duffing, R"({"start":0,"stop":1,"output_step":0.1})", none,
                    R"({"mode":"explicit","count":1,"values":[[0.1,0]]})");
    const auto whole = array(wf::executeJob(full, root).archives[0], "values");
    require(std::abs(whole.back() - childValues.back()) < 1.0e-14, "pause and resume reproduce the uninterrupted solve");
}

void testCoordinatesPresetsShellsPde(const std::filesystem::path& root) {
    auto run = wf::executeJob(job(R"({"name":"rpe","preset":"air_water_20c","parameters":{}})",
        R"({"start":0,"stop":1e-6,"output_step":1e-7,"coordinate_system":"dimensional"})", none,
        R"({"mode":"explicit","count":1,"values":[[1.05e-5,0]]})",
        R"({"method":"DOP853","rtol":1e-9,"atol":1e-15})"), root);
    const auto meta = metadata(run.archives[0]);
    require(meta.at("integration_form").string() == "dimensional", "dimensional RPE archive");

    run = wf::executeJob(job(R"({"name":"nonlinear_pendulum","preset":"driven_damped","parameters":{},"driver":"preset"})",
        shortTime, none, R"({"mode":"default","count":1})"), root);
    require(metadata(run.archives[0]).at("system_driver").at("amplitude").number() == 0.7, "preset driver applied");
    run = wf::executeJob(job(R"({"name":"duffing","preset":"gibson_damped_stabilization_2022","parameters":{}})",
        shortTime, none, R"({"mode":"default","count":1})"), root);
    require(metadata(run.archives[0]).at("preset_reference").string().find("10.1137/21M1456078") != std::string::npos,
            "preset reference archived");

    for (const char* law : {"1", "2"}) {
        run = wf::executeJob(job(std::string(R"({"name":"marmottant","preset":"matlab_cleanup_linear","parameters":{"shell_law":)")
            + law + "}}", R"({"start":0,"stop":1,"output_step":0.1})", none, R"({"mode":"default","count":1})",
            R"({"method":"LSODA"})"), root);
        require(std::isfinite(array(run.archives[0], "values").back()), "Marmottant smoothed shell law solves");
    }

    const auto heat = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"heat",
        "system":{"name":"heat_equation_1d","family":"pde","preset":"periodic_diffusion","parameters":{"points":16}},
        "time":{"start":0,"stop":0.5,"output_step":0.1},"solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},
        "initial_conditions":{"mode":"default","count":1},"controller":{"kind":"none"},
        "campaign":{"enabled":false},"preflight":{"enabled":false}})");
    run = wf::executeJob(heat, root);
    require(metadata(run.archives[0]).at("pde_method").string() == "MOL-RK45", "adaptive method of lines");
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-solve-features-test";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    try {
        testDesignedInputs(root);
        testPairsAndInitialConditions(root);
        testRejectAndResample(root / "reject");
        testCoordinateResample(root / "resample");
        testDelayLatticeCheckpoints(root / "checkpoints");
        testPauseResumeContinue(root / "pause");
        testCoordinatesPresetsShellsPde(root);
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_solve_features_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_solve_features_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
