// State-space controllers (LQR, SDRE, MPC, PID) through the native solve, and rejection of a
// state-dependent classical controller without a defined control clock.
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}  // namespace

void verifyStateControllers(const std::filesystem::path& root) {
    const auto verifyStateController = [&](const std::string& kind,
                                           const std::string& settings,
                                           const std::string& evaluation) {
        const std::string controllerJson =
            "{\"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"solve\","
            "\"label\":\"native-" + kind + "-test\","
            "\"system\":{\"name\":\"duffing\",\"preset\":\"double_well_reference\",\"parameters\":{}},"
            "\"time\":{\"start\":0,\"stop\":0.02,\"output_step\":0.005,\"coordinate_system\":\"automatic\"},"
            "\"initial_conditions\":{\"mode\":\"explicit\",\"count\":1,\"clusters\":1,\"seed\":0,"
            "\"values\":[[-1,0.05]],\"ranges\":[]},"
            "\"controller\":{\"kind\":\"" + kind + "\",\"settings\":"
            + settings + ",\"target\":{\"kind\":\"constant\",\"state\":[-1,0]}},"
            "\"solver\":{\"method\":\"RK4\",\"step_size\":0.001,"
            "\"control_evaluation\":\"automatic\"},\"archive\":{\"enabled\":true}}";
        const auto run = dynsyslab::workflow::executeSolve(
            dynsyslab::core::json::parse(controllerJson), root,
            root / (kind + "-progress.json"));
        require(run.archives.size() == 1, "state controller archive");
        const auto record = dynsyslab::core::readNpz(run.archives.front());
        const auto metadata = record.metadataJson();
        require(metadata.find("\"controller_family\":\"state_space\"")
                != std::string::npos, "state controller family metadata");
        require(metadata.find("\"control_evaluation\":\"" + evaluation + "\"")
                != std::string::npos, "state controller clock metadata");
        const auto controls = record.at("control_values").doubles();
        require(!controls.empty(), "state controller history");
        for (const double value : controls)
            require(std::isfinite(value), "state controller finite history");
    };
    verifyStateController("lqr",
        R"({"Q":[10,1],"R":1,"linearization_state":[-1,0]})",
        "continuous_rhs");
    verifyStateController("sdre",
        R"({"Q":[10,1],"R":1,"sample_time":0.005})", "sampled");
    verifyStateController("mpc",
        R"({"Q":[10,1],"R":1,"horizon":4,"sample_time":0.005})",
        "sampled");
    verifyStateController("pid",
        R"({"Kp":[2,0.5],"Ki":[0.1,0],"Kd":[0,0],"sample_time":0.005})",
        "sampled");

    auto unsafeSdre = dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"unsafe-sdre",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":0.01,"output_step":0.005,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,
                            "values":[[-1,0]],"ranges":[]},
      "controller":{"kind":"sdre","settings":{"Q":1,"R":1},
                    "target":{"kind":"constant","state":[-1,0]}},
      "solver":{"method":"RK4","step_size":0.001,
                "control_evaluation":"continuous_rhs"},
      "archive":{"enabled":true}
    })");
    bool unsafeClockRejected = false;
    try {
        (void)dynsyslab::workflow::executeSolve(
            unsafeSdre, root, root / "unsafe-sdre-progress.json");
    } catch (const std::invalid_argument&) {
        unsafeClockRejected = true;
    }
    require(unsafeClockRejected, "state-dependent classical control requires a defined clock");
}
