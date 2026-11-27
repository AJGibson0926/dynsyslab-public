// Solve history features: the controller call log, partial history kept on an aborting
// failure, and the frozen dimensional-baseline RPE system.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace core = dynsyslab::core;
namespace json = dynsyslab::core::json;
namespace wf = dynsyslab::workflow;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value job(const std::string& system, const std::string& time, const std::string& controller,
                const std::string& initial, const std::string& solver, const std::string& history) {
    return json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"history",
        "system":)" + system + R"(,"time":)" + time + R"(,"controller":)" + controller
        + R"(,"initial_conditions":)" + initial + R"(,"solver":)" + solver + R"(,"history":)" + history
        + R"(,"campaign":{"enabled":false},"preflight":{"enabled":false}})");
}

const std::string duffing = R"({"name":"duffing","preset":"double_well_reference","parameters":{}})";
const std::string lqr = R"({"kind":"lqr","settings":{"A":[0,1,-1,-1],"B":[0,1],"Q":1,"R":1},"target":{"kind":"constant","state":[1,0]}})";

void controllerCallLog(const std::filesystem::path& root) {
    for (const std::string method : {"RK4", "RK45"}) {
        const auto run = wf::executeJob(job(duffing, R"({"start":0,"stop":2,"output_step":0.1})", lqr,
            R"({"mode":"explicit","count":1,"values":[[0.8,0]]})", R"({"method":")" + method + R"(","rtol":1e-3})",
            R"({"controller_call_log":true})"), root);
        const auto archive = core::readNpz(run.archives.at(0));
        const auto time = archive.at("controller_call_time").doubles();
        const auto kind = archive.at("controller_call_kind").doubles();
        const auto status = archive.at("controller_call_status").doubles();
        const auto states = archive.at("controller_call_state").doubles();
        require(time.size() > 21 && kind.size() == time.size() && states.size() == 2 * time.size(),
                "every controller call is logged with its state");
        std::set<double> kinds(kind.begin(), kind.end()), statuses(status.begin(), status.end());
        require(kinds.count(0.0) && kinds.count(2.0), "right-hand-side and output-sample calls are both logged");
        for (const auto value : statuses) require(value == 1.0 || value == 0.0 || value == -1.0, "status codes");
        if (method == "RK4") require(statuses == std::set<double>{1.0}, "a fixed-step solve accepts every call");
        const auto meta = json::parse(archive.at("metadata_json").string());
        require(meta.at("controller_call_log").at("entries").number() == static_cast<double>(time.size()),
                "the log is described in the metadata");
    }
}

void partialHistoryOnFailure(const std::filesystem::path& root) {
    const auto spec = job(duffing, R"({"start":0,"stop":5,"output_step":0.1})",
        R"({"kind":"none","settings":{}})", R"({"mode":"explicit","count":2,"values":[[0.5,0],[0.1,0]]})",
        R"({"method":"DOP853","on_failure":"error","state_bounds":[[-0.9,0.9],null]})",
        R"({"keep_partial_history_on_failure":true})");
    bool failed = false;
    try { (void)wf::executeJob(spec, root); }
    catch (const std::exception& error) {
        failed = std::string(error.what()).find("partial history kept in") != std::string::npos;
    }
    require(failed, "the solve still fails and names its partial history");
    std::filesystem::path partial;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root))
        if (entry.path().filename().string().find(".partial.npz") != std::string::npos) partial = entry.path();
    require(!partial.empty(), "a partial archive is kept");
    const auto archive = core::readNpz(partial);
    const auto meta = json::parse(archive.at("metadata_json").string());
    const auto rows = archive.at("time").doubles().size();
    require(meta.at("kept_on_failure").boolean() && rows > 1 && rows < 51
            && meta.at("completed_samples").number() == static_cast<double>(rows),
            "the partial archive holds the samples computed before the failure");
}

void dimensionalBaseline(const std::filesystem::path& root) {
    // The Python baseline right-hand side, evaluated independently (air-water 20 C).
    const double rho = 998.0, mu = 1.0e-3, sigma = 0.0725, pv = 2330.0, kappa = 4.0 / 3.0;
    const double r0 = 10.0e-6, pinf = 101325.0, pg0 = pinf - pv + 2.0 * sigma / r0;
    core::SystemConfiguration request;
    request.identifier = "rayleigh-plesset-dimensional-baseline";
    request.driver.amplitude = 40000.0;
    request.driver.frequencyHz = 2.0e5;
    request.driver.phaseRadians = 0.3;
    auto system = core::createSystem(request);
    require(system->descriptor().identifier == "rayleigh-plesset-dimensional-baseline"
            && system->descriptor().preferredCoordinates == "dimensional"
            && system->descriptor().stateNames == std::vector<std::string>({"radius", "radial_velocity"}),
            "baseline descriptor");
    const auto initial = system->initialState();
    require(std::abs(initial[0] - r0) < 1.0e-18 && initial[1] == 0.0, "baseline starts at rest at R0");
    for (const auto& [t, r, v] : std::vector<std::tuple<double, double, double>>{
             {0.0, 1.0e-5, 0.0}, {1.3e-6, 1.2e-5, 0.4}, {2.0e-6, 0.8e-5, -1.5}}) {
        const double drive = request.driver.amplitude * std::cos(2.0 * 3.141592653589793 * request.driver.frequencyHz * t
                                                                 + request.driver.phaseRadians);
        const double pg = pg0 * std::pow(r0 / r, 3.0 * kappa);
        const double jump = pg + pv - (pinf - drive) - 2.0 * sigma / r - 4.0 * mu * v / r;
        const double expected = (jump / rho - 1.5 * v * v) / r;
        const double state[2]{r, v};
        double derivative[2]{};
        system->rhs(t, state, derivative);
        require(derivative[0] == v && std::abs(derivative[1] - expected) <= 1.0e-9 * std::abs(expected),
                "baseline right-hand side equals the frozen Python form");
    }
    const auto run = wf::executeJob(job(R"({"name":"rayleigh-plesset-dimensional-baseline","parameters":{}})",
        R"({"start":0,"stop":2e-6,"output_step":1e-7})", R"({"kind":"none","settings":{}})",
        R"({"mode":"explicit","count":1,"values":[[1.05e-5,0]]})", R"({"method":"DOP853","rtol":1e-9,"atol":1e-15})",
        "{}"), root);
    const auto values = core::readNpz(run.archives.at(0)).at("values").doubles();
    require(std::abs(values[0] - 1.05e-5) < 1.0e-18 && std::isfinite(values.back()), "baseline solve job runs in SI units");
}

// A partitioned heat solve with rank-local halo exchange and adaptive rebalancing gives the
// serial numbers and records its communication, policy and repartitions.
void partitionedPdeJob(const std::filesystem::path& root) {
    const auto heat = [&](const std::string& partitioning) {
        return json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"heat",
            "system":{"name":"heat_equation_1d","family":"pde","preset":"periodic_diffusion","parameters":{"points":24})"
            + partitioning + R"(},"time":{"start":0,"stop":0.2,"output_step":0.05},"solver":{"method":"RK4"},
            "initial_conditions":{"mode":"default","count":1},"controller":{"kind":"none"},
            "campaign":{"enabled":false},"preflight":{"enabled":false}})");
    };
    const auto serial = wf::executeJob(heat(""), root / "serial");
    const auto parallel = wf::executeJob(heat(R"(,"pde_partitioning":{"enabled":true,"partition_shape":[3],
        "local_workers":2,"communication":"halo_exchange","adaptive":{"enabled":true,"rebalance_imbalance":1.0,
        "rebalance_interval_calls":3}})"), root / "parallel");
    const auto a = core::readNpz(serial.archives.at(0)).at("values").doubles();
    const auto b = core::readNpz(parallel.archives.at(0)).at("values").doubles();
    require(a.size() == b.size(), "partitioned PDE output grid");
    for (std::size_t index = 0; index < a.size(); ++index)
        require(std::abs(a[index] - b[index]) < 1.0e-14, "partitioned halo-exchange PDE equals the serial solve");
    const auto meta = json::parse(core::readNpz(parallel.archives.at(0)).at("metadata_json").string());
    require(meta.at("pde_communication").string() == "halo_exchange"
            && meta.at("pde_adaptive_policy").at("enabled").boolean()
            && !meta.at("pde_repartition_records").isNull(), "PDE parallel metadata");
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-solve-history-test";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    std::filesystem::create_directories(root);
    try {
        controllerCallLog(root / "log");
        partialHistoryOnFailure(root / "partial");
        dimensionalBaseline(root / "baseline");
        partitionedPdeJob(root / "pde");
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_solve_history_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_solve_history_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
