// Bilinear K-ROMs in the estimator and controller, and the end-to-end input_map -> certificate ->
// KLQE/KLQR workflow.
#include "BilinearTestSupport.hpp"

#include <chrono>
#include <iostream>

namespace {
using namespace bilinear_test;

json::Value closedLoopJob(const fs::path& krom, const fs::path& certificate, double x0) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"certified-klqr",
      "system":{"name":"duffing","parameters":{"alpha":2,"beta":0.2,"delta":1}},
      "time":{"start":0,"stop":10,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","values":[[0,0]],"count":1,"clusters":1},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_evaluation":"sampled","control_sample_time":0.02},
      "controller":{"kind":"klqr","target":{"space":"state","kind":"constant","state":[0,0]},
        "settings":{"control_representation":"bilinear","Q":1,"R":1}},
      "archive":{"enabled":true}})");
    job.at("initial_conditions").at("values") = json::Value::Array{json::Value::Array{x0, 0.0}};
    auto& settings = job.at("controller").at("settings").object();
    settings["krom"] = krom.generic_u8string();
    settings["bilinear_certificate"] = certificate.generic_u8string();
    return job;
}

void estimatorAndControllerUseBilinear() {
    c::SystemConfiguration config = duffing(0.0);
    auto model = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {1., -2.});
    model.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
    model.setReconstruction("state", {{1. / 3, -1. / 3, 1. / 3, 2. / 3}, 2, 2, 0, 0, 0});
    c::BilinearControlMap map;
    map.modes = 2; map.inputs = 1; map.offsets = {1.0, 1.0}; map.operators = {0.4, 0.1, -0.2, 0.3};
    map.provenance = "analytic test map";
    model.setControlMap("bilinear", map);

    // KLQE: on a stable spectrum, negligible process noise makes the Kalman gain vanish, so the
    // prediction is exp((G + u N) dt) z + ZOH(b0) u. (An unstable mode would keep a nonzero gain.)
    auto stable = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {-1., -2.});
    stable.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
    stable.setReconstruction("state", {{1. / 3, -1. / 3, 1. / 3, 2. / 3}, 2, 2, 0, 0, 0});
    stable.setControlMap("bilinear", map);
    const auto design = c::designKoopmanLQEFromMeasurements(stable, "sensor", {1e-8, 0, 0, 1e-8}, {1e8});
    c::KoopmanLQEObserver observer(stable, design, "bilinear", {0.5, -0.25});
    (void)observer.observe(0.0, {0.25}, {0.7});
    const auto next = observer.observe(0.05, {0.25}, {0.7});
    auto generator = stable.generatorSplit();
    const auto split = c::continuousBilinearSplit(stable, "bilinear");
    for (std::size_t e = 0; e < 4; ++e) generator[e] += 0.7 * split.B[0][e];
    const auto held = dynsyslab::numerics::zeroOrderHold(generator, split.B0, 2, 1, 0.05);
    for (std::size_t r = 0; r < 2; ++r) {
        const double exact = held.transition[r * 2] * 0.5 + held.transition[r * 2 + 1] * -0.25 + held.response[r] * 0.7;
        check(std::abs(next.lifted[r] - exact) < 1e-6, "KLQE does not propagate the bilinear term");
    }

    // KLQR: gain is scheduled on B(z) = b0 + N z.
    auto shared = std::make_shared<const c::KoopmanModel>(model);
    const std::vector<double> Q{1, 0, 0, 1}, R{1};
    c::KoopmanLQRController scheduled(shared, Q, R, {}, "bilinear", c::KoopmanTargetSpace::realLiftedValues);
    const std::vector<Complex> lifted{0.5, -0.25};
    const double u = scheduled.valueFromLifted(0.0, lifted).front();
    const std::vector<double> B{1.0 + 0.4 * 0.5 + 0.1 * -0.25, 1.0 - 0.2 * 0.5 + 0.3 * -0.25};
    const auto solved = c::solveContinuousLQR(model.generatorSplit(), B, Q, R, 2, 1);
    check(std::abs(u + solved.gain[0] * 0.5 + solved.gain[1] * -0.25) < 1e-9, "KLQR ignores the bilinear map");
    check(scheduled.gainSolveCount() == 1, "KLQR gain solve count");

    // Certified gain: u = u* + K (z - z*), no Riccati solve, and it refuses augmentation.
    c::KoopmanLQRController certified(shared, Q, R, {}, "bilinear", c::KoopmanTargetSpace::realLiftedValues);
    certified.setCertifiedFeedback({{-2.0, 0.5}, {0.1, -0.05}, {0.2}, 0.02});
    const double v = certified.valueFromLifted(0.0, lifted).front();
    check(std::abs(v - (0.2 - 2.0 * (0.5 - 0.1) + 0.5 * (-0.25 + 0.05))) < 1e-14, "certified sign or target input");
    check(certified.gainSolveCount() == 0, "certified gain ran a Riccati solve");
    rejects([&] { certified.setBalancedReduction({}); }, "balanced truncation on a certified gain");
    c::KoopmanLQRController moved(shared, Q, R, [](double) { return std::vector<double>{0.9, 0.9}; },
                                  "bilinear", c::KoopmanTargetSpace::realLiftedValues);
    moved.setCertifiedFeedback({{-2.0, 0.5}, {0.1, -0.05}, {0.2}, 0.02});
    rejects([&] { (void)moved.valueFromLifted(0.0, lifted); }, "controller target that leaves the certificate");
}

void endToEnd(const fs::path& root) {
    const auto config = duffing(0.2);
    auto krom = c::KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {1., -2.});
    krom.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
    krom.setReconstruction("state", {{1. / 3, -1. / 3, 1. / 3, 2. / 3}, 2, 2, 0, 0, 0});
    const auto kromPath = root / "duffing-krom.npz";
    c::saveKoopmanModel(kromPath, krom);
    const auto inputOnly = forcedDuffing(root / "forced", 0.2, 1.0);
    const auto forced = forcedWithHolds(root / "forced-holds");
    check(forced.size() == 128, "forced records with holds");

    json::Value::Object mapSettings{{"krom", kromPath.generic_u8string()}, {"representation", "bilinear"}, {"fit", "held_input"}};
    const auto mapped = w::executePostprocess(postJob("input_map", mapSettings, forced), root).archives.front();
    const auto mappedBytes = bytes(mapped);
    check(c::loadKoopmanModel(mapped).generatorSplit() == krom.generatorSplit(), "input_map changed A");

    const auto certify = [&](json::Value::Object extra, std::vector<fs::path> records = {}) {
        if (records.empty()) records = forced;
        json::Value::Object s{{"krom", mapped.generic_u8string()}, {"target", json::Value::Object{{"state", json::Value::Array{0.0, 0.0}}}},
            {"input_bounds", json::Value::Array{0.75}}, {"state_radius", 0.5}, {"contraction", 0.99}, {"inflation", 1.25}, {"tolerance", 1e-5}};
        for (auto& [k, v] : extra) s[k] = std::move(v);
        return w::executePostprocess(postJob("bilinear_certificate", s, records), root).archives.front();
    };
    const auto path = certify({});
    check(bytes(mapped) == mappedBytes, "certificate modified the K-ROM");
    const auto archive = c::readNpz(path);
    const auto meta = json::parse(archive.metadataJson());
    const auto& record = meta.at("bilinear_certificate");
    std::cout << "certificate status=" << record.at("status").string()
              << " c_x=" << record.at("sector").at("c_x").number()
              << " c_u=" << record.at("sector").at("c_u").number()
              << " inscribed=" << record.at("inscribed_radius").number(0) << '\n';
    check(record.at("status").string() == "certified", "Duffing certificate failed: " + record.at("message").string());
    check(record.at("scope").string().find("not a guarantee for unseen states") != std::string::npos, "scope statement missing");
    const auto exact = dynsyslab::numerics::matrixExponential(krom.generatorSplit(), 2, 0.02);
    check(archive.at("frozen_A").doubles() == exact && archive.at("model_A").doubles() == exact, "certificate altered A");

    const auto heldOut = heldOutSectorCheck(archive, record, 0.02);
    std::cout << "held-out grid: points=" << heldOut.points << " violations=" << heldOut.violations
              << " worst r/bound=" << heldOut.worstRatio << " worst zero-input r/bound=" << heldOut.worstZeroInputRatio << '\n';
    check(record.at("sector").at("c_x").number() > 0.0 && record.at("sector").at("c_u").number() > 0.0, "sector is one-sided");
    check(heldOut.violations == 0, "the sector does not bound the true residual on the held-out grid");

    // An input box the records do not reach is refused rather than extrapolated.
    const auto wideBox = c::readNpz(certify({{"input_bounds", json::Value::Array{1.6}}}));
    check(json::parse(wideBox.metadataJson()).at("bilinear_certificate").at("status").string() == "sector_unsupported"
          && !wideBox.contains("certificate_K"), "input box beyond the data was certified");

    // Input-only forced records cannot identify c_x, so they yield no gain and cannot be deployed.
    const auto unsupportedPath = certify({}, {inputOnly});
    const auto unsupported = c::readNpz(unsupportedPath);
    const auto unsupportedMeta = json::parse(unsupported.metadataJson());
    check(unsupportedMeta.at("bilinear_certificate").at("status").string() == "sector_unsupported", "input-only data was not refused");
    check(!unsupportedMeta.at("bilinear_certificate").at("certified").boolean(true) && !unsupported.contains("certificate_K"), "unsupported archive carries a gain");
    auto unsupportedRun = closedLoopJob(mapped, unsupportedPath, 0.01);
    rejects([&] { (void)w::executeJob(unsupportedRun, root / "unsupported"); }, "KLQR deployed a sector_unsupported archive");

    // Closed loop from inside the certified ellipsoid.
    const double inscribed = record.at("inscribed_radius").number();
    const double x0 = 0.4 * inscribed / std::sqrt(5.0);
    const auto run = c::readNpz(w::executeJob(closedLoopJob(mapped, path, x0), root).archives.front());
    const auto x = run.at("values").doubles();
    check(std::hypot(x[x.size() - 2], x[x.size() - 1]) < 1e-3 * std::max(x0, 1e-3) + 1e-6, "certified KLQR did not settle");
    const auto runMeta = json::parse(run.metadataJson());
    const auto evidence = json::dump(runMeta);
    check(evidence.find("not a guarantee for unseen states") != std::string::npos, "run metadata lacks the certificate scope");
    const auto K = archive.at("certificate_K").doubles();
    const auto u = run.at("control_values").doubles();
    const auto z0 = krom.evaluateSplit(x.data());
    check(std::abs(u[0] - (K[0] * z0[0] + K[1] * z0[1])) < 1e-12, "applied input is not u = K z");

    // Output feedback: KLQE on position, same K-ROM and map.
    auto feedback = closedLoopJob(mapped, path, x0);
    feedback.at("controller").at("settings").object()["estimator"] = json::parse(
        R"({"measurement":"sensor","measurement_observable":"position","output":"state",
            "process_covariance":1,"measurement_covariance":0.001,"initial_lifted":[0,0]})");
    const auto estimated = c::readNpz(w::executeJob(feedback, root / "klqe").archives.front());
    const auto xe = estimated.at("values").doubles();
    check(std::hypot(xe[xe.size() - 2], xe[xe.size() - 1]) < 1e-3, "KLQE + certified gain did not settle");
    check(json::dump(json::parse(estimated.metadataJson())).find("not certified: the design assumes z is known") != std::string::npos,
          "output-feedback limit not reported");

    // Refused deployments.
    auto badClock = closedLoopJob(mapped, path, x0);
    badClock.at("solver").at("control_sample_time") = 0.01;
    rejects([&] { (void)w::executeJob(badClock, root / "bad-clock"); }, "wrong control clock");
    auto badCap = closedLoopJob(mapped, path, x0);
    badCap.at("controller").at("settings").object()["control_limit"] = 0.5;
    rejects([&] { (void)w::executeJob(badCap, root / "bad-cap"); }, "actuator cap inside the certified box");
    auto wrongKrom = closedLoopJob(kromPath, path, x0);
    rejects([&] { (void)w::executeJob(wrongKrom, root / "bad-krom"); }, "K-ROM without the certified map");

    // A target that is not an equilibrium yields an honest archive without a gain.
    const auto off = c::readNpz(certify({{"target", json::Value::Object{{"state", json::Value::Array{0.3, 0.0}}}}}));
    const auto offMeta = json::parse(off.metadataJson());
    check(offMeta.at("bilinear_certificate").at("status").string() == "target_not_equilibrium", "constant term not reported");
    check(!off.contains("certificate_K"), "uncertified archive carries a gain");
    auto uncertified = closedLoopJob(mapped, certify({{"target", json::Value::Object{{"state", json::Value::Array{0.3, 0.0}}}}}), x0);
    rejects([&] { (void)w::executeJob(uncertified, root / "uncertified"); }, "uncertified archive deployed");

    // Fitted operators keep A too.
    const auto fitted = c::readNpz(certify({{"operators", "fit"}}));
    check(fitted.at("model_A").doubles() == exact, "operator fit altered A");
}
}  // namespace

int main(int argc, char** argv) {
    const auto root = argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path()
        / ("dynsyslab-bilinear-workflow-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        fs::create_directories(root);
        estimatorAndControllerUseBilinear();
        endToEnd(root);
        std::cout << "bilinear workflow tests passed: " << checks << " checks\n";
        if (argc <= 1) { std::error_code ignored; fs::remove_all(root, ignored); }
        return 0;
    } catch (const std::exception& failure) {
        std::cerr << failure.what() << '\n';
        return 1;
    }
}
