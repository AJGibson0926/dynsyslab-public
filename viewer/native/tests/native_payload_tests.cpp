#include "../src/ModeWeightPayload.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ComparisonPayload.hpp"
#include "dynsyslab_viewer/RecordWindow.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void trajectoryViews() {
    const auto path = std::filesystem::temp_directory_path() / "dynsyslab-native-payload-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({3}, {0.0, 0.5, 1.0});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles({3, 2},
        {1.0, 0.0, 1.1, 0.2, 1.2, 0.1});
    archive.arrays["control_values"] = dynsyslab::core::NpyArray::fromDoubles({3, 1},
        {0.0, 0.3, 0.1});
    archive.arrays["control_target_values"] = dynsyslab::core::NpyArray::fromDoubles({3, 2},
        {1.0, 0.0, 1.15, 0.15, 1.2, 0.0});
    archive.arrays["true_states"] = dynsyslab::core::NpyArray::fromDoubles({3,2},
        {1.0,0.0,1.05,0.1,1.1,0.2});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"native trajectory","integration_form":"dimensionless"})");
    dynsyslab::core::writeNpzAtomic(path, archive);

    const auto payload = dsl::buildNativePayload(path, "control");
    require(payload.figures.size() == 4, "control payload figure count");
    require(payload.figures.front().title == "Control history" &&
        payload.figures.front().traces.front().role == "control",
        "control history must precede state plots");
    require(payload.figures[1].kind == "time" && payload.figures[2].kind == "time",
        "state time series must follow control history");
    require(!payload.figures.back().square && payload.figures.back().kind == "phase", "phase portrait must fill its panel");
    require(payload.figures.back().traces.size() == 3, "target and truth phase overlays");
    require(payload.figures.back().traces[1].lineStyle == "dash", "target line style");
    require(payload.figures.back().traces[2].role=="truth","truth trace role");
    const auto combined = dsl::buildNativePayload(path, "time_frequency");
    require(combined.figures.size() == 5
            && combined.figures[0].title == "Control history"
            && combined.figures[1].kind == "time"
            && combined.figures[2].kind == "time"
            && combined.figures[3].kind == "fourier"
            && combined.figures[4].kind == "fourier",
            "combined view must order control, time series, then Fourier spectra");
    const auto controlSpectrum = dsl::buildNativePayload(path, "control_fourier");
    require(controlSpectrum.figures.size() == 1
            && controlSpectrum.figures[0].kind == "control_fourier"
            && controlSpectrum.figures[0].xLabel == "frequency (cycles / tau)"
            && controlSpectrum.figures[0].yLabel ==
                "control amplitude",
            "dimensionless control Fourier payload");
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"legacy controlled trajectory","control_coordinate_system":"dimensionless"})");
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto taggedControlSpectrum = dsl::buildNativePayload(path, "control_fourier");
    require(taggedControlSpectrum.figures[0].xLabel ==
                "frequency (cycles / tau)" &&
            taggedControlSpectrum.figures[0].yLabel ==
                "control amplitude",
            "legacy control-coordinate tag must determine Fourier units");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}
void legacyDimensionlessArchiveCoordinates() {
    const auto path = std::filesystem::temp_directory_path() /
        "dynsyslab-native-rpe-display-coordinate-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
        {2}, {0.0, 3.0e-6});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {2, 2}, {2.0e-5, 0.0, 1.5e-5, 2.0});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"legacy dimensionless RPE","system":"rayleigh-plesset","integration_form":"dimensionless","time_series_format":"dynsyslab.TimeSeries/v1","system_parameters":{"equilibrium_radius_m":1e-5},"nondimensionalization":{"time_scale_s":3e-6},"execution_backend":{"name":"python-reference"}})");
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto payload = dsl::buildNativePayload(path, "data_overview");
    require(payload.figures.size() == 3, "RPE display payload count");
    require(payload.figures[0].kind == "time" && payload.figures[1].kind == "time" &&
            payload.figures[2].kind == "phase",
            "time series must precede phase space");
    require(std::abs(payload.figures[2].traces[0].x[0] - 2.0) < 1.0e-12,
            "archived physical radius was not converted to r/R0");
    require(std::abs(payload.figures[2].traces[0].y[1] - 0.6) < 1.0e-12,
            "archived physical velocity was not converted to dr/dtau");
    require(std::abs(payload.figures[0].traces[0].x[1] - 1.0) < 1.0e-12 &&
            payload.figures[0].xLabel == "\\tau",
            "archived physical time was not converted to nondimensional time");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}
void observableViews() {
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-native-observable-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({3},{0.0,0.5,1.0});
    archive.arrays["values"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{-1.0,0.0,-0.8,0.2,-0.6,0.1});
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"Duffing observables","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{}})");
    dynsyslab::core::writeNpzAtomic(path,archive);
    const auto time=dsl::buildNativePayload(path,"observable_time");
    require(time.figures.size()>=3 && time.figures[0].dashboardPage==0,
            "complete observable time dashboard");
    const auto maps=dsl::buildNativePayload(path,"observable_heatmap");
    require(maps.figures.size()==time.figures.size() && !maps.figures[0].square
        && maps.figures[0].traces[0].zShape==std::vector<int>({101,101}),
        "all observable phase heat maps fill their panels");
    std::error_code ignored;std::filesystem::remove(path,ignored);
}
void transformedPhaseCoordinates() {
    using namespace dynsyslab::core;
    NpzRecord record;
    record.arrays["values"] = NpyArray::fromDoubles({2,2}, {-2,.3,4,-.2});
    record.arrays["time"] = NpyArray::fromDoubles({2}, {1.275,2.});
    const auto metadata = dsl::json::parse(R"({"system":"rpe","value_space":"history-coordinates","coordinate_names":["echo","trailing echo integral"],"archive_coordinate_system":"dimensionless","integration_form":"dimensionless","system_parameters":{"equilibrium_radius_m":1e-5},"nondimensionalization":{"time_scale_s":3e-6}})");
    const auto payload = dsl::buildNativePayload("echo.npz", "data_overview", record, metadata);
    const auto& phase = payload.figures.back();
    require(phase.xLabel=="echo" && phase.yLabel=="trailing echo integral" &&
            payload.figures.front().yLabel==phase.xLabel, "transformed axis names");
    require(phase.traces.front().x==std::vector<double>({-2,4}) &&
            phase.traces.front().y==std::vector<double>({.3,-.2}) && !phase.alternateCoordinates,
            "physical radius/velocity conversion offered on echo data");
    require(std::none_of(phase.references.begin(),phase.references.end(),
            [](const auto& point){return point.equilibrium;}), "physical equilibrium in transformed plane");
}
void observableControlTarget() {
    using namespace dynsyslab::core;
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-echo-target-payload.npz";
    NpzRecord archive;
    archive.arrays["time"]=NpyArray::fromDoubles({2},{0,.01});
    archive.arrays["values"]=NpyArray::fromDoubles({2,2},{1,0,1.01,.02});
    archive.arrays["control_target_values"]=NpyArray::fromDoubles({2,2},{1,0,1.02,.04});
    archive.arrays["control_observable_target_values"]=NpyArray::fromDoubles({2,1},{.3,.4});
    archive.arrays["control_observable_target_names"]=NpyArray::fromString("far_field_acoustic_echo");
    archive.arrays["metadata_json"]=NpyArray::fromString(
        R"({"system":"rpe","archive_coordinate_system":"dimensionless","control_target_space":"observable"})");
    writeNpzAtomic(path,archive);
    const auto payload=dsl::buildNativePayload(path,"observable_time");
    std::size_t targets=0;
    for(const auto& figure:payload.figures)for(const auto& trace:figure.traces)
        if(trace.role=="control_target") {
            ++targets;require(trace.y==std::vector<double>({.3,.4}),"Viewer must use the supplied echo target");
            require(figure.observableIdentity.find(":far_field_acoustic_echo:0:0")!=std::string::npos,
                "Viewer attached echo target to another observable");
        }
    require(targets==1,"only the echo plot has an observable objective");
    std::filesystem::remove(path);
}
void spectrumView() {
    const auto path = std::filesystem::temp_directory_path() / "dynsyslab-native-spectrum-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({2},
        {{-0.1, 1.2}, {-0.1, -1.2}});
    archive.arrays["continuous_eigenvalues"] = archive.arrays["eigenvalues"];
    archive.arrays["discrete_eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({2},
        {{0.9, 0.1}, {0.9, -0.1}});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"native KROM","eigenvalue_domain":"continuous-time"})");
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto payload = dsl::buildNativePayload(path, "koopman_spectrum");
    require(payload.figures.size() == 1, "spectrum figure count");
    require(payload.figures[0].traces[0].x.size() == 2, "spectrum point count");
    require(payload.figures[0].xLabel == "Re(mu)" &&
        std::abs(payload.figures[0].traces[0].x[0] - 0.9) < 1.0e-12,
        "retained discrete spectrum was mislabeled or replaced by generator eigenvalues");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}
void modeWeightView() {
    using namespace dynsyslab::core;
    const auto dir = std::filesystem::temp_directory_path();
    const auto parentPath = dir / "dynsyslab-mode-weight-parent.npz";
    const auto childPath = dir / "dynsyslab-mode-weight-child.npz";
    // Parent: three modes, two observables; phi_2 has rms 2, so its weight doubles.
    NpzRecord parent;
    parent.arrays["discrete_eigenvalues"] = NpyArray::fromComplex({3}, {{0.9, 0.1}, {0.9, -0.1}, {0.5, 0.0}});
    parent.arrays["continuous_eigenvalues"] = NpyArray::fromComplex({3}, {{-0.1, 1.0}, {-0.1, -1.0}, {-0.7, 0.0}});
    parent.arrays["eigenfunction_values"] = NpyArray::fromComplex({2, 3}, {{1, 0}, {1, 0}, {2, 0}, {1, 0}, {1, 0}, {2, 0}});
    parent.arrays["reconstruction__radius"] = NpyArray::fromComplex({1, 3}, {{0.5, 0}, {0.5, 0}, {0.1, 0}});
    parent.arrays["reconstruction__echo"] = NpyArray::fromComplex({1, 3}, {{0, 0.3}, {0, -0.3}, {0.4, 0}});
    parent.arrays["reconstruction__state"] = NpyArray::fromComplex({2, 3}, {{1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}});
    parent.arrays["reconstruction_names"] = NpyArray::fromString("echo\nradius\nstate");
    parent.arrays["metadata_json"] = NpyArray::fromString(R"({"label":"parent","observables":["radius","echo"]})");
    writeNpzAtomic(parentPath, parent);
    auto index = readNpzIndex(parentPath);
    require(dsl::hasModeWeights(index), "parent has mode weights");
    const auto alone = dsl::buildNativePayload(parentPath, "mode_weights");
    // Three toggleable pages (full unit circle, zoomed mu, lambda), each total + one
    // panel per scalar observable (the multi-row state map is skipped).
    require(alone.figures.size() == 9 && alone.pageLabels.size() == 3, "three planes of three panels");
    require(alone.figures[1].title.find("radius") != std::string::npos
        && alone.figures[2].title.find("echo") != std::string::npos, "observables in training order");
    const auto& stems = alone.figures[2].traces.back();   // echo, after the unit circle
    require(alone.figures[2].traces.size() == 2 && alone.figures[2].traces[0].label == "unit circle"
        && std::abs(stems.z[7] - 0.8) < 1e-12 && std::abs(stems.x[0] - 0.9) < 1e-12,
        "weight is |c| times rms(phi) on a linear height axis, based at mu");
    require(alone.figures[0].zLabel == "|c_j| rms(phi_j)" && alone.figures[0].xLabel == "Re(mu)"
        && alone.figures[0].dashboardPage == 0, "first page: discrete eigenvalues");
    require(alone.figures[3].dashboardPage == 1 && alone.figures[3].xLabel == "Re(mu)", "second page: zoomed mu");
    require(alone.figures[6].dashboardPage == 2 && alone.figures[6].xLabel == "Re(lambda)"
        && alone.figures[8].traces[0].label == "Re = 0"
        && std::abs(alone.figures[8].traces.back().z[7] - 0.8) < 1e-12, "third page: same weights on lambda");
    // Child extracted from the parent, keeping the conjugate pair.
    NpzRecord child;
    child.arrays["continuous_eigenvalues"] = NpyArray::fromComplex({2}, {{-0.1, 1.0}, {-0.1, -1.0}});
    child.arrays["reconstruction__radius"] = NpyArray::fromComplex({1, 2}, {{0.5, 0}, {0.5, 0}});
    child.arrays["reconstruction__echo"] = NpyArray::fromComplex({1, 2}, {{0, 0.3}, {0, -0.3}});
    child.arrays["metadata_json"] = NpyArray::fromString(
        R"({"label":"child","method":"krom-extract","selected_modes":[0,1],"source_krom":")" + parentPath.generic_u8string() + R"("})");
    writeNpzAtomic(childPath, child);
    const auto payload = dsl::buildNativePayload(childPath, "mode_weights");
    require(payload.figures.size() == 9 && payload.warnings.empty(), "child panels from parent, all three planes");
    const auto& radius = payload.figures[1];
    require(radius.title.find("radius") != std::string::npos, "child panels follow the parent's training order");
    require(radius.traces.size() == 3 && radius.traces[1].label == "parent modes"
        && radius.traces[2].label == "kept modes" && radius.traces[2].color == "#ef4444", "parent blue, kept red");
    require(radius.traces[1].x.size() == 9 && radius.traces[2].x.size() == 6, "three parent stems, two kept");
    require(radius.title.find("2 of 3 modes kept") != std::string::npos, "kept count in title");
    std::error_code ignored; std::filesystem::remove(parentPath, ignored); std::filesystem::remove(childPath, ignored);
}
void kromDashboardView() {
    const auto path = std::filesystem::temp_directory_path() /
        "dynsyslab-native-krom-dashboard-test.npz";
    dynsyslab::core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues = {{-0.1, 1.0}, {-0.1, -1.0}};
    fit.coefficients = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.2},
                        {0.0, 0.0}, {1.0, 0.0}, {0.0, -0.2}};
    fit.exponents = {0, 0, 1, 0, 0, 1};
    fit.center = {0.0, 0.0}; fit.scale = {1.0, 1.0};
    fit.stateDimension = 2; fit.dynamicalStateDimension = 2;
    fit.coordinateObservableNames = {"position", "total_energy"};
    fit.featureCount = 3;
    dynsyslab::core::SystemConfiguration system;
    system.identifier = "duffing"; system.preset = "canonical";
    dynsyslab::core::KoopmanModel model(std::move(fit), system);
    const double dt = 0.1;
    std::map<std::string, dynsyslab::core::NpyArray> extra;
    extra["training_time"] = dynsyslab::core::NpyArray::fromDoubles(
        {3}, {0.0, dt, 2.0 * dt});
    extra["training_states"] = dynsyslab::core::NpyArray::fromDoubles(
        {3, 2}, {-1.0, 0.0, -0.8, 0.3, -0.6, 0.1});
    std::vector<std::complex<double>> discrete;
    for (const auto value : model.eigenfunctions().continuousEigenvalues)
        discrete.push_back(std::exp(value * dt));
    extra["discrete_eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({2}, discrete);
    extra["eigenvalues"] = extra["discrete_eigenvalues"];
    dynsyslab::core::savePolynomialKoopmanModel(path, model,
        R"({"label":"two-mode Duffing KROM","delta_t":0.1})", std::move(extra));

    const auto payload = dsl::buildNativePayload(path, "krom_dashboard");
    require(payload.layout == "krom" && payload.figures.size() == 5, "default selected KROM workspace");
    require(payload.figures[0].dashboardPage == -2 && payload.figures[0].traces.size() == 3,
        "persistent retained spectrum");
    require(payload.figures[0].traces[0].lineStyle == "dash", "dashed unit circle");
    require(payload.figures[1].kind == "krom_heatmap" && payload.figures[2].kind == "krom_time" &&
        payload.figures[3].kind == "krom_heatmap" && payload.figures[4].kind == "krom_time",
        "each phase-space component is paired with its matching training history");
    require(std::abs(payload.figures[2].traces[0].y[0] + 1) < 1e-12 &&
        std::abs(payload.figures[4].traces[0].y[0] + .05) < 1e-12,
        "eigenfunctions use the native position/energy coordinate evaluator");
    const auto& map = payload.figures[3].traces[0];
    const double point[] = {map.x[40],map.y[60]};
    require(std::abs(map.z[60*101+40]-model.evaluate(point)[0].imag()) < 1e-12,
        "eigenfunction heat map evaluates native state coordinates");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}
void synchronizedComparisonView() {
    dsl::Payload first{"first.npz", "time", "First", 1};
    first.pageLabels = {"State", "Error"};
    dsl::Figure history{"time", "Radius", "time", "radius", false};
    history.dashboardPage = 0;
    history.traces.push_back({"actual", "", "#3bbfb8", "solid", "", false,
                              false, {0.0, 1.0}, {1.0, 2.0}});
    dsl::Figure error{"control_error", "Error", "x", "y", true};
    error.dashboardPage = 1;
    dsl::Trace field; field.label = "relative"; field.heatmap = true;
    field.x = {-1.0, 1.0}; field.y = {-1.0, 1.0}; field.z = {1.0, 2.0, 3.0, 4.0};
    field.zShape = {2, 2}; error.traces.push_back(field);
    first.figures = {history, error};
    dsl::Payload second = first; second.source = "second.npz"; second.tileName = "Second";
    second.figures[0].traces[0].y = {1.1, 1.9};
    const auto compared = dsl::comparePayload(first, "First", second, "Second");
    require(compared.pageLabels == std::vector<std::string>({"State", "Error"}),
            "comparison page synchronization");
    require(compared.figures.size() == 3 && compared.figures[0].traces.size() == 2,
            "compatible curves were not overlaid");
    require(compared.figures[0].traces[1].lineStyle == "dash"
        && compared.figures[0].traces[1].label.find("Second") == 0,
        "comparison trace identity is not visible");
    require(compared.figures[1].dashboardPage == 1
        && compared.figures[2].dashboardPage == 1
        && compared.figures[1].traces[0].heatmap,
        "non-overlayable heat maps were not paired side by side");
}

// Closed-loop estimator archive: output grid (time) plus a finer controller clock.
dynsyslab::core::NpzRecord estimatorArchive() {
    using dynsyslab::core::NpyArray;
    dynsyslab::core::NpzRecord archive;
    std::vector<double> t, x, target, te, applied, estimate, measured, predicted, innovation;
    for (int i = 0; i <= 10; ++i) {
        t.push_back(i); x.push_back(1.0 + 0.01 * i); x.push_back(0.01);
        target.push_back(1.0); target.push_back(0.0);
    }
    for (int i = 0; i <= 100; ++i) {
        te.push_back(0.1 * i); applied.push_back(i % 2 ? 0.5 : -0.5);
        estimate.push_back(1.0 + 0.001 * i); estimate.push_back(0.0);
        measured.push_back(std::sin(0.1 * i)); predicted.push_back(std::sin(0.1 * i) + 0.01);
        innovation.push_back(-0.01);
    }
    archive.arrays["time"] = NpyArray::fromDoubles({11}, t);
    archive.arrays["values"] = NpyArray::fromDoubles({11, 2}, x);
    archive.arrays["control_target_values"] = NpyArray::fromDoubles({11, 2}, target);
    archive.arrays["control_values"] = NpyArray::fromDoubles({11, 1}, std::vector<double>(11, 0.0));
    archive.arrays["estimator_time"] = NpyArray::fromDoubles({101}, te);
    archive.arrays["applied_input"] = NpyArray::fromDoubles({101, 1}, applied);
    archive.arrays["state_estimate"] = NpyArray::fromDoubles({101, 2}, estimate);
    archive.arrays["measurement_values"] = NpyArray::fromDoubles({101, 1}, measured);
    archive.arrays["measurement_estimate"] = NpyArray::fromDoubles({101, 1}, predicted);
    archive.arrays["innovation"] = NpyArray::fromDoubles({101, 1}, innovation);
    archive.arrays["estimator_gain"] = NpyArray::fromDoubles({4, 1}, {1, 2, 3, 4});
    archive.arrays["metadata_json"] = NpyArray::fromString(
        R"({"label":"echo-only loop","system":"duffing","controller_kind":"klqr"})");
    return archive;
}

void recordWindowView() {
    const auto archive = estimatorArchive();
    const auto windowed = dsl::windowRecord(archive, {2.0, 5.0});
    require(windowed.at("time").shape[0] == 4 && windowed.at("values").shape == std::vector<std::size_t>({4, 2})
        && windowed.at("values").doubles().front() == 1.02, "output-grid arrays follow the window");
    require(windowed.at("estimator_time").shape[0] == 31 && windowed.at("applied_input").shape[0] == 31
        && std::abs(windowed.at("estimator_time").doubles().front() - 2.0) < 1e-12,
        "controller-clock arrays follow their own clock");
    require(windowed.at("estimator_gain").shape[0] == 4, "unaligned arrays are kept whole");
    bool rejected = false;
    try { (void)dsl::windowRecord(archive, {20.0, 30.0}); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "an empty window is rejected");
}

void controllerClockAndEstimatorView() {
    const auto path = std::filesystem::temp_directory_path() / "dynsyslab-native-estimator-view-test.npz";
    dynsyslab::core::writeNpzAtomic(path, estimatorArchive());
    const auto time = dsl::buildNativePayload(path, "time");
    const auto& control = time.figures.front();
    require(control.title == "Control history" && control.traces.front().x.size() == 101
        && control.traces.front().y[1] == 0.5, "control history uses the applied input on the controller clock");
    const auto& radius = time.figures[1];
    require(std::any_of(radius.traces.begin(), radius.traces.end(), [](const dsl::Trace& trace) {
        return trace.role == "estimate" && trace.x.size() == 101; }), "state panels overlay the estimate");
    const auto phase = dsl::buildNativePayload(path, "phase");
    require(phase.figures.front().traces.back().role == "estimate", "phase space overlays the estimate");
    const auto estimator = dsl::buildNativePayload(path, "estimator");
    require(estimator.figures.size() == 2 && estimator.figures[0].traces.size() == 2
        && estimator.figures[0].traces[1].role == "estimate" && estimator.figures[1].title == "Innovation",
        "estimator measurement and innovation panels");
    const auto spectrum = dsl::buildNativePayload(path, "control_fourier");
    // 101 samples at 0.1: last one-sided bin 50 / 10.1 = 4.95, far above the output grid's 0.5.
    require(std::abs(spectrum.figures.front().traces.front().x.back() - 50.0 / 10.1) < 1e-9,
        "control spectrum reaches the controller-clock Nyquist frequency");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}

void coloredComparisonView() {
    const auto run = [](double offset, std::vector<double> target) {
        dsl::Payload payload{"run.npz", "time", "run", 1};
        dsl::Figure figure{"time", "r versus time", "t", "r", false};
        figure.traces.push_back({"r", "data", "", "solid", "", false, false, {0, 1, 2}, {offset, offset, offset}});
        figure.traces.push_back({"target", "control_target", "#ef4444", "dash", "", false, false, {0, 1, 2}, target});
        payload.figures = {figure};
        return payload;
    };
    const auto ramped = std::vector<double>{1.0, 1.1, 1.2}, step = std::vector<double>{1.2, 1.2, 1.2};
    const auto compared = dsl::comparePayloads({{"A", run(1.0, ramped), "#004DB8"},
        {"B", run(1.1, ramped), "#D8242A"}, {"C", run(1.2, step), "#DE9B00"}});
    const auto& traces = compared.figures.front().traces;
    require(traces.size() == 5, "identical targets are drawn once per distinct curve");
    require(traces[0].label == "A, B — target" && traces[0].color == "#202020" && traces[0].lineStyle == "solid",
        "the widest shared target is solid dark grey and names its runs");
    require(traces[1].label == "C — target" && traces[1].lineStyle == "dash", "the other target is dashed");
    require(traces[2].color == "#004DB8" && traces[3].color == "#D8242A" && traces[4].color == "#DE9B00"
        && traces[4].lineStyle == "solid", "run colors replace the dash convention");
    const auto same = dsl::comparePayloads({{"A", run(1.0, ramped), "#004DB8"}, {"B", run(1.1, ramped), "#D8242A"}});
    require(same.figures.front().traces.front().label == "target", "a target shared by every run is unprefixed");
    auto open = run(1.2, step); open.figures.front().traces.pop_back();  // an open-loop reference has no target
    const auto mixed = dsl::comparePayloads({{"A", run(1.0, ramped), "#004DB8"}, {"B", run(1.1, ramped), "#D8242A"},
        {"exact", open, "#202020"}});
    require(mixed.figures.front().traces.front().label == "target" && mixed.figures.front().traces.size() == 4,
        "runs without a target do not prefix the shared target");
}

void controlErrorView() {
    const auto path = std::filesystem::temp_directory_path() / "dynsyslab-native-control-error-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["control_error_axis0"] = dynsyslab::core::NpyArray::fromDoubles({3}, {-1,0,1});
    archive.arrays["control_error_axis1"] = dynsyslab::core::NpyArray::fromDoubles({2}, {-2,2});
    archive.arrays["constant_relative_error"] = dynsyslab::core::NpyArray::fromDoubles(
        {2,3}, {0,0.1,0.2,0.3,0.4,0.5});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"native control error"})");
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto payload = dsl::buildNativePayload(path, "control_map_error");
    require(payload.figures.size() == 1, "control error figure count");
    require(payload.figures[0].square && payload.figures[0].traces[0].heatmap,
            "control error heatmap geometry");
    require(payload.figures[0].traces[0].zShape == std::vector<int>({2,3}),
            "control error heatmap shape");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}

void observabilityView() {
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-native-observability-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["observability_singular_values"]=
        dynsyslab::core::NpyArray::fromDoubles({2},{10.0,0.1});
    archive.arrays["observable_directions"]=
        dynsyslab::core::NpyArray::fromDoubles({2,2},{1.0,0.0,0.0,1.0});
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"native observability"})");
    dynsyslab::core::writeNpzAtomic(path,archive);
    const auto payload=dsl::buildNativePayload(path,"observability");
    require(payload.figures.size()==2,"observability figure count");
    require(payload.figures[0].traces[0].y[0]==1.0,"observability log spectrum");
    require(payload.figures[1].traces[0].heatmap,"observability direction heatmap");
    std::error_code ignored; std::filesystem::remove(path,ignored);
}

void klqeTrajectoryView() {
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-native-klqe-view-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({2},{0.0,1.0});
    archive.arrays["values"]=dynsyslab::core::NpyArray::fromDoubles({2,2},{0.0,0.0,1.0,0.2});
    archive.arrays["observability_singular_values"]=dynsyslab::core::NpyArray::fromDoubles({2},{2.0,1.0});
    archive.arrays["observable_directions"]=dynsyslab::core::NpyArray::fromDoubles({2,2},{1,0,0,1});
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"KLQE estimate","method":"klqe","system":"duffing"})");
    dynsyslab::core::writeNpzAtomic(path,archive);
    require(dsl::buildNativePayload(path,"time").figures[0].kind=="time",
            "KLQE time view was masked by observability dispatch");
    require(dsl::buildNativePayload(path,"observability").figures.size()==2,
            "KLQE observability view");
    std::error_code ignored;std::filesystem::remove(path,ignored);
}

void latticeView() {
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-native-lattice-view-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["lattice_initial"]=dynsyslab::core::NpyArray::fromInt64s({2,2},{1,-1,-1,1});
    archive.arrays["lattice_final"]=dynsyslab::core::NpyArray::fromInt64s({2,2},{1,1,1,1});
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"Ising lattice","system":"ising","system_family":"lattice"})");
    dynsyslab::core::writeNpzAtomic(path,archive);
    const auto payload=dsl::buildNativePayload(path,"lattice");
    require(payload.figures.size()==2 && payload.figures[0].square
        && payload.figures[1].traces[0].z.back()==1.0,"native lattice payload");
    std::error_code ignored;std::filesystem::remove(path,ignored);
}

void decompositionViews() {
    const auto root=std::filesystem::temp_directory_path();
    const auto pcaPath=root/"dynsyslab-native-pca-payload-test.npz";
    dynsyslab::core::NpzRecord pca;
    pca.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({3},{0.0,.1,.2});
    pca.arrays["singular_values"]=dynsyslab::core::NpyArray::fromDoubles({2},{4.0,1.0});
    pca.arrays["explained_variance_ratio"]=
        dynsyslab::core::NpyArray::fromDoubles({2},{.8,.2});
    pca.arrays["components"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{1.0,0.0,0.0,1.0});
    pca.arrays["scores"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{0.0,1.0,1.0,0.0,2.0,-1.0});
    pca.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"native PCA"})");
    dynsyslab::core::writeNpzAtomic(pcaPath,pca);
    const auto pcaPayload=dsl::buildNativePayload(pcaPath,"data_overview");
    require(pcaPayload.figures.size()==4,"PCA payload figure count");
    require(pcaPayload.figures[2].traces[0].heatmap,
            "PCA component heatmap");
    const auto robustPath=root/"dynsyslab-native-robust-pca-payload-test.npz";
    dynsyslab::core::NpzRecord robust;
    robust.arrays["time"]=pca.arrays["time"];
    robust.arrays["low_rank"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{1.0,2.0,2.0,4.0,3.0,6.0});
    robust.arrays["sparse"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{0.0,0.0,0.0,5.0,0.0,0.0});
    robust.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"native robust PCA"})");
    dynsyslab::core::writeNpzAtomic(robustPath,robust);
    const auto robustPayload=dsl::buildNativePayload(robustPath,"data_overview");
    require(robustPayload.figures.size()==4,"robust PCA payload figure count");
    require(robustPayload.figures[1].traces[0].role=="sparse",
            "robust PCA sparse trace role");
    std::error_code ignored; std::filesystem::remove(pcaPath,ignored);
    std::filesystem::remove(robustPath,ignored);
}

void freeCampaignCatalog() {
    const auto root = std::filesystem::temp_directory_path()
        / "dynsyslab-native-free-campaign-catalog";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    const auto archiveRoot = root / "archives";
    const auto solution = archiveRoot / "solutions" / "ising" / "sample.npz";
    std::filesystem::create_directories(solution.parent_path());
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({1}, {0.0});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles({1, 2}, {0.0, -2.0});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString("{}");
    dynsyslab::core::writeNpzAtomic(solution, archive);
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "ising";
    configuration.parameters = {{"temperature", 2.0}};
    dynsyslab::core::BatchRecord batch;
    batch.id = "batch"; batch.label = "temperature sweep";
    batch.createdUtc = "now"; batch.dataRole = "campaign";
    batch.campaignId = "campaign"; batch.campaignKind = "temperature_sweep";
    batch.tags = {"free", "campaign"};
    dynsyslab::core::SystemRunRecord system;
    system.id = "system"; system.systemIdentifier = "ising";
    system.systemLabel = "2-D Ising model"; system.parameterSetKey =
        dynsyslab::core::parameterSetKey(configuration);
    system.coordinateSystem = "discrete"; system.stateNames = {"m", "energy"};
    dynsyslab::core::ClusterRecord cluster; cluster.id = "cluster";
    dynsyslab::core::TrajectoryRecord trajectory;
    trajectory.id = "trajectory"; trajectory.label = "sample";
    trajectory.solutionArchive = solution; trajectory.tags = {"free", "campaign"};
    cluster.trajectories.push_back(trajectory);
    system.clusters.push_back(cluster); batch.systems.push_back(system);
    dynsyslab::core::writeHierarchyAtomic(
        archiveRoot / "hierarchies" / "ising" / "campaign.hierarchy.json", batch);
    const auto catalog = dsl::Catalog::scan(root, archiveRoot);
    const int index = catalog.indexOf("trajectory");
    require(index >= 0, "free campaign trajectory catalog entry");
    const auto& item = catalog.items[static_cast<std::size_t>(index)];
    require(std::find(item.tags.begin(), item.tags.end(), "free") != item.tags.end(),
            "free campaign retains free modifier");
    require(std::find(item.tags.begin(), item.tags.end(), "forced") == item.tags.end(),
            "free campaign is not mislabeled forced");
    require(std::find(item.views.begin(), item.views.end(), "control") == item.views.end(),
            "free campaign does not advertise control history");
    std::filesystem::remove_all(root, ignored);
}

void metadataDrivenCatalogRelations() {
    const auto root=std::filesystem::temp_directory_path()/
        "dynsyslab-native-metadata-catalog";
    std::error_code ignored; std::filesystem::remove_all(root,ignored);
    const auto archiveRoot=root/"archives";
    const auto source=archiveRoot/"solutions"/"duffing"/"source.npz";
    const auto derived=archiveRoot/"derived"/"koopman-krom-extract"/"opaque-name.npz";
    const auto legacy=archiveRoot/"derived"/"koopman-krom-extract"/"legacy-name.npz";
    const auto spectrum=archiveRoot/"derived"/"koopman-hankel-dmd"/"spectrum.npz";
    const auto expansion=archiveRoot/"solutions"/"duffing"/"domain-expansion.npz";
    std::filesystem::create_directories(source.parent_path());
    std::filesystem::create_directories(derived.parent_path());
    std::filesystem::create_directories(spectrum.parent_path());
    dynsyslab::core::NpzRecord trajectory;
    trajectory.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({2},{0.0,1.0});
    trajectory.arrays["values"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{0.0,0.0,1.0,1.0});
    trajectory.arrays["control_values"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,1},{0.0,0.25});
    trajectory.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"format":"dynsyslab.solution/v1","label":"Duffing source","system":"duffing","system_label":"Duffing oscillator"})");
    dynsyslab::core::writeNpzAtomic(source,trajectory);
    dynsyslab::core::NpzRecord krom;
    krom.arrays["continuous_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {1},{{-0.1,1.0}});
    krom.arrays["discrete_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {1},{{0.9,0.1}});
    krom.arrays["coefficients"]=dynsyslab::core::NpyArray::fromComplex(
        {1,1},{{1.0,0.0}});
    krom.arrays["polynomial_exponents"]=dynsyslab::core::NpyArray::fromDoubles(
        {1,2},{1.0,0.0});
    krom.arrays["state_center"]=dynsyslab::core::NpyArray::fromDoubles({2},{0.0,0.0});
    krom.arrays["state_scale"]=dynsyslab::core::NpyArray::fromDoubles({2},{1.0,1.0});
    krom.arrays["training_states"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{0.0,0.0,1.0,1.0});
    krom.arrays["control_error_axis0"]=dynsyslab::core::NpyArray::fromDoubles(
        {2},{-1.0,1.0});
    krom.arrays["control_error_axis1"]=dynsyslab::core::NpyArray::fromDoubles(
        {2},{-1.0,1.0});
    krom.arrays["constant_relative_error"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{0.0,0.1,0.2,0.3});
    const std::string metadata="{\"format\":\"dynsyslab.derived/v1\","
        "\"label\":\"Independent minimal KROM\",\"method\":\"krom-extract\","
        "\"system\":\"duffing\",\"system_label\":\"Duffing oscillator\","
        "\"source_archive\":\""+source.generic_u8string()+"\","
        "\"selected_modes\":[7],\"control_maps\":[\"exact\",\"bilinear\"],"
        "\"independent_archive\":true}";
    krom.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(metadata);
    dynsyslab::core::writeNpzAtomic(derived,krom);
    dynsyslab::core::NpzRecord expanded=trajectory;
    expanded.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        "{\"format\":\"dynsyslab.solution/v1\",\"label\":\"Domain expansion\","
        "\"system\":\"duffing\",\"training_preflight\":{\"enabled\":true,"
        "\"source_krom\":\""+derived.generic_u8string()+"\"}}");
    dynsyslab::core::writeNpzAtomic(expansion,expanded);
    dynsyslab::core::NpzRecord legacyKrom;
    legacyKrom.arrays["eigenvalues"]=krom.arrays["discrete_eigenvalues"];
    legacyKrom.arrays["coefficients"]=krom.arrays["coefficients"];
    legacyKrom.arrays["training_states"]=krom.arrays["training_states"];
    legacyKrom.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"analysis":"koopman-krom-extract","label":"Legacy retained KROM","analysis_metadata":{"krom":{"delta_t":0.1}}})");
    dynsyslab::core::writeNpzAtomic(legacy,legacyKrom);
    dynsyslab::core::NpzRecord hdmd;
    hdmd.arrays["continuous_eigenvalues"]=krom.arrays["continuous_eigenvalues"];
    hdmd.arrays["eigenfunction_values"]=dynsyslab::core::NpyArray::fromComplex(
        {2,1},{{1.0,0.0},{0.9,0.1}});
    hdmd.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"format":"dynsyslab.derived/v1","label":"H-DMD spectrum","method":"hankel-dmd","system":"duffing"})");
    dynsyslab::core::writeNpzAtomic(spectrum,hdmd);
    const auto catalog=dsl::Catalog::scan(root,archiveRoot);
    const dsl::CatalogItem* sourceItem=nullptr; const dsl::CatalogItem* kromItem=nullptr;
    const dsl::CatalogItem* spectrumItem=nullptr;
    const dsl::CatalogItem* legacyItem=nullptr;
    const dsl::CatalogItem* expansionItem=nullptr;
    for (const auto& item : catalog.items) {
        if (item.label=="Duffing source") sourceItem=&item;
        if (item.label=="Independent minimal KROM") kromItem=&item;
        if (item.label=="H-DMD spectrum") spectrumItem=&item;
        if (item.label=="Legacy retained KROM") legacyItem=&item;
        if (item.label=="Domain expansion") expansionItem=&item;
    }
    require(sourceItem && kromItem,"metadata-driven catalog objects are missing");
    require(std::find(sourceItem->views.begin(),sourceItem->views.end(),
        "observable_time")!=sourceItem->views.end()
        && std::find(sourceItem->views.begin(),sourceItem->views.end(),
        "observable_heatmap")!=sourceItem->views.end(),
        "two-state trajectory did not advertise native observable views");
    require(std::find(sourceItem->views.begin(),sourceItem->views.end(),
        "control_fourier")!=sourceItem->views.end(),
        "controlled trajectory did not advertise control spectra");
    require(spectrumItem && spectrumItem->kind=="koopman_spectrum"
        && spectrumItem->category=="Koopman spectra",
        "H-DMD spectrum was incorrectly exposed as a control-ready KROM");
    require(legacyItem && legacyItem->kind=="krom" && !legacyItem->views.empty()
        && legacyItem->views.front()=="krom_dashboard",
        "legacy retained KROM lost its native dashboard route");
    require(kromItem->kind=="krom" && kromItem->category=="KROMs",
        "metadata-driven KROM classification");
    require(!kromItem->views.empty() && kromItem->views.front()=="krom_dashboard",
        "renderable KROM does not default to the standard dashboard");
    require(std::find(kromItem->views.begin(),kromItem->views.end(),
        "control_map_error")!=kromItem->views.end(),"KROM error-map view discovery");
    require(kromItem->properties.at("source mode indices")=="7"
        && kromItem->properties.at("parent required at runtime")=="no",
        "independent KROM properties");
    require(std::any_of(kromItem->relations.begin(),kromItem->relations.end(),
        [&](const dsl::Relation& relation){return relation.target==sourceItem->id
            && relation.kind=="extracted from";}),
        "archive lineage relation was not resolved");
    require(std::any_of(sourceItem->relations.begin(),sourceItem->relations.end(),
        [&](const dsl::Relation& relation){return relation.target==kromItem->id
            && relation.kind=="has extracted KROM";}),
        "archive reverse lineage relation was not generated");
    require(expansionItem && std::any_of(expansionItem->relations.begin(),
        expansionItem->relations.end(),[&](const dsl::Relation& relation){
            return relation.target==kromItem->id
                && relation.kind=="expands training domain of";}),
        "training-domain source KROM relation was not generated");
    require(std::any_of(kromItem->relations.begin(),kromItem->relations.end(),
        [&](const dsl::Relation& relation){return relation.target==expansionItem->id
            && relation.kind=="has domain-expansion data";}),
        "training-domain reverse relation was not generated");
    const auto update = dynsyslab::core::reconcileArchiveIndex(archiveRoot);
    const auto indexedCatalog = dsl::Catalog::indexed(root, archiveRoot);
    require(update.indexed >= 5 && indexedCatalog.indexOf(kromItem->id) >= 0 &&
            indexedCatalog.indexOf(sourceItem->id) >= 0,
            "central archive index did not reproduce catalog metadata and lineage");
    std::filesystem::remove_all(root,ignored);
}

void pdeFieldView() {
    const auto path = std::filesystem::temp_directory_path()
        / "dynsyslab-native-pde-field-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
        {3}, {0.0, 0.1, 0.2});
    archive.arrays["space"] = dynsyslab::core::NpyArray::fromDoubles(
        {4}, {0.0, 0.25, 0.5, 0.75});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {3, 4}, {0, 1, 0, -1, 0, 0.8, 0, -0.8, 0, 0.6, 0, -0.6});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"heat field","system_family":"pde"})");
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto payload = dsl::buildNativePayload(path, "field");
    require(payload.figures.size() == 2, "PDE field and profile figures");
    require(payload.figures[0].traces[0].heatmap,
            "PDE space-time heatmap");
    require(payload.figures[0].traces[0].zShape == std::vector<int>({3, 4}),
            "PDE heatmap shape");
    require(payload.figures[1].traces[0].y.back() == -0.6,
            "PDE final profile values");
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
}

}  // namespace

int main() {
    try {
        trajectoryViews(); legacyDimensionlessArchiveCoordinates(); transformedPhaseCoordinates(); observableViews();
        spectrumView(); modeWeightView(); kromDashboardView(); synchronizedComparisonView();
        controlErrorView(); observabilityView(); klqeTrajectoryView(); latticeView();
        recordWindowView(); controllerClockAndEstimatorView(); coloredComparisonView();
        decompositionViews();
        freeCampaignCatalog(); metadataDrivenCatalogRelations(); pdeFieldView();
        std::cout << "dynsyslab_native_payload_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_native_payload_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
