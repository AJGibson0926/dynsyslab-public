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
void dimensionalArchiveCoordinates() {
    const auto path = std::filesystem::temp_directory_path() /
        "dynsyslab-native-rpe-display-coordinate-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
        {2}, {0.0, 3.0e-6});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {2, 2}, {2.0e-5, 0.0, 1.5e-5, 2.0});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"dimensional RPE archive","system":"rayleigh-plesset","integration_form":"dimensionless","archive_coordinate_system":"dimensional","system_parameters":{"equilibrium_radius_m":1e-5},"nondimensionalization":{"time_scale_s":3e-6}})");
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
    auto domain = dynsyslab::core::FunctionDomain::fromCatalog({"position", "total_energy"}, 2);
    domain.physicalStateDimension = 2;
    auto fit = dynsyslab::core::polynomialEigenfunctions(2, {0, 0, 1, 0, 0, 1},
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.2}, {0.0, 0.0}, {1.0, 0.0}, {0.0, -0.2}},
        dynsyslab::core::Spectrum::fromContinuousRates({{-0.1, 1.0}, {-0.1, -1.0}}),
        {0.0, 0.0}, {1.0, 1.0}, std::move(domain));
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
    for (const auto value : model.spectrum().continuousRates())
        discrete.push_back(std::exp(value * dt));
    extra["discrete_eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({2}, discrete);
    extra["eigenvalues"] = extra["discrete_eigenvalues"];
    dynsyslab::core::saveKoopmanModel(path, model,
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

}  // namespace

int main() {
    try {
        trajectoryViews(); dimensionalArchiveCoordinates(); transformedPhaseCoordinates(); observableViews();
        spectrumView(); modeWeightView(); kromDashboardView(); synchronizedComparisonView();
        controlErrorView(); observabilityView(); klqeTrajectoryView();
        recordWindowView(); controllerClockAndEstimatorView(); coloredComparisonView();
        std::cout << "dynsyslab_native_payload_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_native_payload_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
