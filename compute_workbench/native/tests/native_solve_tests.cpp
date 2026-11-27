#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Campaign.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>

#include <filesystem>
#include <cmath>
#include <iostream>
#include <stdexcept>

void verifyStateControllers(const std::filesystem::path& root);  // native_solve_controller_tests.cpp

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::string jsonEscape(const std::string& value) {
    std::string result;
    for (char item : value) { if (item == '\\' || item == '"') result.push_back('\\'); result.push_back(item); }
    return result;
}
void solveAndArchive() {
    const auto root = std::filesystem::temp_directory_path() / "dynsyslab-native-solve-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    const auto job = dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"native-test",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":0.1,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,
                            "values":[[0.1,0]],"ranges":[]},
      "controller":{"kind":"waveform","settings":{"amplitude":0.2,"frequency":1},
                    "target":{"kind":"constant","state":[-1,0]}},
      "solver":{"method":"BDF","rtol":1e-8,"atol":1e-10},
      "archive":{"enabled":true}
    })");
    const auto result = dynsyslab::workflow::executeSolve(job, root, root / "progress.json");
    require(result.archives.size() == 1, "one trajectory archive");
    const auto archive = dynsyslab::core::readNpz(result.archives.front());
    require(archive.at("values").shape == std::vector<std::size_t>({11, 2}), "state matrix shape");
    require(archive.at("control_values").shape == std::vector<std::size_t>({11, 1}), "control history shape");
    require(!archive.contains("control_target_values"),
            "prescribed waveform must not create a target history");
    require(archive.metadataJson().find("native-cpp") != std::string::npos, "native provenance");
    const auto solutionMetadata = dynsyslab::core::json::parse(archive.metadataJson());
    require(solutionMetadata.at("control_coordinate_system").string() == "dimensionless",
            "control coordinate metadata");
    require(solutionMetadata.at("control_units").string() == "dimensionless system input",
            "control unit metadata");
    require(solutionMetadata.at("control_scale_to_dimensional").number() == 1.0,
            "control scale metadata");
    require(solutionMetadata.at("controller_settings_coordinates").string() == "model",
            "controller-setting coordinate metadata");
    require(std::filesystem::is_regular_file(root / "progress.json"), "progress record");
    require(std::filesystem::is_regular_file(result.hierarchyManifest), "hierarchy manifest");
    const auto hierarchy = dynsyslab::core::json::load(result.hierarchyManifest.wstring());
    require(hierarchy.at("schema").string() == "dynsyslab.hierarchy/v2", "hierarchy schema");
    require(hierarchy.at("systems").array().size() == 1, "one system run");
    require(hierarchy.at("systems").array().front().at("clusters").array().size() == 1,
            "one trajectory cluster");
    const auto& systemRun = hierarchy.at("systems").array().front();
    require(!systemRun.at("equations_fingerprint").string().empty(),
            "hierarchy equation fingerprint");
    require(systemRun.at("resolved_parameters").at("delta").number() == 0.2,
            "hierarchy resolved parameters");
    const auto& waveformSeries = systemRun.at("clusters").array().front()
        .at("trajectories").array().front().at("time_series").array();
    require(waveformSeries.size() == 2
            && waveformSeries[0].at("role").string() == "state"
            && waveformSeries[1].at("role").string() == "control",
            "prescribed waveform hierarchy must contain state and input only");

    const auto delayJob = dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"native-delay-test",
      "system":{"name":"mackey_glass","family":"dde","preset":"classic_chaos",
                "parameters":{"history_value":1.2}},
      "time":{"start":0,"stop":1,"output_step":0.05,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":2,"clusters":1,"seed":0,
                            "values":[[1.2],[1.0]],"ranges":[]},
      "controller":{"kind":"none","settings":{},
                    "target":{"kind":"constant","state":[0]}},
      "solver":{"method":"method_of_steps","step_size":0.05},
      "archive":{"enabled":true,"data_role":"free_training"}
    })");
    const auto delayRun = dynsyslab::workflow::executeJob(
        delayJob, root, root / "delay-progress.json");
    require(delayRun.archives.size() == 2, "two delay trajectory archives");
    const auto delayArchive = dynsyslab::core::readNpz(delayRun.archives.front());
    require(delayArchive.at("values").shape == std::vector<std::size_t>({21, 1}),
            "delay state matrix shape");
    require(delayArchive.metadataJson().find("\"system_family\":\"dde\"")
            != std::string::npos, "delay family provenance");
    require(std::filesystem::is_regular_file(delayRun.manifest),
            "delay hierarchy manifest");
    const auto delayHierarchy = dynsyslab::core::json::load(delayRun.manifest.wstring());
    require(delayHierarchy.at("systems").array().front().at("resolved_parameters")
            .at("beta").number() == 0.2, "delay resolved hierarchy parameters");

    const auto latticeJob = dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"native-lattice-test",
      "system":{"name":"ising","family":"lattice","preset":"critical",
                "parameters":{"size":4,"temperature":2.269}},
      "time":{"start":0,"stop":3,"output_step":1,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"default","count":2,"clusters":1,"seed":42,
                            "values":[],"ranges":[]},
      "controller":{"kind":"none","settings":{},
                    "target":{"kind":"constant","state":[0,0]}},
      "solver":{"method":"metropolis"},
      "archive":{"enabled":true,"data_role":"free_training"}
    })");
    const auto latticeRun = dynsyslab::workflow::executeJob(
        latticeJob, root, root / "lattice-progress.json");
    require(latticeRun.archives.size() == 2, "two lattice trajectory archives");
    const auto latticeArchive = dynsyslab::core::readNpz(latticeRun.archives.front());
    require(latticeArchive.at("values").shape == std::vector<std::size_t>({3, 2}),
            "lattice observable matrix shape");
    require(latticeArchive.at("lattice_initial").shape
            == std::vector<std::size_t>({4, 4}), "initial lattice snapshot shape");
    require(latticeArchive.at("lattice_final").type == dynsyslab::core::ArrayType::int64,
            "integer final lattice snapshot");
    require(latticeArchive.metadataJson().find("\"system_family\":\"lattice\"")
            != std::string::npos, "lattice family provenance");
    require(std::filesystem::is_regular_file(latticeRun.manifest),
            "lattice hierarchy manifest");
    const auto latticeHierarchy = dynsyslab::core::json::load(latticeRun.manifest.wstring());
    require(latticeHierarchy.at("systems").array().front().at("resolved_parameters")
            .at("size").number() == 4.0, "lattice resolved hierarchy parameters");

    const auto pdeJob = dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"native-pde-test",
      "system":{"name":"heat_equation_1d","family":"pde",
                "preset":"periodic_diffusion",
                "parameters":{"points":16,"diffusivity":0.001},
                "pde_partitioning":{"enabled":true,"partition_shape":[4],
                    "ghost_width":1,"local_workers":2}},
      "time":{"start":0,"stop":0.02,"output_step":0.01,
              "coordinate_system":"automatic"},
      "initial_conditions":{"mode":"default","count":2,"clusters":1,
                            "seed":0,"values":[],"ranges":[]},
      "controller":{"kind":"none","settings":{},
                    "target":{"kind":"constant","state":[]}},
      "solver":{"method":"RK4","step_size":0.005},
      "archive":{"enabled":true,"data_role":"free_training"}
    })");
    const auto pdeRun = dynsyslab::workflow::executeJob(
        pdeJob, root, root / "pde-progress.json");
    require(pdeRun.archives.size() == 2, "two PDE trajectory archives");
    const auto pdeArchive = dynsyslab::core::readNpz(pdeRun.archives.front());
    require(pdeArchive.at("values").shape == std::vector<std::size_t>({3, 16}),
            "PDE field history shape");
    require(pdeArchive.at("space").shape == std::vector<std::size_t>({16}),
            "PDE coordinate shape");
    require(pdeArchive.metadataJson().find("\"system_family\":\"pde\"")
            != std::string::npos, "PDE family provenance");
    require(pdeArchive.metadataJson().find("\"pde_execution_mode\":\"local-partitioned\"")
            != std::string::npos, "PDE partitioned execution provenance");
    const auto pdeHierarchy = dynsyslab::core::json::load(pdeRun.manifest.wstring());
    require(pdeHierarchy.at("systems").array().front().at("system_family").string()
            == "pde", "PDE hierarchy family");
    dynsyslab::core::ExecutionContext limitedPdeExecution;
    limitedPdeExecution.limits.maximumThreads = 1;
    bool threadLimitRejected = false;
    try {
        (void)dynsyslab::workflow::executeJob(pdeJob, root,
            root / "limited-pde-progress.json", &limitedPdeExecution);
    } catch (const dynsyslab::core::ExecutionLimitExceeded&) {
        threadLimitRejected = true;
    }
    require(threadLimitRejected, "PDE partition workers respect execution limit");
    auto invalidPdeJob = pdeJob;
    invalidPdeJob.at("initial_conditions").at("mode") =
        dynsyslab::core::json::Value("uniform");
    invalidPdeJob.at("initial_conditions").at("ranges") =
        dynsyslab::core::json::parse("[-0.25,0.25]");
    invalidPdeJob.at("initial_conditions").at("seed") =
        dynsyslab::core::json::Value(-1.0);
    bool invalidSeedRejected = false;
    try {
        (void)dynsyslab::workflow::executeJob(
            invalidPdeJob, root, root / "invalid-pde-progress.json");
    } catch (const std::invalid_argument&) {
        invalidSeedRejected = true;
    }
    require(invalidSeedRejected, "PDE random IC rejects negative seed");

    const auto campaignJob = dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"native-campaign-test",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,
                            "values":[],"ranges":[]},
      "controller":{"kind":"waveform","settings":{"amplitude":0.1,"frequency":1},
                    "target":{"kind":"constant","state":[0,0]}},
      "campaign":{"enabled":true,"kind":"amplitude_sweep","mode":"cartesian",
                  "axes":[{"path":"controller.settings.amplitude","values":[0.1,0.2]}]},
      "solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},
      "archive":{"enabled":true,"data_role":"control_training"}
    })");
    const auto campaign = dynsyslab::workflow::executeCampaign(
        campaignJob, root, root / "campaign-progress.json");
    require(campaign.points.size() == 2, "two Cartesian campaign points");
    require(std::filesystem::is_regular_file(campaign.manifest), "campaign manifest");
    require(dynsyslab::core::json::load(campaign.manifest.wstring()).at("points").array().size() == 2,
            "campaign manifest point count");

    auto delayCampaignJob = delayJob;
    delayCampaignJob.at("label") = dynsyslab::core::json::Value("delay-parameter-campaign");
    delayCampaignJob.at("initial_conditions").at("count") =
        dynsyslab::core::json::Value(1.0);
    delayCampaignJob.at("initial_conditions").at("values") =
        dynsyslab::core::json::parse("[[1.2]]");
    delayCampaignJob.at("campaign") = dynsyslab::core::json::parse(R"({
      "enabled":true,"kind":"parameter_sweep","mode":"cartesian",
      "axes":[{"path":"system.parameters.beta","values":[0.2,0.3]}]
    })");
    double lastCampaignFraction = 0.0;
    bool monotonicCampaignProgress = true;
    dynsyslab::core::ExecutionContext campaignExecution;
    campaignExecution.progress = [&](const dynsyslab::core::ExecutionProgress& event) {
        if (event.stage.rfind("campaign", 0) != 0) return;
        monotonicCampaignProgress = monotonicCampaignProgress
            && event.fraction + 1.0e-12 >= lastCampaignFraction;
        lastCampaignFraction = event.fraction;
    };
    const auto delayCampaign = dynsyslab::workflow::executeJob(
        delayCampaignJob, root, root / "delay-campaign-progress.json",
        &campaignExecution);
    require(delayCampaign.kind == "campaign", "DDE campaign dispatcher");
    require(delayCampaign.archives.size() == 2, "DDE campaign archive count");
    require(monotonicCampaignProgress && lastCampaignFraction == 1.0,
            "campaign progress maps child execution monotonically");
    const auto delayCampaignManifest = dynsyslab::core::json::load(
        delayCampaign.manifest.wstring());
    require(delayCampaignManifest.at("schema").string() == "dynsyslab.campaign/v3",
            "generic campaign schema");
    const auto firstDelayHierarchy = dynsyslab::core::json::load(
        std::filesystem::u8path(delayCampaignManifest.at("points").array()[0]
            .at("hierarchy_manifest").string()).wstring());
    const auto secondDelayHierarchy = dynsyslab::core::json::load(
        std::filesystem::u8path(delayCampaignManifest.at("points").array()[1]
            .at("hierarchy_manifest").string()).wstring());
    require(firstDelayHierarchy.at("systems").array()[0].at("resolved_parameters")
            .at("beta").number() == 0.2
        && secondDelayHierarchy.at("systems").array()[0].at("resolved_parameters")
            .at("beta").number() == 0.3,
        "DDE campaign preserves distinct equation parameters");
    const auto firstDelayCampaignArchive = dynsyslab::core::readNpz(
        delayCampaign.archives.front());
    const auto firstDelayCampaignMetadata = dynsyslab::core::json::parse(
        firstDelayCampaignArchive.metadataJson());
    require(firstDelayCampaignMetadata.at("campaign_id").string()
            == delayCampaignManifest.at("id").string()
        && firstDelayCampaignMetadata.at("campaign_point_index").number() == 0.0,
        "DDE trajectory carries campaign identity");

    auto latticeCampaignJob = latticeJob;
    latticeCampaignJob.at("label") = dynsyslab::core::json::Value(
        "lattice-temperature-campaign");
    latticeCampaignJob.at("initial_conditions").at("count") =
        dynsyslab::core::json::Value(1.0);
    latticeCampaignJob.at("time").at("stop") = dynsyslab::core::json::Value(1.0);
    latticeCampaignJob.at("campaign") = dynsyslab::core::json::parse(R"({
      "enabled":true,"kind":"temperature_sweep","mode":"cartesian",
      "axes":[{"path":"system.parameters.temperature","values":[1.5,3.5]}]
    })");
    const auto latticeCampaign = dynsyslab::workflow::executeJob(
        latticeCampaignJob, root, root / "lattice-campaign-progress.json");
    require(latticeCampaign.kind == "campaign", "lattice campaign dispatcher");
    require(latticeCampaign.archives.size() == 2, "lattice campaign archive count");
    const auto latticeCampaignManifest = dynsyslab::core::json::load(
        latticeCampaign.manifest.wstring());
    const auto firstLatticeCampaignMetadata = dynsyslab::core::json::parse(
        dynsyslab::core::readNpz(latticeCampaign.archives.front()).metadataJson());
    require(firstLatticeCampaignMetadata.at("campaign_id").string()
            == latticeCampaignManifest.at("id").string()
        && firstLatticeCampaignMetadata.at("campaign_point_index").number() == 0.0,
        "lattice trajectory carries campaign identity");

    auto pdeCampaignJob = pdeJob;
    pdeCampaignJob.at("label") = dynsyslab::core::json::Value(
        "pde-diffusivity-campaign");
    pdeCampaignJob.at("initial_conditions").at("count") =
        dynsyslab::core::json::Value(1.0);
    pdeCampaignJob.at("campaign") = dynsyslab::core::json::parse(R"({
      "enabled":true,"kind":"parameter_sweep","mode":"cartesian",
      "axes":[{"path":"system.parameters.diffusivity","values":[0.001,0.002]}]
    })");
    const auto pdeCampaign = dynsyslab::workflow::executeJob(
        pdeCampaignJob, root, root / "pde-campaign-progress.json");
    require(pdeCampaign.archives.size() == 2, "PDE campaign archive count");
    const auto pdeCampaignManifest = dynsyslab::core::json::load(
        pdeCampaign.manifest.wstring());
    const auto firstPdeCampaignMetadata = dynsyslab::core::json::parse(
        dynsyslab::core::readNpz(pdeCampaign.archives.front()).metadataJson());
    require(firstPdeCampaignMetadata.at("campaign_id").string()
            == pdeCampaignManifest.at("id").string()
        && firstPdeCampaignMetadata.at("campaign_point_index").number() == 0.0,
        "PDE trajectory carries campaign identity");

    // Eigenfunction/spectrum methods train on UNFORCED data (the forced run feeds only the
    // control regression below); the pipeline refuses forced sources for them.
    auto freeJob = job;
    freeJob.at("controller") = dynsyslab::core::json::parse(R"({"kind":"none","settings":{}})");
    freeJob.at("label") = dynsyslab::core::json::Value("native-test-free");
    const auto freeResult = dynsyslab::workflow::executeSolve(freeJob, root, root / "free-progress.json");
    const std::string postJson = "{\"schema\":\"dynsyslab.compute-job/v1\","
        "\"kind\":\"postprocess\",\"label\":\"native-post-test\",\"sources\":[\""
        + jsonEscape(freeResult.archives.front().u8string()) + "\"],\"methods\":["
        "{\"name\":\"dmd\",\"enabled\":true,\"settings\":{}},"
        "{\"name\":\"hdmd\",\"enabled\":true,\"settings\":{\"matrix_size\":4,"
        "\"sampling_window\":6,\"observable\":\"position\"}},"
        "{\"name\":\"sindy\",\"enabled\":true,\"settings\":{\"polynomial_degree\":2,"
        "\"threshold\":1e-6,\"max_iterations\":2,\"compute_training_values\":true,"
        "\"coordinate_observables\":[\"position\",\"velocity\"]}},"
        "{\"name\":\"edmd\",\"enabled\":true,\"settings\":{\"polynomial_degree\":2,"
        "\"regularization\":1e-10}},"
        "{\"name\":\"havok\",\"enabled\":true,\"settings\":{\"delays\":4,"
        "\"rank\":3,\"observable\":\"position\"}},"
        "{\"name\":\"svd\",\"enabled\":true,\"settings\":{}}],"
        "\"archive\":{\"enabled\":true}}";
    const auto post = dynsyslab::workflow::executePostprocess(dynsyslab::core::json::parse(postJson), root,
                                                              root / "post-progress.json");
    require(post.archives.size() == 6, "native post-processing archives");  // DMD, H-DMD, SINDy, EDMD (kroms), HAVOK, SVD
    const auto spectrum = dynsyslab::core::readNpz(post.archives.front());
    require(spectrum.contains("discrete_eigenvalues"), "DMD discrete spectrum");
    require(spectrum.contains("continuous_eigenvalues"), "DMD continuous spectrum");
    const auto sindy = dynsyslab::core::readNpz(post.archives[2]);
    const auto& sindyReadouts = sindy;  // C belongs to the krom
    require(sindy.contains("expansion_coefficients"), "SINDy expansion coefficients");
    require(!sindy.contains("training_states"), "SINDy copies no training states");
    {
        const auto set = dynsyslab::core::readKoopmanModel(sindy);
        require(set.representation().function->coordinates().names
                    == std::vector<std::string>{"position", "velocity"}, "SINDy observable coordinates");
    }
    require(sindyReadouts.at("reconstruction_names").string()=="state"
            && sindyReadouts.at("reconstruction__state").shape[0]==2,
            "SINDy observable reconstruction C");
    const auto hdmd=dynsyslab::core::readNpz(post.archives[1]);
    require(hdmd.contains("eigenfunction_values") && hdmd.contains("reconstruction_matrix")
            && hdmd.contains("reconstruction_names"), "H-DMD krom with its readouts");
    const auto edmd=dynsyslab::core::readNpz(post.archives[3]);
    require(edmd.contains("eigenfunction_coefficients")&&edmd.contains("operator"),
            "native EDMD contract");
    const auto havok=dynsyslab::core::readNpz(post.archives[4]);
    require(havok.contains("A")&&havok.contains("B")&&havok.contains("forcing"),
            "native HAVOK contract");
    const auto controlledCopy=root / "controlled-copy.npz";
    std::filesystem::copy_file(result.archives.front(),controlledCopy,
        std::filesystem::copy_options::overwrite_existing);
    const std::string regressionJson = "{\"schema\":\"dynsyslab.compute-job/v1\","
        "\"kind\":\"postprocess\",\"label\":\"native-regression-test\",\"sources\":[\""
        + jsonEscape(result.archives.front().u8string()) + "\"],\"methods\":[{"
        "\"name\":\"control_regression\",\"enabled\":true,\"settings\":{\"krom\":\""
        + jsonEscape(post.archives[2].u8string()) + "\",\"additional_sources\":[\""
        + jsonEscape(controlledCopy.u8string())
        + "\"],\"grid_columns\":9,\"grid_rows\":7}}],"
        "\"archive\":{\"enabled\":true}}";
    const auto regression = dynsyslab::workflow::executePostprocess(
        dynsyslab::core::json::parse(regressionJson), root, root / "regression-progress.json");
    require(regression.archives.size() == 1, "native control regression archive");
    const auto regressed = dynsyslab::core::readNpz(regression.archives.front());
    require(regressed.contains("control_constant_matrix"), "constant control regression");
    require(regressed.contains("control_constant_discrete_matrix")
            &&regressed.contains("control_constant_sample_time"),
            "multistep regression preserves discrete/continuous conversion evidence");
    require(regressed.contains("control_bilinear_operators"), "bilinear control regression");
    require(regressed.at("bilinear_relative_error").shape ==
            std::vector<std::size_t>({7,9}), "control regression error map");
    require(!regressed.contains("reconstruction__state")
            && dynsyslab::core::loadKoopmanModel(regression.archives.front()).hasReconstruction("state"),
            "control regression is a child krom: it keeps its parent's readouts by reference");
    const auto regressionMetadata=dynsyslab::core::json::parse(regressed.metadataJson());
    require(regressionMetadata.at("constant_regression_method").string()=="multistep-output",
            "control regression method provenance");
    require(regressionMetadata.at("controlled_trajectory_count").number()==2.0,
            "control regression pools trajectory archives");

    auto fit = dynsyslab::core::polynomialEigenfunctions(2, {0,1}, {{1.0,0.0}},
        dynsyslab::core::Spectrum::fromContinuousRates({{-1.0,0.0}}), {0,0}, {1,1});
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing"; configuration.preset = "double_well_reference";
    const auto kromPath = root / "velocity-krom.npz";
    dynsyslab::core::KoopmanModel velocityModel(std::move(fit), configuration);
    velocityModel.setReconstruction("velocity",{{{1.0,0.0}},1,1,0.0,0.0,0.0});
    velocityModel.setReconstruction("state",{{{0.0,0.0},{1.0,0.0}},2,1,0.0,0.0,0.0});
    dynsyslab::core::saveKoopmanModel(kromPath,velocityModel);
    const std::string estimationJson = "{\"schema\":\"dynsyslab.compute-job/v1\","
        "\"kind\":\"postprocess\",\"label\":\"native-estimation-test\",\"sources\":[\""
        + jsonEscape(result.archives.front().u8string()) + "\"],\"methods\":[{"
        "\"name\":\"klqe\",\"enabled\":true,\"settings\":{\"krom\":\""
        + jsonEscape(kromPath.u8string()) + "\",\"measurement\":\"velocity\","
        "\"state_reconstruction\":\"state\",\"process_variance\":1,"
        "\"measurement_variance\":0.1,"
        "\"observability_method\":\"continuous-infinite-gramian\"}}],"
        "\"archive\":{\"enabled\":true}}";
    const auto estimation=dynsyslab::workflow::executePostprocess(
        dynsyslab::core::json::parse(estimationJson),root,root/"estimation-progress.json");
    require(estimation.archives.size()==1,"native KLQE archive");
    const auto estimated=dynsyslab::core::readNpz(estimation.archives.front());
    require(estimated.contains("observability_singular_values"),"native observability output");
    require(estimated.contains("observability_gramian"),"native observability Gramian output");
    require(estimated.metadataJson().find(
            "\"observability_method\":\"continuous-infinite-gramian\"")
            !=std::string::npos,"native observability method provenance");
    require(estimated.contains("estimator_gain") && estimated.contains("observer_poles"),
            "native LQE design output");
    require(estimated.contains("values") && estimated.contains("innovation"),
            "native KLQE run output");
    const std::string klqrJson = "{\"schema\":\"dynsyslab.compute-job/v1\","
        "\"kind\":\"solve\",\"label\":\"native-klqr-test\","
        "\"system\":{\"name\":\"duffing\",\"preset\":\"double_well_reference\",\"parameters\":{}},"
        "\"time\":{\"start\":0,\"stop\":0.02,\"output_step\":0.01,\"coordinate_system\":\"automatic\"},"
        "\"initial_conditions\":{\"mode\":\"explicit\",\"count\":1,\"seed\":0,"
        "\"values\":[[0.1,0.1]],\"ranges\":[]},\"controller\":{\"kind\":\"klqr\","
        "\"settings\":{\"krom\":\"" + jsonEscape(kromPath.u8string())
        + "\",\"Q\":1,\"R\":1,\"control_representation\":\"exact\"},"
        "\"target\":{\"kind\":\"constant\",\"state\":[0,0]}},"
        "\"solver\":{\"method\":\"RK45\",\"rtol\":1e-8,\"atol\":1e-10},"
        "\"archive\":{\"enabled\":true}}";
    const auto klqr = dynsyslab::workflow::executeSolve(
        dynsyslab::core::json::parse(klqrJson), root, root / "klqr-progress.json");
    require(klqr.archives.size() == 1, "native KLQR archive");
    const auto controlled = dynsyslab::core::readNpz(klqr.archives.front());
    require(controlled.contains("control_values"), "native KLQR control history");
    require(controlled.metadataJson().find("dimensionless system input") != std::string::npos,
            "native KLQR control coordinate metadata");
    require(controlled.metadataJson().find("\"controller_family\":\"koopman\"")
            != std::string::npos, "KLQR remains in the Koopman controller family");
    require(controlled.metadataJson().find("\"control_evaluation\":\"continuous_rhs\"")
            != std::string::npos, "KLQR retains continuous adaptive RHS evaluation");

    verifyStateControllers(root);
    std::filesystem::remove_all(root, ignored);
}

}  // namespace

int main() {
    try {
        solveAndArchive();
        std::cout << "dynsyslab_native_solve_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_native_solve_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
