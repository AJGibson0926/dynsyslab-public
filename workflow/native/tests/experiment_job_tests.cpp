// Experiment job mechanics: merge patch, parameters, extends/before/overrides,
// reference planning, and an executed two-stage solve -> postprocess chain.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Experiment.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
namespace json = dynsyslab::core::json;
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
template <class F> void rejects(F f) {
    bool failed = false; try { f(); } catch (const std::exception&) { failed = true; }
    require(failed, "invalid experiment accepted");
}

json::Value solveJob(double stop) {
    return json::parse(R"({"kind":"solve","system":{"name":"duffing","parameters":{"alpha":-1,"beta":1,"delta":0.2}},
      "time":{"start":0,"stop":)" + std::to_string(stop) + R"(,"output_step":0.01},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[1.0,0.0]]},
      "solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},"archive":{"enabled":true}})");
}

json::Value base() {
    json::Value pod = json::parse(R"({"kind":"postprocess","sources":["${train.archive}"],
      "methods":[{"name":"svd","settings":{"rank":"${params.rank}"}}],"archive":{"enabled":true}})");
    return json::Value::Object{
        {"schema", json::Value("dynsyslab.compute-job/v1")}, {"kind", json::Value("experiment")},
        {"label", json::Value("experiment-test")},
        {"parameters", json::parse(R"({"rank":2})")},
        {"stages", json::Value::Array{
            json::Value::Object{{"id", json::Value("train")}, {"job", solveJob(2.0)}},
            json::Value::Object{{"id", json::Value("reduce")}, {"job", pod}}}}};
}

// Training data -> coordinate change -> post-processing methods -> artefact.
// Which training data, coordinates and methods are all parameters.
json::Value plugAndPlay() {
    return json::Value::Object{
        {"schema", json::Value("dynsyslab.compute-job/v1")}, {"kind", json::Value("experiment")},
        {"label", json::Value("plug-and-play")},
        {"parameters", [] {
            auto p = json::parse(R"({"coordinates":{"kind":"identity"},"methods":["hdmd","sindy"],
              "method_settings":{
                "hdmd":{"observable":0,"matrix_size":20,"sampling_window":100,"rank":"${params.hdmd_rank}"},
                "sindy":{"polynomial_degree":3,"threshold":0,"max_iterations":1,"gamma_method":"system_rhs","reconstruction":"none"},
                "mpedmd":{"rank":4,"dictionary":{"kind":"polynomial","degree":2,"normalize":true},"assume_measure_preserving":true}},
              "hdmd_rank":6,"artefact":{"mode_selection":{"rule":"leading","count":2}}})");
            p.at("training") = solveJob(2.0);
            return p; }()},
        {"stages", json::parse(R"([
            {"id":"train","job":"${params.training}"},
            {"id":"coords","job":{"kind":"postprocess","sources":["${train.archive}"],
              "methods":[{"name":"coordinate_transform","settings":"${params.coordinates}"}]}},
            {"id":"methods","sequence":"${params.methods}","sources":["${coords.archive}"]},
            {"id":"rom","job":{"kind":"postprocess","sources":["${methods.archive}"],
              "methods":[{"name":"krom_extract","settings":"${params.artefact}"}]}}])")}};
}

json::Value stageJob(const json::Value& expanded, const std::string& id) {
    for (const auto& stage : expanded.at("stages").array()) if (stage.at("id").string() == id) return stage.at("job");
    throw std::runtime_error("missing stage " + id);
}

std::string lastMethod(const dynsyslab::workflow::RunResult& run) { return run.methods.at(run.methods.size() - 2); }

}  // namespace

int main() {
    using namespace dynsyslab::workflow;
    try {
        // RFC 7386: objects merge, null deletes, arrays replace.
        const auto patched = mergePatch(json::parse(R"({"a":{"b":1,"c":2},"d":[1,2]})"),
                                        json::parse(R"({"a":{"c":null,"e":3},"d":[9]})"));
        require(patched.at("a").at("b").number() == 1 && patched.at("a").at("c").isNull()
                && patched.at("a").at("e").number() == 3 && patched.at("d").array().size() == 1,
                "merge patch semantics");

        // Parameters substitute typed values; schema and labels are defaulted.
        const auto expanded = expandExperiment(base());
        const auto& reduce = expanded.at("stages").array()[1].at("job");
        require(reduce.at("methods").array()[0].at("settings").at("rank").number() == 2,
                "typed parameter substitution");
        require(reduce.at("schema").string() == "dynsyslab.compute-job/v1", "stage schema default");

        // extends + parameter override + before insertion + merge-patch override.
        const auto root = std::filesystem::temp_directory_path() / ("dynsyslab-experiment-test-"
            + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        { std::ofstream out(root / "base.json"); out << json::dump(base()); }
        auto derived = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"experiment",
          "label":"derived","parameters":{"rank":1},
          "stages":[{"id":"extra","before":"reduce","job":{"kind":"solve"}}],
          "overrides":{"train":{"time":{"stop":1.5}}}})");
        derived.at("extends") = json::Value((root / "base.json").generic_u8string());
        const auto variant = expandExperiment(derived);
        const auto& stages = variant.at("stages").array();
        require(stages.size() == 3 && stages[1].at("id").string() == "extra", "before insertion");
        require(stages[0].at("job").at("time").at("stop").number() == 1.5
                && stages[0].at("job").at("time").at("output_step").number() == 0.01, "override merge");
        require(stages[2].at("job").at("methods").array()[0].at("settings").at("rank").number() == 1,
                "parameter override");
        require(variant.at("extends_chain").array().size() == 1, "extends chain recorded");

        // Command-line inputs: parameter, indexed stage path, --until truncation.
        const auto inputs = applyJobInputs(base(), {{"rank", "3"}},
            {{"train.time.stop", "1.25"}, {"reduce.methods.0.settings.extra", "text"}, {"__until", "reduce"}});
        const auto& inputStages = inputs.at("stages").array();
        require(inputStages[1].at("job").at("methods").array()[0].at("settings").at("rank").number() == 3,
                "--param substitution");
        require(inputStages[0].at("job").at("time").at("stop").number() == 1.25, "--set stage path");
        require(inputStages[1].at("job").at("methods").array()[0].at("settings").at("extra").string() == "text",
                "--set array index and text fallback");
        require(inputs.at("command_inputs").at("until").string() == "reduce", "inputs recorded");
        const auto truncated = applyJobInputs(base(), {}, {{"__until", "train"}});
        require(truncated.at("stages").array().size() == 1, "--until truncation");
        rejects([&] { (void)applyJobInputs(base(), {}, {{"missing.time.stop", "1"}}); });

        // Forward/unknown references and unknown override targets fail at plan time.
        auto forward = base();
        std::swap(forward.at("stages").array()[0], forward.at("stages").array()[1]);
        rejects([&] { (void)planJob(forward); });
        auto unknown = base(); unknown.at("overrides") = json::parse(R"({"missing":{}})");
        rejects([&] { (void)expandExperiment(unknown); });
        auto undefined = base(); undefined.at("parameters") = json::Value(json::Value::Object{});
        rejects([&] { (void)expandExperiment(undefined); });
        // The old kind name is refused with a message naming the new one.
        auto renamed = base(); renamed.at("kind") = json::Value("pipeline");
        try { (void)planJob(renamed); require(false, "kind \"pipeline\" accepted"); }
        catch (const std::invalid_argument& error) {
            require(std::string(error.what()).find("\"experiment\"") != std::string::npos, "rename message");
        }

        // Execute: the postprocess stage consumes the solve stage's archive.
        const auto result = executeJob(base(), root / "archives");
        require(result.kind == "experiment" && result.archives.size() == 2, "experiment outputs");
        require(std::filesystem::is_regular_file(result.manifest), "experiment manifest");
        std::ifstream manifest(result.manifest); std::string text((std::istreambuf_iterator<char>(manifest)), {});
        const auto record = json::parse(text);
        require(record.at("schema").string() == "dynsyslab.experiment/v1"
                && record.at("stages").array().size() == 2, "manifest contents");
        const auto source = record.at("stages").array()[1].at("job").at("sources").array()[0].string();
        require(source.find("${") == std::string::npos
                && std::filesystem::path(source) == result.archives[0], "reference resolved to stage output");
        // Plug and play: the methods list expands into one stage per method;
        // SINDy is fed the preceding method's spectrum; ${methods.*} is the last.
        const auto composed = expandExperiment(plugAndPlay());
        require(stageJob(composed, "methods_1").at("methods").array()[0].at("settings").at("rank").number() == 6,
                "method settings from parameters, nested parameters resolved");
        require(stageJob(composed, "methods_2").at("methods").array()[0].at("settings").at("spectrum_archive").string()
                == "${methods_1.archive}", "spectrum chained between methods");
        (void)planExperiment(plugAndPlay());
        const auto dotted = applyJobInputs(plugAndPlay(), {{"method_settings.hdmd.rank", "4"}, {"training.time.stop", "1.5"}}, {});
        require(stageJob(dotted, "methods_1").at("methods").array()[0].at("settings").at("rank").number() == 4
                && stageJob(dotted, "train").at("time").at("stop").number() == 1.5, "dotted --param paths");
        require(applyJobInputs(plugAndPlay(), {}, {{"__until", "methods"}}).at("stages").array().size() == 4,
                "--until a sequence keeps all its methods");
        require(applyJobInputs(plugAndPlay(), {}, {{"__refresh", "methods"}}).at("refresh_stages").array().size() == 2,
                "--refresh a sequence refreshes every method");
        const auto viaSindy = executeJob(plugAndPlay(), root / "plug");
        require(viaSindy.methods.back() == "rom:krom_extract", "K-ROM from H-DMD + SINDy");
        // Same experiment, different methods: K-ROM straight from H-DMD, and from mpEDMD.
        for (const auto* methods : {R"(["hdmd"])", R"(["mpedmd"])"}) {
            const auto swapped = executeJob(applyJobInputs(plugAndPlay(), {{"methods", methods}}, {}), root / "plug");
            require(swapped.methods.back() == "rom:krom_extract" && lastMethod(swapped).rfind("methods_1:", 0) == 0,
                    "swapped post-processing method feeds the same artefact stage");
            (void)dynsyslab::core::loadKoopmanModel(swapped.archives.back());
        }
        // Different training data: an existing archive instead of a solve.
        const auto existing = json::Value(json::Value::Object{{"kind", json::Value("archive")},
            {"archives", json::Value(json::Value::Array{json::Value(viaSindy.archives.front().generic_u8string())})}});
        const auto reusedTraining = executeJob(applyJobInputs(plugAndPlay(), {{"training", json::dump(existing)}}, {}),
                                               root / "plug");
        require(reusedTraining.archives.front() == viaSindy.archives.front(), "archive stage passes its file through");
        // Nested parameter fields: ${params.NAME.FIELD}.
        {
            auto nested = base();
            nested.at("parameters") = json::parse(R"({"rank":2,"B":{"representation":"constant","depth":{"value":3}}})");
            nested.at("stages").array()[1].at("job").at("methods").array()[0].at("settings") =
                json::parse(R"({"rank":"${params.B.depth.value}","name":"${params.B.representation}"})");
            const auto resolved = expandExperiment(nested);
            const auto& s = resolved.at("stages").array()[1].at("job").at("methods").array()[0].at("settings");
            require(s.at("rank").number() == 3 && s.at("name").string() == "constant", "nested parameter fields");
            auto missing = nested;
            missing.at("stages").array()[1].at("job").at("methods").array()[0].at("settings") =
                json::parse(R"({"rank":"${params.B.nothing}"})");
            rejects([&] { (void)expandExperiment(missing); });
            // A parameter present with the value null is defined (toggles such as
            // control_limit or target_ramp are switched off with null); only a missing key is an error.
            auto toggled = nested;
            toggled.at("parameters").object()["control_limit"] = json::Value();
            toggled.at("stages").array()[1].at("job").at("methods").array()[0].at("settings") =
                json::parse(R"({"rank":2,"limit":"${params.control_limit}"})");
            const auto off = expandExperiment(toggled);
            require(off.at("stages").array()[1].at("job").at("methods").array()[0].at("settings").at("limit").isNull(),
                    "null experiment parameter resolves to null");
        }
        // Training data roles: eigenfunction methods refuse forced data; DMDc accepts it.
        {
            auto forced = solveJob(4.0);
            forced.at("schema") = json::Value("dynsyslab.compute-job/v1");
            forced.at("controller") = json::parse(R"({"kind":"random_pulse","settings":{"model_amplitude":0.5,"model_pulse_duration":0.2,"seed":3}})");
            const auto forcedRun = executeJob(forced, root / "roles");
            const auto source = json::Value(forcedRun.archives.front().generic_u8string());
            auto hdmd = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"roles",
              "methods":[{"name":"hdmd","settings":{"observable":0,"matrix_size":20,"sampling_window":100,"rank":4}}]})");
            hdmd.at("sources") = json::Value(json::Value::Array{source});
            rejects([&] { (void)executeJob(hdmd, root / "roles"); });
            auto dmdc = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"roles",
              "methods":[{"name":"dmdc","settings":{}}]})");
            dmdc.at("sources") = json::Value(json::Value::Array{source});
            (void)executeJob(dmdc, root / "roles");
        }
        // KLQE pieces: joint H-DMD -> observable_energy extraction with an intercept refit ->
        // estimator-only solve whose traces are archived.
        {
            auto trainJob = solveJob(30.0); trainJob.at("schema") = json::Value("dynsyslab.compute-job/v1");
            const auto train = executeJob(trainJob, root / "klqe");
            auto hdmd = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"klqe",
              "methods":[{"name":"hdmd","settings":{"observables":["position","velocity"],"matrix_size":200,
                "sampling_window":200,"snapshot_stride":2,"window_alignment":"end","rank":0,"time_window":[5,30],
                "reconstruction_observables":["state"]}}]})");
            hdmd.at("sources") = json::Value(json::Value::Array{json::Value(train.archives.front().generic_u8string())});
            const auto fitted = executeJob(hdmd, root / "klqe");
            auto extract = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"klqe",
              "methods":[{"name":"krom_extract","settings":{
                "mode_selection":{"rule":"observable_energy","reconstructions":["position","velocity"],
                  "thresholds":[0.99,0.99],"caps":[6,6],"max_real_part":0,"exclude_constant":true},
                "refit_reconstructions":{"observables":["state","position"],"intercept":true}}}]})");
            // H-DMD publishes its krom, readouts included.
            const auto& hdmdKrom = fitted.archives.back();
            extract.at("sources") = json::Value(json::Value::Array{json::Value(hdmdKrom.generic_u8string())});
            const auto krom = executeJob(extract, root / "klqe").archives.back();
            const auto model = dynsyslab::core::loadKoopmanModel(krom);
            require(model.hasReconstruction("state") && model.reconstruction("state").offset.size() == 2,
                    "refit reconstruction keeps its intercept through the archive");
            auto estimate = solveJob(5.0);
            estimate.at("schema") = json::Value("dynsyslab.compute-job/v1");
            estimate.at("estimator") = json::parse(R"({"kind":"klqe","measurement":"position","output":"state",
              "process_covariance":1,"measurement_covariance":1e-4,"sample_time":0.01})");
            estimate.at("estimator").at("krom") = json::Value(krom.generic_u8string());
            const auto run = executeJob(estimate, root / "klqe");
            const auto archive = dynsyslab::core::readNpz(run.archives.front());
            require(archive.contains("estimator_time") && archive.contains("state_estimate")
                    && archive.contains("innovation") && !archive.contains("control_values"),
                    "estimator-only solve archives estimator traces and stays unforced");
        }
        std::error_code ignored; std::filesystem::remove_all(root, ignored);  // best effort on Windows
        std::cout << "experiment job tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "experiment job tests failed: " << error.what() << '\n';
        return 1;
    }
}
