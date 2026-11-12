#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/KromValidationPlan.hpp>

#include <chrono>
#include <complex>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

core::KoopmanModel constantKrom() {
    core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues = {{0.0, 0.0}};
    fit.coefficients = {{1.0, 0.0}};
    fit.exponents = {0, 0};
    fit.center = {0.0, 0.0};
    fit.scale = {1.0, 1.0};
    fit.stateDimension = 2;
    fit.dynamicalStateDimension = 2;
    fit.featureCount = 1;
    core::SystemConfiguration system;
    system.identifier = "duffing";
    system.preset = "double_well_reference";
    return core::KoopmanModel(std::move(fit), std::move(system));
}

json::Value validationJob(const std::filesystem::path& trajectory,
                          const std::filesystem::path& krom) {
    return json::Value::Object{
        {"schema", "dynsyslab.compute-job/v1"},
        {"kind", "postprocess"},
        {"label", "constant-mode-validation"},
        {"sources", json::Value::Array{trajectory.generic_u8string()}},
        {"methods", json::Value::Array{json::Value::Object{
            {"name", "koopman_validation"}, {"enabled", true},
            {"settings", json::Value::Object{
                {"krom", krom.generic_u8string()}, {"dataset_role", "test"},
                {"horizons", json::Value::Array{0.1, 0.2}}}}}}},
        {"archive", json::Value::Object{{"enabled", true}}}};
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path() / (
            "dynsyslab-krom-validation-" + std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto krom = root / "constant-krom.npz";
        std::map<std::string, core::NpyArray> evidence;
        evidence["training_time"] = core::NpyArray::fromDoubles(
            {3}, {0.0, 0.1, 0.2});
        evidence["training_states"] = core::NpyArray::fromDoubles(
            {3, 2}, {-1.0, -0.5, 0.0, 0.0, 1.0, 0.5});
        auto original=constantKrom();
        original.addDatasetReference({"forced-training.npz","training","forced","recorded-input.npz"});
        original.addDatasetReference({"forced-test.npz","test","forced",{}});
        core::savePolynomialKoopmanModel(krom, original,
            R"({"format":"dynsyslab.derived/v1","label":"Constant KROM"})",
            std::move(evidence));

        const auto trajectory = root / "test-trajectory.npz";
        core::NpzRecord source;
        source.arrays["time"] = core::NpyArray::fromDoubles(
            {4}, {0.0, 0.1, 0.2, 0.3});
        source.arrays["values"] = core::NpyArray::fromDoubles(
            {4, 2}, {-0.5, 0.0, 0.0, 0.1, 0.5, 0.0, 1.2, 0.0});
        source.arrays["metadata_json"] = core::NpyArray::fromString(
            R"({"format":"dynsyslab.trajectory/v1","label":"Validation trajectory","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{},"integration_form":"dimensionless","forcing_type":"free","controller_kind":"none"})");
        core::writeNpzAtomic(trajectory, source);

        const auto result = dynsyslab::workflow::executeJob(
            validationJob(trajectory, krom), root / "archives");
        require(result.archives.size() == 1,
                "validation did not publish exactly one archive");
        const auto archive = core::readNpz(result.archives.front());
        require(archive.at("one_step_metrics").shape ==
                    std::vector<std::size_t>({1, 6}),
                "one-step metric shape is wrong");
        require(archive.at("horizon_relative_l2").shape ==
                    std::vector<std::size_t>({1, 2}),
                "finite-horizon metric shape is wrong");
        require(archive.at("one_step_metrics").doubles()[3] == 0.0,
                "constant eigenfunction should have zero one-step defect");
        const auto metadata = json::parse(archive.metadataJson());
        require(metadata.at("method").string() == "koopman-validation"
                    && metadata.at("dataset_role").string() == "test",
                "validation metadata is incomplete");
        require(metadata.at("outside_training_samples").number() == 1.0,
                "training-domain extrapolation was not recorded");
        require(archive.at("one_step_change_defined").doubles()[0]==0.,"constant change denominator should be undefined");
        const auto linked=core::loadKoopmanModel(std::filesystem::u8path(metadata.at("validated_krom").string()));
        require(linked.datasetReferences().size()==3,"linked KROM lost training/test/validation pointers");
        require(core::loadKoopmanModel(krom).datasetReferences().size()==2,"original KROM was overwritten");
        require(linked.extract({0}).datasetReferences().size()==3,"extraction lost dataset pointers");
        require(linked.datasetReferences()[0].dynamics=="forced"&&linked.datasetReferences()[0].inputPath=="recorded-input.npz",
                "forced training metadata did not survive round trip");
        require(!core::readNpz(std::filesystem::u8path(metadata.at("validated_krom").string())).contains("training_states"),
                "linked KROM unnecessarily duplicated training data");
        const auto plan=dynsyslab::workflow::planAutonomousKoopmanValidation(linked,krom,root/"archives");
        require(plan.at("matched_count").number()==1.,"reuse-first lookup did not find the linked free data");
        auto contaminated=linked;contaminated.addDatasetReference({trajectory.generic_u8string(),"training","unforced",{}});
        require(dynsyslab::workflow::planAutonomousKoopmanValidation(contaminated,krom,root/"archives")
            .at("matched_count").number()==0.,"training data leaked into automatic held-out selection");
        const auto empty=dynsyslab::workflow::planAutonomousKoopmanValidation(original,krom,root/"empty");
        require(empty.at("status").string()=="simulation-offer"&&!empty.at("simulation_executed").boolean(),
                "missing data must produce an offer, not an automatic solve");
        auto rejects=[&](const core::NpzRecord& invalid){
            const auto path=root/"invalid-trajectory.npz";core::writeNpzAtomic(path,invalid);
            try {(void)dynsyslab::workflow::executeJob(validationJob(path,krom),root/"archives");}
            catch(const std::invalid_argument&){return;}
            throw std::runtime_error("invalid autonomous-validation dataset was accepted");
        };
        auto forced=source;forced.arrays["control_values"]=core::NpyArray::fromDoubles({4,1},{1,1,1,1});rejects(forced);
        auto different=source;auto params=json::parse(different.metadataJson());
        params.object()["resolved_system_parameters"]=json::Value::Object{{"delta",0.7}};
        different.arrays["metadata_json"]=core::NpyArray::fromString(json::dump(params));rejects(different);
        // Dataset roles and forcing are independent, including legacy imports.
        auto legacy=core::koopmanDatasetReferencesFromMetadata(json::parse(
            R"({"method":"edmdc","training_sources":[{"states":"x.npz","inputs":"u.npz"}],"validation_archives":["v.npz"]})"));
        require(legacy.size()==2&&legacy[0].dynamics=="unknown","legacy forcing classification must not be invented");
        auto sampledFit=core::KoopmanSINDyResult::fromDiscreteCoefficients(2,{0,0},{{1.,0.}},{{1.,0.}},.2);
        const auto sampledPath=root/"sampled-krom.npz";
        core::saveKoopmanModel(sampledPath,core::KoopmanModel(sampledFit,original.systemConfiguration()),"{}");
        auto fine=source;fine.arrays["time"]=core::NpyArray::fromDoubles({7},{0.,.1,.2,.3,.4,.5,.6});
        fine.arrays["values"]=core::NpyArray::fromDoubles({7,2},std::vector<double>(14,0.));
        const auto finePath=root/"fine.npz";core::writeNpzAtomic(finePath,fine);
        const auto decimated=dynsyslab::workflow::executeJob(validationJob(finePath,sampledPath),root/"archives");
        const auto decimatedMeta=json::parse(core::readNpz(decimated.archives.front()).metadataJson());
        require(decimatedMeta.at("sample_count").number()==4.&&decimatedMeta.at("validation_delta_t").number()==.2,
                "fitted-dt exact sample selection failed");
        auto invalidStep=validationJob(finePath,sampledPath);
        invalidStep.at("methods").array()[0].at("settings").object()["validation_delta_t"]=.15;
        bool stepRejected=false;
        try {(void)dynsyslab::workflow::executeJob(invalidStep,root/"archives");}
        catch(const std::invalid_argument&){stepRejected=true;}
        require(stepRejected,"incommensurate validation dt was accepted");
        std::filesystem::remove_all(root);
        std::cout << "dynsyslab_krom_validation_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        if (!root.empty()) {
            std::error_code ignored;
            std::filesystem::remove_all(root, ignored);
        }
        std::cerr << "dynsyslab_krom_validation_tests: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
