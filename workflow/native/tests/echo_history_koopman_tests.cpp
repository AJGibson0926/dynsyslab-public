#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path source(const std::filesystem::path& root,
                             std::size_t run) {
    constexpr std::size_t samples = 720;
    constexpr double step = 0.01;
    constexpr double frequency = 5.0;
    std::vector<double> time(samples), states(samples * 2), echo(samples);
    const double phase = 0.45 * static_cast<double>(run);
    for (std::size_t sample = 0; sample < samples; ++sample) {
        time[sample] = step * static_cast<double>(sample);
        states[sample * 2] = std::cos(time[sample] + phase);
        states[sample * 2 + 1] = -std::sin(time[sample] + phase);
        echo[sample] = (1.0 + 0.15 * static_cast<double>(run))
            * std::sin(frequency * time[sample] + phase);
    }
    core::NpzRecord archive;
    archive.arrays["time"] = core::NpyArray::fromDoubles({samples}, time);
    archive.arrays["values"] = core::NpyArray::fromDoubles({samples, 2}, states);
    archive.arrays["echo"] = core::NpyArray::fromDoubles({samples}, echo);
    archive.arrays["metadata_json"] = core::NpyArray::fromString(
        "{\"system\":\"duffing\",\"system_preset\":\"double_well_reference\","
        "\"forcing_type\":\"free\",\"controller_kind\":\"none\","
        "\"data_role\":\"free_training\",\"integration_form\":\"dimensionless\"}");
    const auto path = root / ("training-" + std::to_string(run) + ".npz");
    core::writeNpzAtomic(path, archive);
    return path;
}

json::Value job(const std::vector<std::filesystem::path>& sources) {
    json::Value::Array encoded;
    for (const auto& path : sources) encoded.emplace_back(path.generic_u8string());
    json::Value::Object settings{
        {"observation_array","echo"},{"decimation",1.0},{"delays",80.0},
        {"hankel_rank",2.0},{"coordinate_rank",2.0},
        {"hankel_oversampling",2.0},{"hankel_subspace_iterations",6.0},
        {"polynomial_degree",1.0},{"pair_stride",2.0},
        {"fit_start",0.79},{"fit_end",7.19},
        {"minimum_continuous_imaginary",4.5},
        {"maximum_continuous_imaginary",5.5}};
    return json::Value::Object{
        {"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
        {"label","native-echo-history-test"},{"sources",std::move(encoded)},
        {"methods",json::Value::Array{json::Value::Object{
            {"name","echo_history_koopman"},{"enabled",true},
            {"settings",std::move(settings)}}}},
        {"archive",json::Value::Object{{"enabled",true}}}};
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path() / (
            "dynsyslab-echo-history-" + std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        std::vector<std::filesystem::path> sources;
        for (std::size_t run = 0; run < 3; ++run) sources.push_back(source(root, run));
        const auto result = dynsyslab::workflow::executeJob(
            job(sources), root / "archives", {});
        std::vector<std::filesystem::path> kroms;  // the stage publishes one krom
        for(const auto& path:result.archives)
            if(core::dataObjectKind(core::readNpz(path).metadataJson())==core::object_kind::krom)kroms.push_back(path);
        require(result.archives.size() == 1 && kroms.size() == 1,
                "pooled echo-history stage must publish one krom");
        const auto archive = core::readNpz(kroms.front());
        require(archive.at("left_singular_vectors").shape
                == std::vector<std::size_t>({80,2})
            && archive.at("training_trajectory_offsets").elementCount() == 3,
                "native echo-history archive contract is incomplete");
        const auto metadata = json::parse(archive.metadataJson());
        require(metadata.at("native_history_pipeline").boolean()
            && !metadata.at("execution_backend").at("python_required").boolean(true),
                "echo-history native provenance is incomplete");
        const auto model = core::loadKoopmanModel(kroms.front());
        require(model.requiresExternalCoordinates()
            && model.coordinateDimension() == 2 && model.liftedDimension() == 1,
                "history KROM did not load through the canonical model contract");
        bool rejected = false;
        const double state[]{0.0,0.0};
        try { (void)model.evaluate(state); }
        catch (const std::logic_error&) { rejected = true; }
        require(rejected && !model.evaluateCoordinates(state).empty(),
                "history KROM did not distinguish state from supplied coordinates");
        std::filesystem::remove_all(root);
        std::cout << "native echo-history workflow tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        if (!root.empty()) {
            std::error_code ignored;
            std::filesystem::remove_all(root, ignored);
        }
        std::cerr << "echo-history workflow test failure: " << error.what() << '\n';
        return 1;
    }
}
