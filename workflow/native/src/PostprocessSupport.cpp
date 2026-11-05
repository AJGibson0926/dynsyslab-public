#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"

#include <dynsyslab/core/Text.hpp>

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow::postprocess_detail {

std::string escaped(const std::string& value) {
    std::string result;
    for (const char item : value) {
        if (item == '"' || item == '\\') result.push_back('\\');
        if (item == '\n') result += "\\n"; else if (item != '\r') result.push_back(item);
    }
    return result;
}

std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t raw = std::chrono::system_clock::to_time_t(now); std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out; out << std::put_time(&utc, "%Y%m%dT%H%M%S") << 'Z'; return out.str();
}

std::string safeName(std::string value) {
    for (char& item : value) if (!std::isalnum(static_cast<unsigned char>(item))
        && item != '-' && item != '_') item = '-';
    return value;
}

std::string safeArrayKey(std::string value) {
    for (char& item : value) if (!std::isalnum(static_cast<unsigned char>(item))
        && item!='_' && item!='-') item='_';
    return value;
}

void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned, const std::string& message) {
    if (path.empty()) return;
    std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    out << "{\"schema\":\"dynsyslab.compute-progress/v1\",\"determinate\":true,\"percent\":"
        << fraction * 100.0 << ",\"estimated_remaining_seconds\":-1,\"completed_trajectories\":"
        << complete << ",\"planned_trajectories\":" << planned
        << ",\"state\":\"" << (fraction >= 1.0 ? "complete" : "running")
        << "\",\"message\":\"" << escaped(message) << "\"}";
    out.close();
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish post-processing progress");
#else
    std::filesystem::rename(temporary, path);
#endif
}

SourceData loadSource(const std::filesystem::path& path) {
    const auto archive = dynsyslab::core::readNpz(path);
    const auto& states = archive.contains("values") ? archive.at("values") : archive.at("states");
    const auto& time = archive.at("time");
    if (states.shape.size() != 2 || time.shape.size() != 1 || states.shape[0] != time.shape[0])
        throw std::invalid_argument("post-processing source must be one samples-by-state trajectory");
    std::vector<double> inputs; std::size_t inputDimension = 0;
    const char* controlName = archive.contains("control_interval_values")
        ? "control_interval_values" : "control_values";
    if (archive.contains(controlName)) {
        const auto& control = archive.at(controlName);
        const bool intervalAligned = std::string(controlName) == "control_interval_values";
        if (control.shape.size() != 2
            || control.shape[0] != states.shape[0] - (intervalAligned ? 1 : 0))
            throw std::invalid_argument("controlled source input history has incompatible shape");
        inputs = control.doubles(); inputDimension = control.shape[1];
    }
    SourceData source{path, time.doubles(), states.doubles(), states.shape[0], states.shape[1],
            std::move(inputs), inputDimension, archive.metadataJson()};
    normalizeSourceCoordinates(source);
    return source;
}

double uniformDeltaTime(const SourceData& source) {
    if (source.time.size() < 3) throw std::invalid_argument("post-processing needs at least three samples");
    const double step = source.time[1] - source.time[0];
    if (!(step > 0.0)) throw std::invalid_argument("post-processing source time must increase");
    for (std::size_t index = 2; index < source.time.size(); ++index)
        if (std::abs((source.time[index] - source.time[index - 1]) - step)
            > 1.0e-10 * std::max(1.0, std::abs(step)))
            throw std::invalid_argument("DMD-family methods require uniform archived output time");
    return step;
}

std::string metadata(const SourceData& source, const std::string& method,
                     const std::string& label, const std::string& details) {
    auto result=json::parse("{"+details+"}").object();
    result["format"]="dynsyslab.derived/v1";
    result["label"]=label;
    result["method"]=method;
    result["source_archive"]=std::filesystem::absolute(source.path).lexically_normal().generic_u8string();
    result["execution_backend"]=json::Value::Object{
        {"name","native-cpp"},{"compiled",true},{"python_required",false}};
    if (!source.metadataJson.empty()) {
        const auto sourceMetadata=json::parse(source.metadataJson);
        const char* copied[]{"system","system_label","system_preset",
            "resolved_system_parameters","parameter_units","dimensionless_groups",
            "dimensionless_group_units","equations_version","equations_fingerprint",
            "integration_form","parameter_set_key","archive_coordinate_system",
            "source_archive_coordinate_system","source_coordinate_transform"};
        for (const auto key : copied)
            if (!sourceMetadata.at(key).isNull()) result[key]=sourceMetadata.at(key);
        if (!sourceMetadata.at("system_configuration").isNull()) {
            result["system_configuration"]=sourceMetadata.at("system_configuration");
        } else if (!sourceMetadata.at("system").string().empty()) {
            json::Value::Object configuration;
            configuration["identifier"]=sourceMetadata.at("system");
            configuration["preset"]=sourceMetadata.at("system_preset");
            configuration["parameters"]=sourceMetadata.at("resolved_system_parameters");
            configuration["parameter_units"]=sourceMetadata.at("parameter_units");
            configuration["dimensionless_groups"]=sourceMetadata.at("dimensionless_groups");
            configuration["dimensionless_group_units"]=
                sourceMetadata.at("dimensionless_group_units");
            configuration["equations_version"]=sourceMetadata.at("equations_version");
            configuration["initial_state"]=sourceMetadata.at("initial_state");
            result["system_configuration"]=json::Value(std::move(configuration));
        }
    }
    if(method=="hankel-dmd" || method=="koopman-with-sindy" || method=="edmdc" || method=="dmdc"
        || method=="boundary-eigenfunctions" || method=="echo-history-koopman") {
        auto refs=dynsyslab::core::koopmanDatasetReferencesFromMetadata(json::Value(result));
        dynsyslab::core::mergeKoopmanDatasetReference(refs,{
            std::filesystem::absolute(source.path).lexically_normal().generic_u8string(),
            "training",datasetDynamics(source),{}});
        result["dataset_references"]=dynsyslab::core::koopmanDatasetReferencesJson(refs);
    }
    return json::dump(json::Value(std::move(result)));
}

std::filesystem::path publish(const std::filesystem::path& root, const SourceData& source,
                              const std::string& method, const std::string&,
                              dynsyslab::core::NpzRecord archive) {
    const auto directory = root / L"derived" / dynsyslab::core::wide(method);
    const auto stem=safeName(source.path.stem().u8string()+"__"+method+"__"+timestamp());
    for (std::size_t suffix=0;suffix<10000;++suffix) {
        const auto path=directory/dynsyslab::core::wide(stem
            +(suffix?"__"+std::to_string(suffix):"")+".npz");
        if(std::filesystem::exists(path))continue;
        try { dynsyslab::core::writeNpzAtomicNew(path,archive);return path; }
        catch(...) { if(!std::filesystem::exists(path))throw; }
    }
    throw std::runtime_error("unable to allocate an unused post-processing archive path");
}

}  // namespace dynsyslab::workflow::postprocess_detail
