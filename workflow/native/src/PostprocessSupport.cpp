#include <cstdio>
#include <dynsyslab/core/OutputPath.hpp>
#include <cstdint>
#include <limits>
#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"
#include "ObjectPublishing.hpp"
#include "TimeSeriesSource.hpp"

#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>

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
    // Derived names append their parent's name, so chains (coordinate change ->
    // H-DMD -> SINDy -> K-ROM ...) grow past Windows' 260-character path limit.
    // Keep the original run (start) and this step (end), and a stable FNV-1a
    // hash of the whole name between them; full provenance is in the metadata.
    constexpr std::size_t limit = 112, head = 56, tail = 40;
    if (value.size() <= limit) return value;
    std::uint64_t hash = 1469598103934665603ull;
    for (unsigned char c : value) { hash ^= c; hash *= 1099511628211ull; }
    char digest[17];
    std::snprintf(digest, sizeof digest, "%016llx", static_cast<unsigned long long>(hash));
    return value.substr(0, head) + "~" + std::string(digest, 12) + "~" + value.substr(value.size() - tail);
}

std::string safeArrayKey(std::string value) {
    for (char& item : value) if (!std::isalnum(static_cast<unsigned char>(item))
        && item!='_' && item!='-') item='_';
    return value;
}

void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned, const std::string& message) {
    if (path.empty()) return;
    (void)dynsyslab::core::prepareOutputFile(path, "progress file");
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

dynsyslab::core::TimeSeries inputSeries(const std::vector<double>& time, std::vector<double> inputs,
                                        std::size_t inputDimension, bool onIntervals) {
    if (!inputDimension && inputs.empty()) return {};
    const std::size_t rows = inputDimension ? inputs.size() / inputDimension : 0;
    if (onIntervals ? rows + 1 != time.size() && !(rows == 0 && time.empty()) : rows != time.size())
        throw std::invalid_argument("input history does not match its record's samples");
    return {std::vector<double>(time.begin(), time.begin() + std::ptrdiff_t(rows)), std::move(inputs), inputDimension};
}

SourceData withSeries(const SourceData& like, std::vector<double> time, std::vector<double> values,
                      std::size_t dimension) {
    SourceData result = like;
    // Channel identities describe the old components; they survive only when the width does.
    if (dimension != like.dimension()) result.record.channels.clear();
    result.record.series = dynsyslab::core::TimeSeries(std::move(time), std::move(values), dimension);
    return result;
}

SourceData loadSource(const std::filesystem::path& path, bool archiveCoordinates) {
    const auto archive = dynsyslab::core::readNpz(path);
    const auto& states = archive.contains("values") ? archive.at("values") : archive.at("states");
    const auto& time = archive.at("time");
    if (states.shape.size() != 2 || time.shape.size() != 1 || states.shape[0] != time.shape[0])
        throw std::invalid_argument("post-processing source must be one samples-by-state trajectory");
    SourceData source{path, archive.metadataJson(), {}};
    auto times = time.doubles();
    const char* controlName = archive.contains("control_interval_values")
        ? "control_interval_values" : "control_values";
    if (archive.contains(controlName)) {
        const auto& control = archive.at(controlName);
        const bool intervalAligned = std::string(controlName) == "control_interval_values";
        if (control.shape.size() != 2
            || control.shape[0] != states.shape[0] - (intervalAligned ? 1 : 0))
            throw std::invalid_argument("controlled source input history has incompatible shape");
        source.record.inputs = inputSeries(times, control.doubles(), control.shape[1], intervalAligned);
        source.record.inputsOnIntervals = intervalAligned;
    }
    source.record.series = dynsyslab::core::TimeSeries(std::move(times), states.doubles(), states.shape[1]);
    if (archiveCoordinates) {  // analysis on the stored values; nothing is rescaled
        auto meta = source.metadataJson.empty() ? json::Value(json::Value::Object{}) : json::parse(source.metadataJson);
        meta.object()["analysis_coordinates"] = "archive";
        source.metadataJson = json::dump(meta);
    } else {
        normalizeSourceCoordinates(source);
    }
    describeRecord(source);
    return source;
}

double uniformDeltaTime(const SourceData& source) {
    if (source.time().size() < 3) throw std::invalid_argument("post-processing needs at least three samples");
    const double step = source.time()[1] - source.time()[0];
    if (!(step>0.) || !std::isfinite(step)) throw std::invalid_argument("post-processing source time must increase");
    for (std::size_t i=2;i<source.time().size();++i) {
        const double actual=source.time()[i]-source.time()[i-1];
        const double clockScale=std::max(std::abs(source.time()[i]),std::abs(source.time()[i-1]));
        const double tolerance=1e-8*step+8*std::numeric_limits<double>::epsilon()*clockScale;
        if (!(actual>0.) || !std::isfinite(actual) || std::abs(actual-step)>tolerance)
            throw std::invalid_argument("DMD-family methods require uniform archived output time");
    }
    return step;
}

std::string metadata(const SourceData& source, const std::string& method,
                     const std::string& label, const std::string& details) {
    auto result=json::parse("{"+details+"}").object();
    // A method that names its own contract (e.g. dynsyslab.dmdc/v1) keeps it; others are generic.
    if (!result.count("schema")) result["schema"]="dynsyslab.derived/v1";
    result["label"]=label;
    // The step's own name (e.g. "hdmd"), the same as data_object.produced_by.method.
    const auto stage=currentStageMethod();
    result["method"]=stage.empty()?method:stage;
    result["source_archive"]=std::filesystem::absolute(source.path).lexically_normal().generic_u8string();
    result["execution_backend"]=json::Value::Object{
        {"name","native-cpp"},{"compiled",true},{"python_required",false}};
    if (!source.metadataJson.empty()) {
        const auto sourceMetadata=json::parse(source.metadataJson);
        const char* copied[]{"system","system_label","system_preset",
            "resolved_system_parameters","parameter_units","dimensionless_groups",
            "dimensionless_group_units","equations_version","equations_fingerprint",
            "integration_form","parameter_set_key","archive_coordinate_system",
            "source_archive_coordinate_system","source_coordinate_transform","analysis_coordinates"};
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
    const auto recorded=result["method"].string();  // the step's name
    if(recorded=="hdmd" || recorded=="hankel_dmd" || recorded=="sindy" || recorded=="edmdc" || recorded=="dmdc"
        || recorded=="boundary_eigenfunctions" || recorded=="echo_history_koopman") {
        auto refs=dynsyslab::core::koopmanDatasetReferencesFromMetadata(json::Value(result));
        dynsyslab::core::mergeKoopmanDatasetReference(refs,{
            std::filesystem::absolute(source.path).lexically_normal().generic_u8string(),
            "training",datasetDynamics(source),{}});
        result["dataset_references"]=dynsyslab::core::koopmanDatasetReferencesJson(refs);
    }
    return json::dump(json::Value(std::move(result)));
}

std::filesystem::path publishObject(const std::filesystem::path& root, const SourceData& source,
                                    const std::string& kind, const std::string& method,
                                    dynsyslab::core::NpzRecord archive) {
    // Filed by what the archive is (its declared kind, core/DataObject.hpp), stamped with how it was made.
    stampDataObject(archive, kind, method, {source.path});
    const auto directory = dynsyslab::core::objectDirectory(root, kind);
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

std::filesystem::path publishKrom(const std::filesystem::path& root, const SourceData& source,
                                  const dynsyslab::core::KoopmanModel& model, const std::string& method,
                                  const std::string& metadataJson,
                                  std::map<std::string, dynsyslab::core::NpyArray> extraArrays,
                                  const std::vector<dynsyslab::core::ChannelDescriptor>& inputChannels) {
    const auto stamped=withDataObject(metadataJson,dynsyslab::core::object_kind::krom,{source.path});
    return publishObject(root,source,dynsyslab::core::object_kind::krom,method,
                         dynsyslab::core::makeKoopmanArchive(model,stamped,std::move(extraArrays),inputChannels));
}

}  // namespace dynsyslab::workflow::postprocess_detail
