#include "KromTraining.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <cctype>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dsl {
namespace {
const dynsyslab::core::NpyArray* array(const dynsyslab::core::NpzRecord& a,
    std::initializer_list<const char*> names) {
    for (const auto name : names) { const auto i = a.arrays.find(name); if (i != a.arrays.end()) return &i->second; }
    return nullptr;
}
void trainingPaths(const json::Value& value, std::vector<std::string>& paths, int depth = 0) {
    if (depth > 8) return;
    for (const auto& [key, item] : value.object()) {
        if (key == "trajectory_archives" || key == "training_archives")
            for (const auto& entry : item.array()) if (!entry.string().empty()) paths.push_back(entry.string());
        if (key == "primary_source" || key == "source_training_archive" || key == "source_archive")
            if (!item.string().empty()) paths.push_back(item.string());
        if (key != "eigenfunctions" && key != "eigenfunction_records") trainingPaths(item, paths, depth + 1);
    }
}
bool aligned(const std::vector<double>& saved, const std::vector<double>& candidate, std::size_t start) {
    if (start + candidate.size() > saved.size()) return false;
    for (std::size_t i = 0; i < candidate.size(); ++i)
        if (!std::isfinite(candidate[i]) || std::abs(saved[start+i]-candidate[i]) >
            1e-8 * std::max({1.0,std::abs(saved[start+i]),std::abs(candidate[i])})) return false;
    return true;
}
}

double kromSampleTime(const json::Value& m) {
    const auto& k = m.at("analysis_metadata").at("krom");
    for (const auto v : {m.at("delta_t").number(), k.at("delta_t").number(),
        m.at("sample_time").number()}) if (std::isfinite(v) && v > 0) return v;
    return 0;
}

dynsyslab::core::KoopmanModel loadKromEigenfunctions(const std::filesystem::path& source,
    const dynsyslab::core::NpzRecord& archive, const json::Value& metadata) {
    if (!metadata.at("system_configuration").at("identifier").string().empty() || !metadata.at("system").string().empty()) {
        return dynsyslab::core::loadKoopmanEigenfunctionModel(source);
    }
    std::vector<std::string> paths; trainingPaths(metadata,paths);
    for (const auto& encoded : paths) {
        auto p = std::filesystem::u8path(encoded); if (p.is_relative()) p = source.parent_path()/p;
        if (p.lexically_normal() == source.lexically_normal()) continue;
        dynsyslab::core::NpzRecord training;
        try { training = dynsyslab::core::readNpz(p); } catch (const std::exception&) { continue; }
        const auto meta = json::parse(training.metadataJson());
        auto name = meta.at("system").string();
        auto token = [](const std::string& text) { std::string r; for (const unsigned char c : text)
            if (std::isalnum(c)) r.push_back(static_cast<char>(std::tolower(c))); return r; };
        for (const auto& descriptor : dynsyslab::core::availableSystems())
            if (token(name) == token(descriptor.identifier) || token(name) == token(descriptor.label)) name = descriptor.identifier;
        if (name.empty()) continue;
        dynsyslab::core::SystemConfiguration config; config.identifier = name;
        const auto& legacy = meta.at("system_parameters");
        config.preset = meta.at("system_preset").string(legacy.at("preset").at("key").string());
        for (const auto* params : {&legacy.at("preset"),&legacy,&meta.at("resolved_system_parameters")})
            for (const auto& [key,value] : params->object())
                if (std::holds_alternative<double>(value.storage())) config.parameters[key] = value.number();
        const auto system = dynsyslab::core::createSystem(config);
        const auto* a = array(training,{"values","states"});
        if (!a || !archive.contains("training_states")) continue;
        const auto conversion = displayCoordinateTransform(meta,system->descriptor().stateDimension,
            system->descriptor().preferredCoordinates);
        if (conversion.coordinates != system->descriptor().preferredCoordinates ||
            !aligned(archive.at("training_states").doubles(),transformDisplayStates(*a,conversion).doubles(),0)) continue;
        return dynsyslab::core::loadKoopmanEigenfunctionModel(source,system->configuration());
    }
    throw std::runtime_error("Legacy KROM system metadata is unavailable in its explicitly recorded training archives");
}

KromTraining readKromTraining(const std::filesystem::path& source,
    const dynsyslab::core::NpzRecord& archive, const json::Value& metadata,
    const dynsyslab::core::KoopmanModel& model) {
    KromTraining result;
    // Coordinate KROMs (e.g. causal-history xi) store their samples as
    // training_coordinates; physical-state KROMs as training_states.
    const auto* states = model.requiresExternalCoordinates()
        ? array(archive, {"training_states", "training_coordinates"})
        : array(archive, {"training_states"});
    const auto expectedDimension = model.requiresExternalCoordinates()
        ? model.coordinateDimension() : model.stateDimension();
    if (!states || states->shape.size() != 2 || !states->shape[0] ||
        states->shape[1] != expectedDimension)
        throw std::runtime_error("KROM requires retained training_states aligned with its native state dimension");
    result.states = *states;
    const auto count = states->shape[0], dimension = states->shape[1];
    const auto native = states->doubles();
    const auto system = dynsyslab::core::createSystem(model.systemConfiguration());
    const auto form = system->descriptor().preferredCoordinates;
    result.metadata = metadata;
    std::vector<std::string> paths; trainingPaths(metadata, paths);
    std::set<std::filesystem::path> seen;
    std::vector<double> sourceTime;
    std::vector<std::size_t> sourceStarts;
    std::vector<std::string> sourceLabels;
    for (const auto& encoded : paths) {
        if (model.requiresExternalCoordinates()) break;
        auto p = std::filesystem::u8path(encoded);
        if (p.is_relative()) p = source.parent_path() / p;
        p = p.lexically_normal();
        if (p == source || !seen.insert(p).second) continue;
        try {
            const auto record = dynsyslab::core::readNpz(p);
            const auto* s = array(record, {"states", "values", "free_states"});
            const auto* t = array(record, {"time", "free_time"});
            if (!s || !t || s->shape.size() != 2 || s->shape[1] != dimension ||
                t->elementCount() != s->shape[0]) continue;
            const auto meta = json::parse(record.metadataJson());
            const auto conversion = displayCoordinateTransform(meta, dimension, form);
            if (conversion.coordinates != form) continue;
            const auto transformed = transformDisplayStates(*s, conversion).doubles();
            if (!aligned(native, transformed, sourceTime.size()*dimension)) continue;
            sourceStarts.push_back(sourceTime.size());
            sourceLabels.push_back(meta.at("label").string(p.stem().u8string()));
            const auto time = transformDisplayTime(t->doubles(), conversion);
            sourceTime.insert(sourceTime.end(), time.begin(), time.end());
            result.metadata = meta;
            if (sourceTime.size() == count) break;
        } catch (const std::exception&) { /* Explicit provenance may have moved; retained states remain authoritative. */ }
    }
    if (const auto* t = array(archive, {"training_time", "time"}); t && t->elementCount() == count)
        result.time = t->doubles();
    if (sourceTime.size() == count) {
        if (result.time.empty()) result.time = sourceTime;
        result.starts = sourceStarts; result.labels = sourceLabels;
    }
    if (result.time.empty()) {
        const double dt = kromSampleTime(metadata);
        result.sampleIndex = !(dt > 0);
        result.time.resize(count);
        for (std::size_t i = 0; i < count; ++i) result.time[i] = i * (dt > 0 ? dt : 1.0);
        result.warnings.push_back(dt > 0 ?
            "Training timestamps unavailable; time is relative to the first retained sample using the archived delta_t." :
            "Training timestamps unavailable; time panels use sample index.");
    }
    if (result.starts.empty()) result.starts.push_back(0);
    if (const auto* offsets = array(archive, {"training_trajectory_offsets", "trajectory_offsets"}))
        for (const auto v : offsets->doubles()) if (v >= 0 && v < count && std::floor(v) == v)
            result.starts.push_back(static_cast<std::size_t>(v));
    for (std::size_t i = 1; i < count; ++i)
        if (!(result.time[i] > result.time[i-1])) result.starts.push_back(i);
    std::sort(result.starts.begin(),result.starts.end());
    result.starts.erase(std::unique(result.starts.begin(),result.starts.end()),result.starts.end());
    if (result.labels.size() != result.starts.size()) {
        result.labels.clear();
        for (std::size_t i = 0; i < result.starts.size(); ++i)
            result.labels.push_back("Training trajectory " + std::to_string(i+1));
    }
    result.starts.push_back(count);
    if (model.requiresExternalCoordinates()) {
        result.coordinates.coordinates = "causal-history";
        result.coordinates.timeLabel = metadata.at("integration_form").string()
            == "dimensionless" ? "\\tau" : "t";
        result.coordinates.stateMultipliers.assign(dimension,1.0);
        const auto& names = model.eigenfunctions().coordinateObservableNames;
        for (std::size_t axis = 0; axis < dimension; ++axis) {
            const auto name = axis < names.size() ? names[axis]
                : "xi" + std::to_string(axis + 1);
            std::string suffix;
            if (name.rfind("xi",0) == 0) suffix = name.substr(2);
            result.coordinates.stateLabels.push_back(suffix.empty()
                ? name : "\\xi_{" + suffix + "}");
        }
        return result;
    }
    auto object = result.metadata.object();
    object["system"] = model.systemConfiguration().identifier;
    object["integration_form"] = form; object["archive_coordinate_system"] = form;
    json::Value::Object parameters;
    for (const auto& [key, value] : model.systemConfiguration().parameters) parameters[key] = value;
    object["resolved_system_parameters"] = std::move(parameters);
    const auto scaling = system->coordinates();
    if (scaling.physicalTimeAvailable) object["time_scale_to_seconds"] = scaling.secondsPerModelTime;
    result.metadata = std::move(object);
    result.coordinates = displayCoordinateTransform(result.metadata,dimension,form);
    if (result.sampleIndex) result.coordinates.timeLabel = "sample index";
    return result;
}
}
