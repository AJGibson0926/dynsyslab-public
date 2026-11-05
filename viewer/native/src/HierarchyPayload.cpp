#include "HierarchyPayload.hpp"
#include "HierarchyCatalog.hpp"
#include "DisplayCoordinates.hpp"
#include "ObservablePayload.hpp"
#include "ArchiveSystem.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <limits>
#include <set>
#include <stdexcept>

namespace dsl {
namespace {

namespace json = dynsyslab::core::json;
using dynsyslab::core::NpyArray;

constexpr std::size_t kMaximumTrajectories = 96;
constexpr std::size_t kSamplesPerTrajectory = 100000;

struct Trajectory {
    std::string label;
    std::filesystem::path source;
    std::vector<double> time;
    std::vector<double> states;
    std::vector<double> controls;
    std::vector<double> analytic;
    json::Value metadata;
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::size_t controlColumns = 0;
    double campaignPoint = 0.0;
};

struct CampaignAxis {
    std::string label{"campaign point"};
    std::map<std::size_t, double> coordinates;
};

const NpyArray* select(const dynsyslab::core::NpzRecord& record,
                       std::initializer_list<const char*> names) {
    for (const char* name : names) {
        const auto found = record.arrays.find(name);
        if (found != record.arrays.end()) return &found->second;
    }
    return nullptr;
}

std::vector<std::size_t> sampleIndices(std::size_t count) {
    const std::size_t samples = std::min(count, kSamplesPerTrajectory);
    std::vector<std::size_t> result(samples);
    for (std::size_t index = 0; index < samples; ++index)
        result[index] = samples == 1 ? 0 : index * (count - 1) / (samples - 1);
    return result;
}

std::vector<double> column(const Trajectory& trajectory, std::size_t selected,
                           const std::vector<std::size_t>& indices) {
    std::vector<double> result;
    if (selected >= trajectory.columns) return result;
    result.reserve(indices.size());
    for (const auto row : indices)
        result.push_back(trajectory.states[row * trajectory.columns + selected]);
    return result;
}

std::vector<double> controlColumn(const Trajectory& trajectory, std::size_t selected,
                                  const std::vector<std::size_t>& indices) {
    std::vector<double> result;
    if (selected >= trajectory.controlColumns) return result;
    result.reserve(indices.size());
    for (const auto row : indices)
        result.push_back(trajectory.controls[row * trajectory.controlColumns + selected]);
    return result;
}

std::vector<double> sampledTime(const Trajectory& trajectory,
                                const std::vector<std::size_t>& indices) {
    std::vector<double> result; result.reserve(indices.size());
    for (const auto row : indices) result.push_back(trajectory.time[row]);
    return result;
}

std::string color(std::size_t index) {
    const double hue = std::fmod(0.618033988749895 * static_cast<double>(index), 1.0);
    const double scaled = hue * 6.0;
    const int sector = static_cast<int>(std::floor(scaled));
    const double fraction = scaled - sector;
    const double value = 0.92, saturation = 0.67;
    const double low = value * (1.0 - saturation);
    const double falling = value * (1.0 - saturation * fraction);
    const double rising = value * (1.0 - saturation * (1.0 - fraction));
    double red = value, green = rising, blue = low;
    if (sector == 1) { red = falling; green = value; }
    else if (sector == 2) { red = low; green = value; blue = rising; }
    else if (sector == 3) { red = low; green = falling; blue = value; }
    else if (sector == 4) { red = rising; green = low; blue = value; }
    else if (sector >= 5) { red = value; green = low; blue = falling; }
    char encoded[8]{};
    std::snprintf(encoded, sizeof(encoded), "#%02x%02x%02x",
        static_cast<int>(red * 255.0), static_cast<int>(green * 255.0),
        static_cast<int>(blue * 255.0));
    return encoded;
}

Trace trace(const Trajectory& trajectory, std::size_t index,
            std::vector<double> x, std::vector<double> y) {
    Trace result; result.label = trajectory.label; result.role = "trajectory";
    result.color = color(index); result.x = std::move(x); result.y = std::move(y);
    return result;
}

void appendArchivedTrajectory(const json::Value& archived,
                              const std::filesystem::path& projectRoot,
                              const std::filesystem::path& archiveRoot,
                              const std::filesystem::path& manifest,
                              std::vector<Trajectory>& trajectories,
                              std::vector<std::string>& warnings) {
    if (trajectories.size() >= kMaximumTrajectories) return;
    const auto source = resolveHierarchySource(projectRoot, archiveRoot,
        manifest, archived.at("solution_archive").string());
    try {
        const auto record = dynsyslab::core::readNpz(source);
        const auto* states = select(record,
            {"values", "states", "controlled_states", "free_training_values"});
        if (!states || states->shape.size() != 2 || states->shape[0] == 0) return;
        Trajectory value;
        value.label = archived.at("label").string(source.stem().u8string());
        value.source = source; value.rows = states->shape[0]; value.columns = states->shape[1];
        json::Value metadata;
        try { metadata = json::parse(record.metadataJson()); } catch (...) {}
        value.metadata = metadata;
        const auto coordinateTransform = displayCoordinateTransform(metadata, value.columns);
        value.states = transformDisplayStates(*states, coordinateTransform).doubles();
        if(record.contains("analytic_values") && record.at("analytic_values").shape==states->shape) {
            value.analytic=transformDisplayStates(record.at("analytic_values"),coordinateTransform).doubles();
            if(record.contains("analytic_valid")) {
                const auto valid=record.at("analytic_valid").doubles();
                if(valid.size()!=value.rows)value.analytic.clear();
                else for(size_t row=0;row<valid.size();++row)if(valid[row]!=1)
                    for(size_t col=0;col<value.columns;++col)value.analytic[row*value.columns+col]=std::numeric_limits<double>::quiet_NaN();
            }
        }
        const auto* time = select(record, {"time", "controlled_time", "free_time"});
        if (time && time->elementCount() >= value.rows)
            value.time = transformDisplayTime(time->doubles(), coordinateTransform);
        else {
            value.time.resize(value.rows);
            for (std::size_t row = 0; row < value.rows; ++row)
                value.time[row] = static_cast<double>(row);
        }
        if (const auto* controls = select(record, {"control_values", "controls",
                "controlled_controls", "controlled_controls_sample_aligned"})) {
            value.controls = controls->doubles();
            value.controlColumns = controls->shape.size() == 2 ? controls->shape[1] : 1;
        }
        if (!coordinateTransform.warning.empty())
            warnings.push_back(value.label + ": " + coordinateTransform.warning);
        value.campaignPoint = metadata.at("campaign_point_index").number(
            static_cast<double>(trajectories.size()));
        trajectories.push_back(std::move(value));
    } catch (...) { warnings.push_back("Could not load " + source.u8string()); }
}

void appendManifest(const CatalogItem& item, const std::filesystem::path& projectRoot,
                    const std::filesystem::path& archiveRoot,
                    const std::filesystem::path& manifest,
                    std::vector<Trajectory>& trajectories,
                    std::vector<std::string>& stateNames,
                    std::string& coordinates, std::vector<std::string>& warnings) {
    json::Value root;
    try { root = json::load(manifest.wstring()); }
    catch (...) { warnings.push_back("Could not read hierarchy " + manifest.u8string()); return; }
    if (root.at("format").string() != "dynsyslab.hierarchy/v2") return;
    if (item.kind == "forcing_campaign" &&
        root.at("campaign_id").string() != item.properties.at("campaign id")) return;
    if (item.kind == "batch" && root.at("id").string() != item.id) return;
    for (const auto& system : root.at("systems").array()) {
        if (item.kind == "system_run" && system.at("id").string() != item.id) continue;
        const auto systemRun = item.properties.find("system run id");
        if (systemRun != item.properties.end() &&
            system.at("id").string() != systemRun->second) continue;
        const auto key = system.at("parameter_set_key").string();
        const auto expected = item.properties.find("parameter set key");
        if (expected != item.properties.end() && key != expected->second) continue;
        if (stateNames.empty())
            for (const auto& name : system.at("state_names").array())
                stateNames.push_back(name.string());
        if (coordinates.empty()) coordinates = system.at("coordinate_system").string();
        for (const auto& cluster : system.at("clusters").array()) {
            if (item.kind == "cluster" && cluster.at("id").string() != item.id) continue;
            for (const auto& archived : cluster.at("trajectories").array()) {
                appendArchivedTrajectory(archived, projectRoot, archiveRoot, manifest,
                    trajectories, warnings);
            }
        }
    }
}

void appendLegacyManifest(const CatalogItem& item,
                          const std::filesystem::path& projectRoot,
                          const std::filesystem::path& archiveRoot,
                          const std::filesystem::path& manifest,
                          std::vector<Trajectory>& trajectories,
                          std::vector<std::string>& stateNames,
                          std::string& coordinates, std::vector<std::string>& warnings) {
    json::Value root;
    try {
        const auto archive = dynsyslab::core::readNpzIndex(manifest);
        root = json::parse(archive.metadataJson);
    } catch (...) {
        warnings.push_back("Could not read legacy hierarchy " + manifest.u8string());
        return;
    }
    if (root.at("format").string() != "dynsyslab.hierarchy/v1") return;
    const auto requestedSystem = item.properties.find("legacy system index");
    const auto requestedCluster = item.properties.find("legacy cluster index");
    std::size_t systemIndex = 0;
    for (const auto& system : root.at("systems").array()) {
        const auto currentSystem = systemIndex++;
        if (requestedSystem != item.properties.end() &&
            std::to_string(currentSystem) != requestedSystem->second) continue;
        const auto identifier = system.at("dynamical_system").at("name").string();
        if (stateNames.empty())
            for (const auto& descriptor : dynsyslab::core::availableSystems())
                if (descriptor.identifier == identifier ||
                    (identifier == "rayleigh-plesset" && descriptor.identifier == "rpe")) {
                    stateNames = descriptor.stateLabels; break;
                }
        std::size_t clusterIndex = 0;
        for (const auto& cluster : system.at("clusters").array()) {
            const auto currentCluster = clusterIndex++;
            if (requestedCluster != item.properties.end() &&
                std::to_string(currentCluster) != requestedCluster->second) continue;
            for (const auto& archived : cluster.at("trajectories").array()) {
                if (archived.at("metadata").at("status").string() != "complete") continue;
                if (coordinates.empty()) coordinates = archived.at("metadata")
                    .at("solver_metadata").at("integration_form").string();
                appendArchivedTrajectory(archived, projectRoot, archiveRoot, manifest,
                    trajectories, warnings);
            }
        }
    }
}

std::vector<Trajectory> loadTrajectories(
    const CatalogItem& item, const std::filesystem::path& projectRoot,
    const std::filesystem::path& archiveRoot, std::vector<std::string>& stateNames,
    std::string& coordinates, std::vector<std::string>& warnings) {
    std::vector<Trajectory> result;
    const auto source = std::filesystem::u8path(item.source);
    if (source.extension() == ".npz") {
        appendLegacyManifest(item, projectRoot, archiveRoot, source, result,
            stateNames, coordinates, warnings);
    } else if (item.kind == "forcing_campaign") {
        const auto root = archiveRoot / "hierarchies";
        for (const auto& entry : std::filesystem::recursive_directory_iterator(
                 root, std::filesystem::directory_options::skip_permission_denied)) {
            if (entry.is_regular_file() && entry.path().extension() == ".json")
                appendManifest(item, projectRoot, archiveRoot, entry.path(), result,
                    stateNames, coordinates, warnings);
            if (result.size() >= kMaximumTrajectories) break;
        }
    } else appendManifest(item, projectRoot, archiveRoot, source, result,
        stateNames, coordinates, warnings);
    if (result.size() >= kMaximumTrajectories)
        warnings.push_back("Aggregate display limited to the first 96 trajectories.");
    return result;
}

CampaignAxis campaignAxis(const CatalogItem& item,
                          const std::filesystem::path& archiveRoot) {
    CampaignAxis result;
    const auto identity = item.properties.find("campaign id");
    if (identity == item.properties.end()) return result;
    const auto root = archiveRoot / "campaigns";
    if (!std::filesystem::is_directory(root)) return result;
    for (const auto& entry : std::filesystem::directory_iterator(root)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".json") continue;
        try {
            const auto campaign = json::load(entry.path().wstring());
            if (campaign.at("format").string() != "dynsyslab.campaign/v3" ||
                campaign.at("id").string() != identity->second) continue;
            const auto& axes = campaign.at("axes").array();
            const auto& points = campaign.at("points").array();
            for (std::size_t axis = 0; axis < axes.size(); ++axis) {
                std::map<std::size_t, double> coordinates;
                for (const auto& point : points) {
                    const auto index = static_cast<std::size_t>(point.at("index").number());
                    const auto& values = point.at("coordinates").array();
                    if (axis >= values.size()) { coordinates.clear(); break; }
                    const double value = values[axis].number(
                        std::numeric_limits<double>::quiet_NaN());
                    if (!std::isfinite(value)) { coordinates.clear(); break; }
                    coordinates[index] = value;
                }
                if (!coordinates.empty()) {
                    result.label = axes[axis].at("path").string("campaign coordinate");
                    result.coordinates = std::move(coordinates); return result;
                }
            }
        } catch (...) { /* Continue to another campaign manifest. */ }
    }
    return result;
}

}  // namespace

std::vector<std::filesystem::path> hierarchyTrajectorySources(const CatalogItem& item,
    const std::filesystem::path& projectRoot,const std::filesystem::path& archiveRoot) {
    std::vector<std::string> names,warnings;std::string coordinates;
    std::vector<std::filesystem::path> sources;
    for(const auto& t:loadTrajectories(item,projectRoot,archiveRoot,names,coordinates,warnings)) sources.push_back(t.source);
    return sources;
}

Payload buildHierarchyPayload(const CatalogItem& item, const std::string& view,
                              const std::filesystem::path& projectRoot,
                              const std::filesystem::path& archiveRoot,
                              const std::optional<PhaseSpaceDomain>& phaseDomain) {
    Payload payload{item.source, view, item.label, 2};
    std::vector<std::string> names; std::string coordinates;
    auto trajectories = loadTrajectories(item, projectRoot, archiveRoot,
        names, coordinates, payload.warnings);
    if (trajectories.empty()) throw std::runtime_error("hierarchy contains no readable trajectories");
    if (view == "observable_time" || view == "observable_workspace") {
        for (std::size_t index = 0; index < trajectories.size(); ++index) {
            const auto record = dynsyslab::core::readNpz(trajectories[index].source);
            const auto metadata = json::parse(record.metadataJson());
            auto observables = buildObservablePayload(trajectories[index].source, view, record, metadata, phaseDomain);
            // Aggregate Main views show the selected solutions, without each solution's auxiliary overlays.
            for (auto& figure : observables.figures) {
                figure.traces.erase(std::remove_if(figure.traces.begin(), figure.traces.end(),
                    [](const auto& trace) { return trace.role != "data"; }), figure.traces.end());
            }
            appendObservablePayload(payload, std::move(observables), trajectories[index].label, color(index));
        }
        return payload;
    }
    const std::size_t dimensions = trajectories.front().columns;
    while (names.size() < dimensions) names.push_back("state " + std::to_string(names.size()));
    const std::string suffix = coordinates == "dimensionless" ? " (dimensionless)" : "";
    const auto axis = campaignAxis(item, archiveRoot);
    for (auto& trajectory : trajectories) {
        const auto found = axis.coordinates.find(
            static_cast<std::size_t>(trajectory.campaignPoint));
        if (found != axis.coordinates.end()) trajectory.campaignPoint = found->second;
    }
    auto analyticTrace = [](const Trajectory& t,const std::vector<std::size_t>& samples,size_t component,bool phase) {
        Trace line;line.role="analytic_solution";line.label=t.label+" — Analytic";line.color="#7b6d9c";line.lineStyle="dash";
        if(t.analytic.size()!=t.states.size())return line;
        for(const auto row:samples){line.x.push_back(phase?t.analytic[row*t.columns]:t.time[row]);line.y.push_back(t.analytic[row*t.columns+component]);}
        return line;
    };
    if (view == "ensemble_phase") {
        if (dimensions < 2) throw std::runtime_error("ensemble phase view requires two state components");
        Figure figure{"phase", "Trajectory ensemble", names[0], names[1], false};
        figure.coordinates = coordinates;
        for (std::size_t index = 0; index < trajectories.size(); ++index) {
            const auto samples = sampleIndices(trajectories[index].rows);
            figure.traces.push_back(trace(trajectories[index], index,
                column(trajectories[index], 0, samples),
                column(trajectories[index], 1, samples)));
            auto analytic=analyticTrace(trajectories[index],samples,1,true);
            if(!analytic.x.empty())figure.traces.push_back(std::move(analytic));
            appendVerifiedEquilibria(figure,trajectories[index].metadata);
        }
        payload.figures.push_back(std::move(figure)); return payload;
    }
    if (view == "ensemble_control") {
        std::size_t controls = 0;
        for (const auto& trajectory : trajectories)
            controls = std::max(controls, trajectory.controlColumns);
        if (!controls) throw std::runtime_error("selected hierarchy has no control histories");
        for (std::size_t component = 0; component < controls; ++component) {
            Figure figure{"time", "Control component " + std::to_string(component),
                          "time" + suffix, "u" + std::to_string(component) + suffix, false};
            for (std::size_t index = 0; index < trajectories.size(); ++index) {
                if (component >= trajectories[index].controlColumns) continue;
                const auto samples = sampleIndices(std::min(trajectories[index].rows,
                    trajectories[index].controls.size() / trajectories[index].controlColumns));
                figure.traces.push_back(trace(trajectories[index], index,
                    sampledTime(trajectories[index], samples),
                    controlColumn(trajectories[index], component, samples)));
            }
            payload.figures.push_back(std::move(figure));
        }
        return payload;
    }
    if (view == "campaign_summary") {
        Figure amplitude{"time", "State-0 response amplitude", axis.label,
                         "half peak-to-peak" + suffix, false};
        Figure terminal{"time", "Terminal state 0", axis.label, names[0] + suffix, false};
        Trace amplitudeTrace; amplitudeTrace.label = "amplitude"; amplitudeTrace.marker = true;
        Trace terminalTrace; terminalTrace.label = "terminal"; terminalTrace.marker = true;
        for (const auto& trajectory : trajectories) {
            const auto values = column(trajectory, 0, sampleIndices(trajectory.rows));
            const auto bounds = std::minmax_element(values.begin(), values.end());
            amplitudeTrace.x.push_back(trajectory.campaignPoint);
            amplitudeTrace.y.push_back(0.5 * (*bounds.second - *bounds.first));
            terminalTrace.x.push_back(trajectory.campaignPoint);
            terminalTrace.y.push_back(trajectory.states[(trajectory.rows - 1) * trajectory.columns]);
        }
        amplitude.traces.push_back(std::move(amplitudeTrace));
        terminal.traces.push_back(std::move(terminalTrace));
        payload.figures.push_back(std::move(amplitude)); payload.figures.push_back(std::move(terminal));
        return payload;
    }
    if (view == "campaign_bifurcation") {
        Figure figure{"campaign_bifurcation", "Asymptotic state-0 response", axis.label,
                      names[0] + suffix, false};
        for (std::size_t index = 0; index < trajectories.size(); ++index) {
            const auto& trajectory = trajectories[index];
            Trace points; points.label = trajectory.label; points.role = "campaign_tail";
            points.color = color(index); points.marker = true;
            const std::size_t tail = std::max<std::size_t>(1, trajectory.rows / 5);
            const std::size_t stride = std::max<std::size_t>(1, tail / 200);
            for (std::size_t row = trajectory.rows - tail; row < trajectory.rows; row += stride) {
                points.x.push_back(trajectory.campaignPoint);
                points.y.push_back(trajectory.states[row * trajectory.columns]);
            }
            figure.traces.push_back(std::move(points));
        }
        payload.figures.push_back(std::move(figure)); return payload;
    }
    for (std::size_t component = 0; component < dimensions; ++component) {
        Figure figure{"time", names[component] + " ensemble", "time" + suffix,
                      names[component] + suffix, false};
        for (std::size_t index = 0; index < trajectories.size(); ++index) {
            if (component >= trajectories[index].columns) continue;
            const auto samples = sampleIndices(trajectories[index].rows);
            figure.traces.push_back(trace(trajectories[index], index,
                sampledTime(trajectories[index], samples),
                column(trajectories[index], component, samples)));
            auto analytic=analyticTrace(trajectories[index],samples,component,false);
            if(!analytic.x.empty())figure.traces.push_back(std::move(analytic));
        }
        payload.figures.push_back(std::move(figure));
    }
    return payload;
}

}  // namespace dsl
