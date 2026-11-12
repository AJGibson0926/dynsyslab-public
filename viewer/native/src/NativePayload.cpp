#include "dynsyslab_viewer/NativePayload.hpp"
#include "KromPayload.hpp"
#include "ObservablePayload.hpp"
#include "LatticePayload.hpp"
#include "KoopmanValidationPayload.hpp"
#include "SindyDynamicsPayload.hpp"
#include "EigenvalueOptimizationPayload.hpp"
#include "KromControlCampaignPayload.hpp"
#include "FourierPayload.hpp"
#include "DisplayCoordinates.hpp"
#include "TrajectoryPresentation.hpp"
#include "AnalyticOverlay.hpp"
#include "FieldPayload.hpp"

#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <iterator>
#include <stdexcept>

namespace dsl {
namespace {

using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;

constexpr std::size_t kDisplaySamples = 2400;

const NpyArray* select(const NpzRecord& record,
                       std::initializer_list<const char*> names,
                       std::size_t dimensions = 0) {
    for (const char* name : names) {
        const auto found = record.arrays.find(name);
        if (found != record.arrays.end() && (!dimensions || found->second.shape.size() == dimensions))
            return &found->second;
    }
    return nullptr;
}

std::vector<std::size_t> displayIndices(std::size_t count, std::size_t maximum = kDisplaySamples) {
    if (count <= maximum) {
        std::vector<std::size_t> result(count);
        for (std::size_t index = 0; index < count; ++index) result[index] = index;
        return result;
    }
    std::vector<std::size_t> result;
    result.reserve(maximum);
    for (std::size_t index = 0; index < maximum; ++index)
        result.push_back(index * (count - 1) / (maximum - 1));
    return result;
}

std::vector<double> sampledVector(const NpyArray* array, std::size_t count,
                                  const std::vector<std::size_t>& indices) {
    std::vector<double> values;
    if (array) values = array->doubles();
    if (values.size() < count) {
        values.resize(count);
        for (std::size_t index = 0; index < count; ++index) values[index] = static_cast<double>(index);
    }
    std::vector<double> result;
    result.reserve(indices.size());
    for (const auto index : indices) result.push_back(values[index]);
    return result;
}

std::vector<double> sampledColumn(const NpyArray& array, std::size_t column,
                                  const std::vector<std::size_t>& indices) {
    if (array.shape.size() != 2 || column >= array.shape[1]) return {};
    const auto values = array.doubles();
    const std::size_t columns = array.shape[1];
    std::vector<double> result;
    result.reserve(indices.size());
    for (const auto row : indices) result.push_back(values[row * columns + column]);
    return result;
}

std::string metadataString(const json::Value& metadata, const char* key,
                           const std::string& fallback = {}) {
    return metadata.at(key).string(fallback);
}

std::string controlAxisSuffix(const json::Value& metadata) {
    const std::string coordinates = metadataString(metadata,
        "control_coordinate_system", metadataString(metadata, "integration_form",
        metadataString(metadata, "coordinate_system")));
    return coordinates == "dimensionless" ? " (dimensionless)" : std::string{};
}

std::string archiveLabel(const json::Value& metadata, const std::filesystem::path& source) {
    return metadataString(metadata, "label", source.stem().u8string());
}

Trace line(std::string label, const std::vector<double>& x, std::vector<double> y,
           std::string role = "data", std::string color = {}, std::string style = "solid") {
    Trace result;
    result.label = std::move(label); result.x = x; result.y = std::move(y);
    result.role = std::move(role); result.color = std::move(color);
    result.lineStyle = std::move(style);
    return result;
}

Figure componentFigure(const std::string& label, const std::string& suffix,
                       const std::string& timeLabel,
                       const std::vector<double>& time, const NpyArray& states,
                       std::size_t component, const std::vector<std::size_t>& indices,
                       const NpyArray* target, const NpyArray* truth) {
    Figure figure{"time", label + " versus time", timeLabel, label + suffix, false};
    figure.traces.push_back(line(label, time, sampledColumn(states, component, indices)));
    if (target && target->shape.size() == 2 && component < target->shape[1])
        figure.traces.push_back(line("target", time, sampledColumn(*target, component, indices),
                                     "control_target", "#ef4444", "dash"));
    if (truth && truth->shape.size()==2 && component<truth->shape[1])
        figure.traces.push_back(line("true state",time,sampledColumn(*truth,component,indices),
                                     "truth","#94a3b8","dash"));
    return figure;
}

Figure phaseFigure(const std::string& suffix, const std::vector<std::string>& labels,
                   const NpyArray& states, const std::vector<std::size_t>& indices,
                   const NpyArray* target, const NpyArray* truth) {
    Figure figure{"phase", "Phase space", labels[0] + suffix, labels[1] + suffix, true};
    figure.traces.push_back(line("solution", sampledColumn(states, 0, indices),
                                 sampledColumn(states, 1, indices)));
    if (target && target->shape.size() == 2 && target->shape[1] >= 2)
        figure.traces.push_back(line("target", sampledColumn(*target, 0, indices),
            sampledColumn(*target, 1, indices), "control_target", "#ef4444", "dash"));
    if (truth && truth->shape.size()==2 && truth->shape[1]>=2)
        figure.traces.push_back(line("true state",sampledColumn(*truth,0,indices),
            sampledColumn(*truth,1,indices),"truth","#94a3b8","dash"));
    return figure;
}

Figure threeDimensionalFigure(const std::string& suffix,
                              const std::vector<std::string>& labels,
                              const NpyArray& states,
                              const std::vector<std::size_t>& indices) {
    Figure figure{"three_d", "Three-dimensional state space", labels[0] + suffix,
                  labels[1] + suffix, true};
    figure.zLabel = labels[2] + suffix;
    Trace trace = line("solution", sampledColumn(states, 0, indices),
                       sampledColumn(states, 1, indices));
    trace.z = sampledColumn(states, 2, indices); figure.traces.push_back(std::move(trace));
    return figure;
}

Figure controlFigure(const std::string& suffix, const std::string& timeLabel,
                     const std::vector<double>& time,
                     const NpyArray& controls, const std::vector<std::size_t>& indices) {
    Figure figure{"time", "Control history", timeLabel, "control" + suffix, false};
    const std::size_t columns = controls.shape.size() == 2 ? controls.shape[1] : 1;
    for (std::size_t column = 0; column < columns; ++column) {
        std::vector<double> values;
        if (controls.shape.size() == 2) values = sampledColumn(controls, column, indices);
        else values = sampledVector(&controls, controls.elementCount(), indices);
        figure.traces.push_back(line("u_{" + std::to_string(column) + "}", time,
                                     std::move(values), "control"));
    }
    return figure;
}

Payload trajectoryPayload(const std::filesystem::path& source, const std::string& view,
                          const NpzRecord& record, const json::Value& metadata) {
    const NpyArray* archivedStates = select(record,
        {"values", "states", "controlled_states", "free_training_values"}, 2);
    if (!archivedStates || archivedStates->shape[0] == 0)
        throw std::runtime_error("archive has no state matrix");
    const auto coordinates = displayCoordinateTransform(metadata, archivedStates->shape[1]);
    const auto transformedStates = transformDisplayStates(*archivedStates, coordinates);
    const NpyArray* states = &transformedStates;
    const auto indices = displayIndices(states->shape[0], 100000);
    const NpyArray* timeArray = select(record,
        {"time", "controlled_time", "control_time", "free_time"}, 1);
    const auto time = transformDisplayTime(
        sampledVector(timeArray, states->shape[0], indices), coordinates);
    const NpyArray* archivedTarget = select(record,
        {"control_target_values", "target_values"}, 2);
    const NpyArray* archivedTruth = select(record,{"true_states","truth_values"},2);
    const auto transformedTarget = archivedTarget
        ? transformDisplayStates(*archivedTarget, coordinates) : NpyArray{};
    const auto transformedTruth = archivedTruth
        ? transformDisplayStates(*archivedTruth, coordinates) : NpyArray{};
    const NpyArray* target = archivedTarget ? &transformedTarget : nullptr;
    const NpyArray* truth = archivedTruth ? &transformedTruth : nullptr;
    const NpyArray* controls = select(record,
        {"control_values", "controls", "controlled_controls", "controlled_controls_sample_aligned"});
    const std::string& suffix = coordinates.suffix;
    const auto& labels = coordinates.stateLabels;

    Payload payload{source.u8string(), view, archiveLabel(metadata, source), 2};
    if (!coordinates.warning.empty()) payload.warnings.push_back(coordinates.warning);
    if (view == "phase") {
        if (states->shape[1] < 2) throw std::runtime_error("phase view requires two state components");
        payload.figures.push_back(phaseFigure(suffix, labels, *states, indices, target, truth));
        appendAnalyticOverlay(payload, record, coordinates, time, indices, states->shape[0], states->shape[1]);
        decorateTrajectory(payload, metadata, coordinates, time);
        return payload;
    }
    if (view == "three_d") {
        if (states->shape[1] < 3) throw std::runtime_error("three-state view requires three components");
        payload.figures.push_back(threeDimensionalFigure(suffix, labels, *states, indices));
        appendAnalyticOverlay(payload, record, coordinates, time, indices, states->shape[0], states->shape[1]);
        decorateTrajectory(payload, metadata, coordinates, time);
        return payload;
    }
    if (view == "control" && !controls)
        payload.warnings.push_back("This archive has no saved control history.");
    if (controls) payload.figures.push_back(controlFigure(controlAxisSuffix(metadata),
        coordinates.timeLabel, time, *controls, indices));
    for (std::size_t component = 0; component < states->shape[1]; ++component)
        payload.figures.push_back(componentFigure(labels[component], suffix,
                                                  coordinates.timeLabel,
                                                  time, *states, component, indices, target, truth));
    if ((view == "control" || view == "data_overview") && states->shape[1] >= 2)
        payload.figures.push_back(phaseFigure(suffix, labels, *states, indices, target, truth));
    appendAnalyticOverlay(payload, record, coordinates, time, indices, states->shape[0], states->shape[1]);
    decorateTrajectory(payload, metadata, coordinates, time);
    return payload;
}

Payload spectrumPayload(const std::filesystem::path& source, const std::string& view,
                        const NpzRecord& record, const json::Value& metadata) {
    const NpyArray* values = select(record, {"discrete_eigenvalues"}, 1);
    bool discrete = values != nullptr;
    if (!values) values = select(record, {"continuous_eigenvalues",
                                          "lambda_values", "lambdas"}, 1);
    if (!values) {
        values = select(record, {"eigenvalues"}, 1);
        discrete = metadataString(metadata, "eigenvalue_domain") != "continuous-time";
    }
    if (!values) throw std::runtime_error("archive has no one-dimensional eigenvalue spectrum");
    std::vector<std::complex<double>> spectrum;
    if (values->type == dynsyslab::core::ArrayType::complex128) spectrum = values->complexDoubles();
    else for (const double value : values->doubles()) spectrum.emplace_back(value, 0.0);
    Trace trace; trace.label = "eigenvalues"; trace.marker = true;
    for (const auto value : spectrum) { trace.x.push_back(value.real()); trace.y.push_back(value.imag()); }
    Figure figure{"koopman_spectrum", discrete ? "Discrete-time Koopman spectrum" :
                  "Continuous-time Koopman spectrum", discrete ? "Re(mu)" : "Re(lambda)",
                  discrete ? "Im(mu)" : "Im(lambda)", true};
    figure.traces.push_back(std::move(trace));
    Payload payload{source.u8string(), view, archiveLabel(metadata, source), 1};
    payload.figures.push_back(std::move(figure));
    return payload;
}

Figure errorFigure(const NpzRecord& record, const char* arrayName,
                   const std::string& title) {
    const auto& values = record.at(arrayName);
    if (values.shape.size() != 2) throw std::invalid_argument("control error field must be a matrix");
    Trace trace; trace.label = title; trace.heatmap = true;
    trace.x = record.at("control_error_axis0").doubles();
    trace.y = record.at("control_error_axis1").doubles();
    trace.z = values.doubles();
    trace.zShape = {static_cast<int>(values.shape[0]), static_cast<int>(values.shape[1])};
    Figure figure{"control_map_error", title, "state 0", "state 1", true};
    figure.traces.push_back(std::move(trace)); return figure;
}

Payload controlErrorPayload(const std::filesystem::path& source, const std::string& view,
                            const NpzRecord& record, const json::Value& metadata) {
    if (!record.contains("control_error_axis0") || !record.contains("control_error_axis1"))
        throw std::runtime_error("archive has no native control-map error grid");
    Payload payload{source.u8string(), view, archiveLabel(metadata, source), 2};
    const std::pair<const char*, const char*> fields[]{
        {"constant_absolute_error", "Constant B absolute error"},
        {"constant_relative_error", "Constant B relative error"},
        {"bilinear_absolute_error", "Bilinear B absolute error"},
        {"bilinear_relative_error", "Bilinear B relative error"}};
    for (const auto& field : fields) if (record.contains(field.first))
        payload.figures.push_back(errorFigure(record, field.first, field.second));
    if (payload.figures.empty()) throw std::runtime_error("archive has no control-map error fields");
    return payload;
}

Payload observabilityPayload(const std::filesystem::path& source,const std::string& view,
                             const NpzRecord& record,const json::Value& metadata) {
    const auto& singular=record.at("observability_singular_values");
    const auto& directions=record.at("observable_directions");
    if (singular.shape.size()!=1 || directions.shape.size()!=2)
        throw std::runtime_error("archive has no native observability analysis");
    Trace spectrum; spectrum.label="observability singular values"; spectrum.marker=true;
    const auto values=singular.doubles();
    for (std::size_t index=0;index<values.size();++index) {
        spectrum.x.push_back(static_cast<double>(index));
        spectrum.y.push_back(std::log10(std::max(values[index],1.0e-300)));
    }
    Figure singularFigure{"observability_spectrum","Observability singular spectrum",
                          "direction","log10(singular value)",false};
    singularFigure.traces.push_back(std::move(spectrum));
    Trace map; map.label="right singular directions"; map.heatmap=true;
    map.z=directions.doubles();
    map.zShape={static_cast<int>(directions.shape[0]),static_cast<int>(directions.shape[1])};
    for (std::size_t i=0;i<directions.shape[1];++i) map.x.push_back(static_cast<double>(i));
    for (std::size_t i=0;i<directions.shape[0];++i) map.y.push_back(static_cast<double>(i));
    Figure directionFigure{"observability_directions","Observable lifted-state directions",
                           "singular direction","lifted component",true};
    directionFigure.traces.push_back(std::move(map));
    Payload payload{source.u8string(),view,archiveLabel(metadata,source),2};
    payload.figures.push_back(std::move(singularFigure));
    payload.figures.push_back(std::move(directionFigure));
    return payload;
}

Figure matrixHeatmap(const NpyArray& values,const std::string& id,
                     const std::string& title,const std::string& xLabel,
                     const std::string& yLabel) {
    if (values.shape.size()!=2)
        throw std::invalid_argument("decomposition heat map must be a matrix");
    Trace trace; trace.label=title; trace.heatmap=true; trace.z=values.doubles();
    trace.zShape={static_cast<int>(values.shape[0]),static_cast<int>(values.shape[1])};
    for (std::size_t column=0;column<values.shape[1];++column)
        trace.x.push_back(static_cast<double>(column));
    for (std::size_t row=0;row<values.shape[0];++row)
        trace.y.push_back(static_cast<double>(row));
    Figure figure{id,title,xLabel,yLabel,false};
    figure.traces.push_back(std::move(trace)); return figure;
}

Figure singularFigure(const NpyArray& singular) {
    Trace trace; trace.label="singular value"; trace.marker=true;
    const auto values=singular.doubles();
    for (std::size_t index=0;index<values.size();++index) {
        trace.x.push_back(static_cast<double>(index)); trace.y.push_back(values[index]);
    }
    Figure figure{"singular_values","Singular spectrum","mode","singular value",false};
    figure.traces.push_back(std::move(trace)); return figure;
}

Payload decompositionPayload(const std::filesystem::path& source,
    const std::string& view,const NpzRecord& record,const json::Value& metadata) {
    Payload payload{source.u8string(),view,archiveLabel(metadata,source),2};
    if (record.contains("low_rank") && record.contains("sparse")) {
        const auto& low=record.at("low_rank"); const auto& sparse=record.at("sparse");
        if (low.shape.size()!=2 || sparse.shape!=low.shape)
            throw std::runtime_error("robust PCA arrays have incompatible geometry");
        const auto indices=displayIndices(low.shape[0]);
        const auto time=sampledVector(select(record,{"time"},1),low.shape[0],indices);
        for (std::size_t component=0;component<low.shape[1];++component) {
            Figure lowFigure{"robust_pca_low","Low-rank component "+std::to_string(component),
                             "time","value",false};
            lowFigure.traces.push_back(line("low rank",time,
                sampledColumn(low,component,indices)));
            payload.figures.push_back(std::move(lowFigure));
            Figure sparseFigure{"robust_pca_sparse","Sparse component "+std::to_string(component),
                                "time","value",false};
            sparseFigure.traces.push_back(line("sparse",time,
                sampledColumn(sparse,component,indices),"sparse","#ef4444"));
            payload.figures.push_back(std::move(sparseFigure));
        }
        return payload;
    }
    if (record.contains("singular_values"))
        payload.figures.push_back(singularFigure(record.at("singular_values")));
    if (record.contains("explained_variance_ratio")) {
        Trace trace; trace.label="explained variance"; trace.marker=true;
        const auto ratios=record.at("explained_variance_ratio").doubles();
        for (std::size_t index=0;index<ratios.size();++index) {
            trace.x.push_back(static_cast<double>(index)); trace.y.push_back(ratios[index]);
        }
        Figure figure{"explained_variance","Explained variance ratio","mode","fraction",false};
        figure.traces.push_back(std::move(trace)); payload.figures.push_back(std::move(figure));
    }
    if (record.contains("components"))
        payload.figures.push_back(matrixHeatmap(record.at("components"),
            "pca_components","Principal components","feature","mode"));
    if (record.contains("scores")) {
        const auto& scores=record.at("scores");
        if (scores.shape.size()!=2) throw std::runtime_error("PCA scores must be a matrix");
        const auto indices=displayIndices(scores.shape[0]);
        const auto time=sampledVector(select(record,{"time"},1),scores.shape[0],indices);
        Figure figure{"pca_scores","Principal-component scores","time","score",false};
        for (std::size_t mode=0;mode<scores.shape[1];++mode)
            figure.traces.push_back(line("mode "+std::to_string(mode),time,
                sampledColumn(scores,mode,indices)));
        payload.figures.push_back(std::move(figure));
    } else if (record.contains("left_singular_vectors")) {
        const auto& left=record.at("left_singular_vectors");
        if (left.shape.size()==2) {
            const auto indices=displayIndices(left.shape[0]);
            const auto time=sampledVector(select(record,{"time"},1),left.shape[0],indices);
            Figure figure{"left_singular_vectors","Temporal singular vectors",
                          "time","amplitude",false};
            for (std::size_t mode=0;mode<left.shape[1];++mode)
                figure.traces.push_back(line("mode "+std::to_string(mode),time,
                    sampledColumn(left,mode,indices)));
            payload.figures.push_back(std::move(figure));
        }
    }
    const NpyArray* right=select(record,{"right_singular_vectors"},2);
    if (right && !record.contains("components"))
        payload.figures.push_back(matrixHeatmap(*right,"right_singular_vectors",
            "Spatial / delay singular vectors","feature","mode"));
    if (payload.figures.empty())
        throw std::runtime_error("archive has no native decomposition payload");
    return payload;
}

}  // namespace

Payload buildNativePayload(const std::filesystem::path& source, const std::string& view,
                           const NpzRecord& record, const json::Value& metadata) {
    if(view=="sindy_dynamics"||view=="sindy_dynamics_transient")
        return buildSindyDynamicsPayload(source,view,record,metadata);
    if (view == "control_map_error" || view == "control_error")
        return controlErrorPayload(source, view, record, metadata);
    if (view == "krom_dashboard")
        return buildKromDashboardPayload(source, record, record.metadataJson());
    if (view == "observable_time" || view == "observable_heatmap" || view == "observable_workspace")
        return buildObservablePayload(source, view, record, metadata);
    if (view == "lattice") return buildLatticePayload(source, view, record, metadata);
    if (view == "koopman_validation")
        return buildKoopmanValidationPayload(source,view,record,metadata);
    if (view == "eigenvalue_optimization")
        return buildEigenvalueOptimizationPayload(source,view,record,metadata);
    if(view=="krom_control_campaign")
        return buildKromControlCampaignPayload(source,view,record,metadata);
    if(view=="fourier") return buildFourierPayload(source,record,metadata);
    if(view=="control_fourier")
        return buildControlFourierPayload(source,record,metadata);
    if(view=="time_frequency") {
        auto payload = trajectoryPayload(source, "time", record, metadata);
        auto fourier = buildFourierPayload(source, record, metadata);
        payload.view = view;
        payload.warnings.insert(payload.warnings.end(),
            fourier.warnings.begin(), fourier.warnings.end());
        payload.figures.insert(payload.figures.end(),
            std::make_move_iterator(fourier.figures.begin()),
            std::make_move_iterator(fourier.figures.end()));
        return payload;
    }
    if (view=="observability")
        return observabilityPayload(source,view,record,metadata);
    if (record.contains("low_rank") || record.contains("singular_values"))
        return decompositionPayload(source,view,record,metadata);
    if (view == "koopman_spectrum" || view == "krom")
        return spectrumPayload(source, view, record, metadata);
    if (view == "field" || metadataString(metadata, "system_family") == "pde")
        return buildFieldPayload(source, view, record, metadata);
    return trajectoryPayload(source, view, record, metadata);
}

Payload buildNativePayload(const std::filesystem::path& source, const std::string& view) {
    const NpzRecord record = dynsyslab::core::readNpz(source);
    json::Value metadata;
    if (!record.metadataJson().empty()) metadata = json::parse(record.metadataJson());
    return buildNativePayload(source, view, record, metadata);
}

}  // namespace dsl
