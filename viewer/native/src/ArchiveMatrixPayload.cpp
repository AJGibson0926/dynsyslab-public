#include "ArchiveMatrixPayload.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace dsl {
namespace {
using dynsyslab::core::NpzRecord;
namespace json = dynsyslab::core::json;

std::string archiveLabel(const json::Value& metadata, const std::filesystem::path& source) {
    return metadata.at("label").string(source.stem().u8string());
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

}  // namespace

Payload buildControlErrorPayload(const std::filesystem::path& source, const std::string& view,
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

Payload buildObservabilityPayload(const std::filesystem::path& source,const std::string& view,
                             const NpzRecord& record,const json::Value& metadata) {
    const auto& singular=record.at("observability_singular_values");
    const auto& directions=record.contains("observable_directions")?record.at("observable_directions")
        :record.at("observability_directions");
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
    if (record.contains("mode_participation")) {  // groups x directions, written by the pipeline
        const auto& participation=record.at("mode_participation");
        Trace heat; heat.label="fractional squared participation"; heat.heatmap=true;
        heat.z=participation.doubles();
        heat.zShape={static_cast<int>(participation.shape[0]),static_cast<int>(participation.shape[1])};
        for (std::size_t i=0;i<participation.shape[1];++i) heat.x.push_back(static_cast<double>(i));
        for (std::size_t i=0;i<participation.shape[0];++i) heat.y.push_back(static_cast<double>(i));
        Figure figure{"observability_participation","Koopman-mode participation by observability direction",
                      "singular direction","nonredundant Koopman mode",true};
        figure.colorLabel="fractional squared participation"; figure.integerXTicks=true;
        figure.traces.push_back(std::move(heat));
        payload.figures.push_back(std::move(figure));
    }
    if (record.contains("observer_poles")) {
        Trace poles; poles.label="observer poles"; poles.marker=true; poles.color="#984EA3";
        for (const auto& pole : record.at("observer_poles").complexDoubles()) {
            poles.x.push_back(pole.real()); poles.y.push_back(pole.imag());
        }
        Figure figure{"observer_poles","Estimator observer poles","real part","imaginary part",true};
        figure.traces.push_back(std::move(poles));
        payload.figures.push_back(std::move(figure));
    }
    if (record.contains("state_error_norm") && record.contains("time")) {
        Trace norm; norm.label="||x - x_hat||_2"; norm.color="#D62728";
        norm.x=record.at("time").doubles(); norm.y=record.at("state_error_norm").doubles();
        Figure figure{"estimation_error","State-estimation error norm","time","Euclidean state error",false};
        figure.traces.push_back(std::move(norm));
        payload.figures.push_back(std::move(figure));
    }
    return payload;
}

Payload buildBifurcationPayload(const std::filesystem::path& source, const std::string& view,
                                const NpzRecord& record, const json::Value& metadata) {
    if (!record.contains("parameter") || !record.contains("maximum"))
        throw std::runtime_error("archive has no bifurcation points");
    Trace points; points.label = "post-transient maxima"; points.marker = true; points.color = "#1F4E79";
    points.x = record.at("parameter").doubles(); points.y = record.at("maximum").doubles();
    Figure figure{"bifurcation", "Bifurcation diagram (peak amplitudes)",
                  metadata.at("parameter_label").string("parameter"), metadata.at("value_label").string("maximum"), false};
    figure.traces.push_back(std::move(points));
    Payload payload{source.u8string(), view, archiveLabel(metadata, source), 1};
    payload.figures.push_back(std::move(figure));
    return payload;
}

Payload buildObservableSpectraPayload(const std::filesystem::path& source, const std::string& view,
                                     const NpzRecord& record, const json::Value& metadata) {
    const auto discrete = record.at("discrete_eigenvalues").complexDoubles();
    const auto continuous = record.at("continuous_eigenvalues").complexDoubles();
    std::vector<std::int64_t> offsets;
    for (const double value : record.at("spectrum_offsets").doubles()) offsets.push_back(static_cast<std::int64_t>(value));
    constexpr const char* colors[] = {"#3478C8", "#E45756", "#54A24B", "#F2CF5B", "#B279A2", "#FF9DA6", "#9D755D"};
    Figure unit{"koopman_spectrum", "H-DMD spectra by scalar observable (discrete time)", "Re(mu)", "Im(mu)", true};
    Figure rates{"koopman_spectrum", "H-DMD spectra by scalar observable (continuous time)", "Re(lambda)", "Im(lambda)", false};
    Trace circle; circle.label = "unit circle"; circle.color = "#777777";
    for (int k = 0; k <= 360; ++k) { circle.x.push_back(std::cos(k * 3.14159265358979323846 / 180)); circle.y.push_back(std::sin(k * 3.14159265358979323846 / 180)); }
    unit.traces.push_back(std::move(circle));
    const auto& names = metadata.at("observable_names").array();
    for (std::size_t k = 0; k + 1 < offsets.size(); ++k) {
        Trace a, b; a.marker = b.marker = true; a.color = b.color = colors[k % 7];
        a.label = b.label = k < names.size() ? names[k].string() : "observable " + std::to_string(k);
        for (auto i = offsets[k]; i < offsets[k + 1]; ++i) {
            a.x.push_back(discrete[i].real()); a.y.push_back(discrete[i].imag());
            b.x.push_back(continuous[i].real()); b.y.push_back(continuous[i].imag());
        }
        unit.traces.push_back(std::move(a)); rates.traces.push_back(std::move(b));
    }
    const auto& distances = record.at("pairwise_discrete_nearest_rms");
    Trace heat; heat.heatmap = true; heat.label = "nearest-spectrum RMS (discrete)";
    heat.z = distances.doubles(); heat.zShape = {static_cast<int>(distances.shape[0]), static_cast<int>(distances.shape[1])};
    for (std::size_t i = 0; i < distances.shape[0]; ++i) { heat.x.push_back(double(i)); heat.y.push_back(double(i)); }
    Figure matrix{"observable_spectrum_distance", "Symmetric nearest-spectrum RMS distance (discrete)", "observable", "observable", true};
    matrix.colorLabel = "RMS distance"; matrix.integerXTicks = true; matrix.traces.push_back(std::move(heat));
    Payload payload{source.u8string(), view, archiveLabel(metadata, source), 2};
    payload.figures = {std::move(unit), std::move(rates), std::move(matrix)};
    return payload;
}

}  // namespace dsl
