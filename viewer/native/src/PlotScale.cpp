#include "dynsyslab_viewer/PlotScale.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dsl {
double ValueScale::forward(double value) const {
    if (kind == ValueScaleKind::linear) return value;
    if (kind == ValueScaleKind::logarithmic)
        return value > 0 ? std::log10(value) : std::numeric_limits<double>::quiet_NaN();
    const double a = std::abs(value) / linearThreshold;
    return std::copysign(a <= 1 ? a : 1 + decadeWeight * std::log10(a),value);
}
double ValueScale::inverse(double value) const {
    if (kind == ValueScaleKind::linear) return value;
    if (kind == ValueScaleKind::logarithmic) return std::pow(10.,std::clamp(value,-307.,307.));
    const double a = std::abs(value);
    return std::copysign(linearThreshold*(a <= 1 ? a : std::pow(10.,std::min((a-1)/decadeWeight,307.))),value);
}
std::string ValueScale::name() const {
    return kind == ValueScaleKind::linear ? "linear" : kind == ValueScaleKind::logarithmic ? "log" : "symlog";
}
bool scalesColorValues(const Figure& figure) {
    return figure.kind == "observable_heatmap" || figure.kind == "krom_heatmap" || figure.kind == "domain_field";
}
bool supportsValueScale(const Figure& figure) {
    if (scalesColorValues(figure))
        return !figure.traces.empty() && figure.traces.front().encoding != "domain-color";
    return figure.kind == "time" || figure.kind == "observable_time" || figure.kind == "krom_time"
        || figure.kind == "stem" || (figure.kind == "three_d" && figure.boxAxes3D)
        || figure.kind == "control" || figure.kind == "klqe"
        // Linear-systems spectra (Hankel singular values, gramians, frequency response).
        || figure.kind == "hankel_spectrum" || figure.kind == "gramian_spectrum"
        || figure.kind == "frequency_response"
        // Fourier amplitude spectra (archive routes and time-series --spectrum).
        || figure.kind == "fourier" || figure.kind == "control_fourier";
}
ValueScale resolveValueScale(const Figure& figure) {
    if (!supportsValueScale(figure)) return {ValueScaleKind::linear,1,"Coordinates or complex hue remain linear"};
    if (figure.valueScale == "linear") return {ValueScaleKind::linear,1,"Manual linear scale"};
    std::vector<double> magnitudes;
    bool nonpositive = false;
    double maximum = 0;
    for (const auto& trace : figure.traces) {
        const auto& values = scalesColorValues(figure) ? trace.z : trace.y;
        const auto stride = std::max<std::size_t>(1,values.size()/4096);
        for (std::size_t i = 0; i < values.size(); ++i) {
            const auto value = values[i];
            if (!std::isfinite(value)) continue;
            nonpositive = nonpositive || value <= 0;
            maximum = std::max(maximum,std::abs(value));
            if (i % stride == 0 && value != 0) magnitudes.push_back(std::abs(value));
        }
    }
    if (magnitudes.empty()) return {ValueScaleKind::linear,1,"No nonzero finite values"};
    std::sort(magnitudes.begin(),magnitudes.end());
    const auto low = magnitudes[(magnitudes.size()-1)*5/100];
    const auto high = magnitudes[(magnitudes.size()-1)*95/100];
    const auto threshold = std::max(low,maximum*1e-6);
    if (figure.valueScale == "symlog") return {ValueScaleKind::symmetric_log,threshold,"Manual symmetric log scale"};
    if (figure.valueScale == "log")
        return nonpositive ? ValueScale{ValueScaleKind::symmetric_log,threshold,"Log uses symmetric log to retain zero and negative values"}
            : ValueScale{ValueScaleKind::logarithmic,threshold,"Manual log scale"};
    if (std::log10(high)-std::log10(low) < 3)
        return {ValueScaleKind::linear,threshold,"Central 90% spans fewer than three orders of magnitude"};
    return {nonpositive ? ValueScaleKind::symmetric_log : ValueScaleKind::logarithmic,threshold,
        nonpositive ? "Wide dynamic range with zero or negative values" : "Positive values span at least three orders of magnitude"};
}
ValueScale ordinateScale(const Figure& figure) {
    return scalesColorValues(figure) ? ValueScale{} : resolveValueScale(figure);
}
void setValueScale(Payload& payload, const std::string& requested) {
    if (requested != "auto" && requested != "linear" && requested != "log" && requested != "symlog")
        throw std::invalid_argument("Value scale must be auto, linear, log, or symlog");
    for (auto& figure : payload.figures) figure.valueScale = requested;
}
}
