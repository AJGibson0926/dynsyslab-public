#include "SpectralMeasurePayload.hpp"

#include "ComplexStemPlot.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <string>

namespace dsl {
namespace {

using dynsyslab::core::ArrayType;
using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;

const NpyArray* complexArray(const NpzRecord& record, const char* name) {
    const auto found = record.arrays.find(name);
    return found != record.arrays.end() && found->second.type == ArrayType::complex128
        ? &found->second : nullptr;
}

const NpyArray* complexMatrix(const NpzRecord& record, const char* name, std::size_t columns) {
    const auto* array = complexArray(record, name);
    return array && array->shape.size() == 2 && array->shape[0] > 0
        && array->shape[1] == columns ? array : nullptr;
}

std::string channelName(const dynsyslab::core::json::Value& metadata, std::size_t channel) {
    const auto& names = metadata.at("observables");
    if (!names.isNull() && channel < names.array().size()) {
        const auto name = names.array()[channel].string();
        if (!name.empty()) return name;
    }
    return "observable " + std::to_string(channel);
}

}  // namespace

std::vector<Figure> spectralMeasureFigures(const NpzRecord& record,
                                           const dynsyslab::core::json::Value& metadata) {
    const auto* discrete = complexArray(record, "discrete_eigenvalues");
    if (!discrete) return {};
    const auto mu = discrete->complexDoubles();
    const auto rank = mu.size();
    const auto* phi = complexMatrix(record, "eigenfunction_values", rank);
    const auto* map = complexMatrix(record, "reconstruction_matrix", rank);
    if (!phi || !map || rank == 0) return {};
    std::vector<std::complex<double>> lambda;
    if (const auto* continuous = complexArray(record, "continuous_eigenvalues"))
        lambda = continuous->complexDoubles();
    const auto phiValues = phi->complexDoubles();  // row 0: first aligned sample
    const auto mapValues = map->complexDoubles();
    std::vector<Figure> figures;
    for (std::size_t channel = 0; channel < map->shape[0]; ++channel) {
        ComplexStemPlot plot;
        plot.heightLabel = "|c_j phi_j(t0)|";
        plot.seriesLabel = "modal weight";
        plot.heights.resize(rank);
        for (std::size_t mode = 0; mode < rank; ++mode)
            plot.heights[mode] = std::abs(mapValues[channel * rank + mode] * phiValues[mode]);
        const auto name = channelName(metadata, channel);
        plot.title = "Spectral measure of " + name + " on mu";
        plot.baseSymbol = "mu"; plot.bases = mu; plot.unitCircle = true;
        figures.push_back(complexPlaneStemFigure(plot));
        if (lambda.size() == rank) {
            plot.title = "Spectral measure of " + name + " on lambda";
            plot.baseSymbol = "lambda"; plot.bases = lambda; plot.unitCircle = false;
            figures.push_back(complexPlaneStemFigure(plot));
            // Zoom near the real axis: |Im lambda| up to 12x the frequency of the
            // heaviest oscillatory mode (all modes remain in the full panel above).
            double dominant = 0.0, heaviest = -1.0;
            for (std::size_t mode = 0; mode < rank; ++mode)
                if (std::abs(lambda[mode].imag()) > 1e-9 && plot.heights[mode] > heaviest) {
                    heaviest = plot.heights[mode]; dominant = std::abs(lambda[mode].imag());
                }
            if (dominant > 0) {
                ComplexStemPlot zoom = plot; zoom.bases.clear(); zoom.heights.clear();
                const double limit = 12.0 * dominant;
                for (std::size_t mode = 0; mode < rank; ++mode)
                    if (std::abs(lambda[mode].imag()) <= limit) {
                        zoom.bases.push_back(lambda[mode]); zoom.heights.push_back(plot.heights[mode]);
                    }
                char text[96]; std::snprintf(text, sizeof text, " (|Im lambda| <= %.3g, %zu modes)",
                                             limit, zoom.bases.size());
                zoom.title = "Spectral measure of " + name + " near the real axis" + text;
                figures.push_back(complexPlaneStemFigure(zoom));
            }
        }
    }
    return figures;
}

}  // namespace dsl
