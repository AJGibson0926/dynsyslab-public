// Post-processing method "bifurcation_points": a peak-amplitude bifurcation diagram from runs.
//
// Each run contributes the post-transient local maxima of one state component
// (core::postTransientMaxima), placed at that run's parameter value; a run with no peak in its tail
// (settled on an equilibrium) contributes its final value. Runs are the job's sources,
// or listed explicitly as a run manifest:
//   {"runs":[{"archive":PATH,"parameter":VALUE}, ...],      (explicit parameter per run), or
//    "parameter":"system_parameters.driver.amplitude"       (dotted metadata path; default
//                                                             "campaign_point_index"),
//    "parameter_scale":1.0, "value_scale":1.0, "array":"values", "component":0,
//    "tail_fraction":0.25}
// The archive holds parameter, maximum and run_index arrays (one entry per point) and the run list.
// References: Strogatz (2015) sec. 10.2; Parlitz & Lauterborn, Phys. Rev. A 36, 1428 (1987).
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Extrema.hpp>
#include <dynsyslab/core/Text.hpp>

#include <optional>
#include <sstream>
#include <variant>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;

double metadataNumber(const json::Value& metadata, const std::string& dotted) {
    const json::Value* value = &metadata;
    std::stringstream path(dotted);
    std::string key;
    while (std::getline(path, key, '.')) value = &value->at(key);
    if (!std::holds_alternative<double>(value->storage()))
        throw std::invalid_argument("bifurcation_points: run metadata has no number at '" + dotted + "'");
    return value->number();
}

}  // namespace

std::filesystem::path runBifurcationPoints(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                           const std::filesystem::path& root, const std::string& label) {
    struct Run { std::filesystem::path archive; std::optional<double> parameter; };
    std::vector<Run> runs;
    for (const auto& run : settings.at("runs").array())
        runs.push_back({std::filesystem::u8path(run.at("archive").string()),
                        run.at("parameter").isNull() ? std::nullopt : std::optional<double>(run.at("parameter").number())});
    if (runs.empty()) for (const auto& source : sources) runs.push_back({source, std::nullopt});
    if (runs.empty()) throw std::invalid_argument("bifurcation_points needs sources or runs");
    const auto parameterPath = settings.at("parameter").string("campaign_point_index");
    const auto arrayName = settings.at("array").string("values");
    const auto component = static_cast<std::size_t>(settings.at("component").number(0.0));
    const double tail = settings.at("tail_fraction").number(0.25);
    const double parameterScale = settings.at("parameter_scale").number(1.0);
    const double valueScale = settings.at("value_scale").number(1.0);
    std::vector<double> parameters, maxima, runIndex;
    json::Value::Array runList;
    for (std::size_t r = 0; r < runs.size(); ++r) {
        const auto record = c::readNpz(runs[r].archive);
        const auto metadata = record.metadataJson().empty() ? json::Value{} : json::parse(record.metadataJson());
        const double parameter = parameterScale
            * (runs[r].parameter ? *runs[r].parameter : metadataNumber(metadata, parameterPath));
        const auto& array = record.at(arrayName);
        const std::size_t columns = array.shape.size() == 2 ? array.shape[1] : 1;
        if (component >= columns) throw std::invalid_argument("bifurcation_points: component is outside the state");
        const auto values = array.doubles();
        std::vector<double> series(values.size() / columns);
        for (std::size_t row = 0; row < series.size(); ++row) series[row] = values[row * columns + component];
        auto peaks = c::postTransientMaxima(series, tail);
        if (peaks.empty() && !series.empty()) peaks.push_back(series.size() - 1);  // settled: its final value
        for (const auto row : peaks) {
            parameters.push_back(parameter);
            maxima.push_back(valueScale * series[row]);
            runIndex.push_back(static_cast<double>(r));
        }
        runList.emplace_back(json::Value::Object{{"archive", json::Value(runs[r].archive.generic_u8string())},
            {"parameter", json::Value(parameter)}, {"points", json::Value(static_cast<double>(peaks.size()))}});
    }
    c::NpzRecord out;
    out.arrays["parameter"] = c::NpyArray::fromDoubles({parameters.size()}, parameters);
    out.arrays["maximum"] = c::NpyArray::fromDoubles({maxima.size()}, maxima);
    out.arrays["run_index"] = c::NpyArray::fromDoubles({runIndex.size()}, runIndex);
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(json::Value::Object{
        {"schema", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("bifurcation-points")},
        {"label", json::Value(label)}, {"parameter", json::Value(parameterPath)},
        {"parameter_label", json::Value(settings.at("parameter_label").string(parameterPath))},
        {"value_label", json::Value(settings.at("value_label").string(arrayName + "[" + std::to_string(component) + "]"))},
        {"component", json::Value(static_cast<double>(component))}, {"tail_fraction", json::Value(tail)},
        {"runs", json::Value(std::move(runList))}, {"point_count", json::Value(static_cast<double>(parameters.size()))}})));
    const auto directory = kindFolder(root, c::object_kind::figure);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    stampDataObject(out, c::object_kind::figure, "bifurcation-points", {});
    c::writeNpzAtomicNew(path, out);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
