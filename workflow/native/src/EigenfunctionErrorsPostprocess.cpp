// Postprocess method eigenfunction_errors: error metrics of an eigenfunction set, swept over
// interpolant methods (and optionally its expansion), to choose the interpolant for an application.
//
//   {"name": "eigenfunction_errors", "settings": {
//      "interpolation_methods": [{"name": "inverse_distance", "neighbors": 8, "power": 2},
//                                {"name": "rbf_thin_plate", "neighbors": 16}],   default: every method
//      "domain": "state",                         default: the interpolation's first domain
//      "include_expansion": true}}                also score the window expansion (needs the input)
//
// Output arrays (rows in the order evaluated): value_error and koopman_residual, rows x modes;
// metadata "rows" names each row.
#include "PostprocessInternal.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/EigenfunctionArchive.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::filesystem::path runEigenfunctionErrors(const std::filesystem::path& sourcePath, const J& settings,
                                             const std::filesystem::path& root, const std::string& label) {
    const auto archive = c::readNpz(sourcePath);
    if (!c::hasEigenfunctionObjects(archive)) throw std::invalid_argument("eigenfunction_errors needs an eigenfunction-objects archive");
    auto set = c::readEigenfunctionSet(archive);
    std::vector<J> methods;
    if (settings.at("interpolation_methods").isNull())
        for (const auto& name : c::interpolationMethods()) methods.push_back(J::Object{{"name", name}});
    else
        for (const auto& item : settings.at("interpolation_methods").array()) methods.push_back(item);
    const auto domain = settings.at("domain").string();
    const std::size_t modes = set.size();
    std::vector<double> value, residual;
    J::Array rows;
    auto record = [&](const c::EigenfunctionErrors& e, J row) {
        value.insert(value.end(), e.valueRelative.begin(), e.valueRelative.end());
        residual.insert(residual.end(), e.koopmanResidual.begin(), e.koopmanResidual.end());
        row.object()["domain"] = e.domain;
        row.object()["samples"] = double(e.samples);
        rows.emplace_back(std::move(row));
    };
    std::size_t defaultNeighbors = c::interpolationOf(set).data().neighbors;
    for (const auto& spec : methods) {
        std::size_t neighbors = defaultNeighbors;
        const auto method = c::interpolationMethodFromJson(spec, &neighbors);
        c::setInterpolationMethod(set, method, neighbors);
        record(set.errors("interpolation", domain),
               J::Object{{"representation", "interpolation"}, {"method", c::interpolationMethodJson(method, neighbors)}});
    }
    if (settings.at("include_expansion").boolean(false)) {
        std::vector<std::filesystem::path> inputs;
        for (const auto& input : set.origin().inputs)
            if (input.role == "time_series") inputs.push_back(std::filesystem::u8path(input.path));
        const auto input = loadTimeSeriesRecords(inputs);
        record(set.errors("expansion", {}, &input), J::Object{{"representation", "expansion"}});
    }
    c::NpzRecord out;
    out.arrays["value_error"] = c::NpyArray::fromDoubles({rows.size(), modes}, value);
    out.arrays["koopman_residual"] = c::NpyArray::fromDoubles({rows.size(), modes}, residual);
    out.arrays["continuous_eigenvalues"] = c::NpyArray::fromComplex({modes}, set.spectrum().continuousRates());
    SourceData source;  // provenance only: the eigenfunction archive is not a time series
    source.path = sourcePath;
    source.metadataJson = archive.metadataJson();
    auto meta = json::parse(metadata(source, "eigenfunction-errors", label, "\"source_eigenfunctions\":\""
        + escaped(std::filesystem::absolute(sourcePath).generic_u8string()) + "\""));
    meta.object()["rows"] = rows;
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
    return publish(root, source, "eigenfunction-errors", label, std::move(out));
}

}  // namespace dynsyslab::workflow::postprocess_detail
