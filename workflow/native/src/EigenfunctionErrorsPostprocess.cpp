// Postprocess method eigenfunction_errors: error metrics of an eigenfunction set, swept over
// interpolant methods (and optionally its expansion), to choose the interpolant for an application.
//
//   {"name": "eigenfunction_errors", "settings": {
//      "interpolation_methods": [{"name": "inverse_distance", "neighbors": 8, "power": 2},
//                                {"name": "rbf_thin_plate", "neighbors": 16}],   default: every method
//      "domain": "state",                         optional check of the interpolation's domain
//      "include_expansion": true,                 also score the time-delay window map (needs the input)
//      "spectral_weight": true}}                  also write each mode's empirical spectral weight
//
// Output arrays (rows in the order evaluated): value_error and koopman_residual, rows x modes;
// metadata "rows" names each row.
#include "PostprocessInternal.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::filesystem::path runEigenfunctionErrors(const std::filesystem::path& sourcePath, const J& settings,
                                             const std::filesystem::path& root, const std::string& label) {
    const auto archive = c::readNpz(sourcePath);
    if (!c::isKrom(archive)) throw std::invalid_argument("eigenfunction_errors needs a krom");
    auto set = c::loadKoopmanModel(sourcePath);
    std::vector<J> methods;
    if (settings.at("interpolation_methods").isNull())
        for (const auto& name : c::interpolationMethods()) methods.push_back(J::Object{{"name", name}});
    else
        for (const auto& item : settings.at("interpolation_methods").array()) methods.push_back(item);
    // The interpolation has one domain, fixed when the eigenfunctions were made; "domain" only checks it.
    const auto domain = settings.at("domain").string();
    const auto own = c::domainName(set.interpolation().coordinates());
    if (!domain.empty() && domain != own)
        throw std::invalid_argument("these eigenfunctions are interpolated on " + own + ", not " + domain);
    const std::size_t modes = set.liftedDimension();
    std::vector<double> value, residual;
    J::Array rows;
    auto record = [&](const c::EigenfunctionErrors& e, J row) {
        value.insert(value.end(), e.valueRelative.begin(), e.valueRelative.end());
        residual.insert(residual.end(), e.koopmanResidual.begin(), e.koopmanResidual.end());
        row.object()["domain"] = e.domain;
        row.object()["samples"] = double(e.samples);
        rows.emplace_back(std::move(row));
    };
    std::size_t defaultNeighbors = set.interpolation().data().neighbors;
    for (const auto& spec : methods) {
        std::size_t neighbors = defaultNeighbors;
        const auto method = c::interpolationMethodFromJson(spec, &neighbors);
        set.setInterpolationMethod(method, neighbors);
        record(set.errors("interpolation"),
               J::Object{{"representation", "interpolation"}, {"method", c::interpolationMethodJson(method, neighbors)}});
    }
    if (settings.at("include_expansion").boolean(false)) {
        std::vector<std::filesystem::path> inputs;
        for (const auto& input : set.origin().inputs)
            if (input.role == "time_series") inputs.push_back(std::filesystem::u8path(input.path));
        const auto input = loadTimeSeriesRecords(inputs);
        record(set.errors("time-delay", &input), J::Object{{"representation", "time-delay"}});
    }
    c::NpzRecord out;
    out.arrays["value_error"] = c::NpyArray::fromDoubles({rows.size(), modes}, value);
    out.arrays["koopman_residual"] = c::NpyArray::fromDoubles({rows.size(), modes}, residual);
    out.arrays["continuous_eigenvalues"] = c::NpyArray::fromComplex({modes}, set.spectrum().continuousRates());
    if (settings.at("spectral_weight").boolean(false))  // samples x modes (KoopmanModel.hpp empiricalSpectralWeight)
        out.arrays["spectral_weight"] = c::NpyArray::fromDoubles({set.samples().samples(), modes}, c::empiricalSpectralWeight(set.samples()));
    SourceData source;  // provenance only: the eigenfunction archive is not a time series
    source.path = sourcePath;
    source.metadataJson = archive.metadataJson();
    auto meta = json::parse(metadata(source, "eigenfunction-errors", label, "\"source_eigenfunctions\":\""
        + escaped(std::filesystem::absolute(sourcePath).generic_u8string()) + "\""));
    meta.object()["rows"] = rows;
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
    return publishObject(root, source, dynsyslab::core::object_kind::report, "eigenfunction-errors", std::move(out));
}

}  // namespace dynsyslab::workflow::postprocess_detail
