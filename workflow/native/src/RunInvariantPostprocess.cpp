// Post-processing method run_invariant: a conserved polynomial of the coordinates learned from run
// identity (core/RunInvariant.hpp). Pools every source of the job; each is one unforced run.
//
//   {"name": "run_invariant", "settings": {"degree": 3, "pair_stride": 1}}
//
// Publishes a report object: the invariant as a polynomial model (invariant_model_json, the
// dictionary-and-coefficients form; coefficients over polynomialExponentTable(d, degree, true)),
// the separation rho, and per run its mean, spread and sampled invariance defect, with its values
// (invariant_<k>).
#include "GeometryPostprocess.hpp"
#include "KromDatasetSupport.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/RunInvariant.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::filesystem::path runRunInvariant(const std::vector<std::filesystem::path>& sources, const J& settings,
                                      const std::filesystem::path& root, const std::string& label) {
    std::vector<c::TimeSeries> runs;
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        if (datasetDynamics(source) == "forced") throw std::invalid_argument("run_invariant needs unforced runs; " + path.u8string() + " is forced");
        runs.push_back(source.record.series);
    }
    const auto degree = std::size_t(settings.at("degree").number(3));
    const auto fit = c::fitRunInvariant(runs, degree, std::size_t(settings.at("pair_stride").number(1)));
    const auto d = fit.dimension;
    const c::PolynomialDictionary dictionary{d, c::polynomialExponentTable(d, degree, true), std::vector<double>(d, 0), std::vector<double>(d, 1)};
    const c::SindyDynamicsModel model(dictionary.asObservable(), fit.coefficients, 1);
    c::NpzRecord archive;
    archive.arrays["invariant_model_json"] = c::NpyArray::fromString(json::dump(model.definition()));
    archive.arrays["coefficients"] = c::NpyArray::fromDoubles({fit.coefficients.size()}, fit.coefficients);
    J::Array perRun;
    for (std::size_t k = 0; k < runs.size(); ++k) {
        std::vector<double> values;
        for (std::size_t s = 0; s < runs[k].sampleCount(); ++s) values.push_back(fit.evaluate(runs[k].valuesData() + s * d));
        archive.arrays["invariant_" + std::to_string(k)] = c::NpyArray::fromDoubles({values.size()}, values);
        const auto& r = fit.runs[k];
        perRun.push_back(J::Object{{"source", sources[k].generic_u8string()}, {"mean", r.mean}, {"standard_deviation", r.standardDeviation},
                                   {"peak_to_peak", r.peakToPeak}, {"sampled_invariance_relative_l2", r.sampledInvariance}});
    }
    const auto first = loadSource(sources.front());
    auto meta = json::parse(metadata(first, "run-invariant", label, "")).object();
    meta["settings"] = settings;
    meta["separation"] = fit.separation;
    meta["runs"] = perRun;
    meta["eigenvalue"] = "continuous 0, multiplier 1 (a conserved quantity)";
    meta["references"] = J::Array{"Fisher 1936, Ann. Eugenics 7:179-188", "Mezic 2005, Nonlinear Dyn. 41:309-325"};
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, first, c::object_kind::report, "run-invariant", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
