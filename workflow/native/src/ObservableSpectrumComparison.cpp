// Post-processing method "observable_spectrum_comparison": H-DMD spectra of different scalar
// observables of the same states, compared pairwise.
//
// Sources are H-DMD spectrum archives fitted with identical settings and different "observable"
// (typically one hdmd stage per observable). {"names":[...] (default each archive's observable),
//  "require_identical_settings":true}. Settings that must agree: matrix_size, sampling_window,
// rank, sample_stride, delta_t, coordinate system. The archive holds the spectra (discrete and
// continuous, concatenated, with spectrum_offsets) and the symmetric nearest-neighbour RMS distance
// matrices in both domains (core/SpectrumComparison.hpp). Hankel DMD: Arbabi & Mezic, SIAM J.
// Appl. Dyn. Syst. 16(4), 2096 (2017).
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SpectrumComparison.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
constexpr const char* kSettings[] = {"matrix_size", "sampling_window", "rank", "sample_stride", "delta_t", "coordinate_system"};
}  // namespace

std::filesystem::path runObservableSpectrumComparison(const std::vector<std::filesystem::path>& sources,
                                                      const json::Value& settings, const std::filesystem::path& root,
                                                      const std::string& label) {
    if (sources.size() < 2) throw std::invalid_argument("observable_spectrum_comparison needs at least two spectra");
    std::vector<std::vector<std::complex<double>>> discrete, continuous;
    json::Value::Array names, sizes, differences;
    json::Value reference;
    for (std::size_t k = 0; k < sources.size(); ++k) {
        const auto record = c::readNpz(sources[k]);
        const auto metadata = record.metadataJson().empty() ? json::Value{} : json::parse(record.metadataJson());
        if (!record.contains("continuous_eigenvalues")) throw std::invalid_argument("source " + std::to_string(k) + " has no spectrum");
        const auto rates = record.at("continuous_eigenvalues").complexDoubles();
        std::vector<std::complex<double>> multipliers;
        if (record.contains("discrete_eigenvalues")) multipliers = record.at("discrete_eigenvalues").complexDoubles();
        else if (record.contains("sample_interval"))
            for (const auto& rate : rates) multipliers.push_back(std::exp(rate * record.at("sample_interval").doubles().at(0)));
        continuous.push_back(rates); discrete.push_back(std::move(multipliers));
        const auto& given = settings.at("names").array();
        names.emplace_back(k < given.size() ? given[k].string() : metadata.at("observable").string("observable " + std::to_string(k)));
        sizes.emplace_back(static_cast<double>(rates.size()));
        if (k == 0) { reference = metadata; continue; }
        for (const char* key : kSettings)
            if (json::dump(metadata.at(key)) != json::dump(reference.at(key)))
                differences.emplace_back(json::Value::Object{{"setting", json::Value(key)}, {"source", json::Value(static_cast<double>(k))},
                    {"value", metadata.at(key)}, {"reference", reference.at(key)}});
    }
    if (!differences.empty() && settings.at("require_identical_settings").boolean(true))
        throw std::invalid_argument("observable_spectrum_comparison: the spectra were fitted with different settings: "
                                    + json::dump(json::Value(differences)));
    for (std::size_t a = 0; a < names.size(); ++a)
        for (std::size_t b = a + 1; b < names.size(); ++b)
            if (names[a].string() == names[b].string()) throw std::invalid_argument("observable names must be unique");
    const std::size_t n = sources.size();
    std::vector<double> discreteRms(n * n, 0.0), continuousRms(n * n, 0.0);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b) {
            discreteRms[a * n + b] = discreteRms[b * n + a] = c::symmetricNearestRms(discrete[a], discrete[b]);
            continuousRms[a * n + b] = continuousRms[b * n + a] = c::symmetricNearestRms(continuous[a], continuous[b]);
        }
    std::vector<std::complex<double>> allDiscrete, allContinuous;
    std::vector<std::int64_t> offsets{0};
    for (std::size_t k = 0; k < n; ++k) {
        allContinuous.insert(allContinuous.end(), continuous[k].begin(), continuous[k].end());
        auto multipliers = discrete[k];
        multipliers.resize(continuous[k].size(), std::complex<double>(std::nan(""), std::nan("")));
        allDiscrete.insert(allDiscrete.end(), multipliers.begin(), multipliers.end());
        offsets.push_back(static_cast<std::int64_t>(allContinuous.size()));
    }
    c::NpzRecord out;
    out.arrays["continuous_eigenvalues"] = c::NpyArray::fromComplex({allContinuous.size()}, allContinuous);
    out.arrays["discrete_eigenvalues"] = c::NpyArray::fromComplex({allDiscrete.size()}, allDiscrete);
    out.arrays["spectrum_offsets"] = c::NpyArray::fromInt64s({offsets.size()}, offsets);
    out.arrays["pairwise_discrete_nearest_rms"] = c::NpyArray::fromDoubles({n, n}, discreteRms);
    out.arrays["pairwise_continuous_nearest_rms"] = c::NpyArray::fromDoubles({n, n}, continuousRms);
    json::Value::Array sourceList;
    for (const auto& source : sources) sourceList.emplace_back(source.generic_u8string());
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(json::Value::Object{
        {"schema", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("observable-spectrum-comparison")},
        {"encoding", json::Value("hdmd-observable-spectrum-comparison/v1")}, {"label", json::Value(label)},
        {"observable_names", json::Value(names)}, {"spectrum_sizes", json::Value(sizes)},
        {"source_archives", json::Value(sourceList)}, {"settings", reference.at("settings")},
        {"setting_differences", json::Value(differences)}, {"same_state_training_data", json::Value(true)}})));
    const auto directory = kindFolder(root, c::object_kind::figure);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    stampDataObject(out, c::object_kind::figure, "observable-spectrum-comparison", {});
    c::writeNpzAtomicNew(path, out);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
