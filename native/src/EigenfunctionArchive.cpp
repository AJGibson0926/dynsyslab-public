#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

json::Value metadataOf(const NpzRecord& archive) {
    const auto text = archive.metadataJson();
    return text.empty() ? json::Value{} : json::parse(text);
}

}  // namespace

json::Value writeEigenfunctionSet(NpzRecord& archive, const EigenfunctionSet& set,
                                  const std::vector<ChannelDescriptor>& inputChannels) {
    const auto& samples = set.samples();
    const std::size_t n = samples.samples(), m = samples.modes;
    archive.arrays["continuous_eigenvalues"] = NpyArray::fromComplex({m}, set.spectrum().continuousRates());
    if (set.spectrum().sampleInterval() > 0) {  // a sampled spectrum; hand-given rates have no interval
        archive.arrays["discrete_eigenvalues"] = NpyArray::fromComplex({m}, set.spectrum().discreteMultipliers());
        archive.arrays["eigenvalues"] = archive.arrays["discrete_eigenvalues"];
        archive.arrays["sample_interval"] = NpyArray::fromDoubles({1}, {set.spectrum().sampleInterval()});
    }
    archive.arrays["eigenfunction_values"] = NpyArray::fromComplex({n, m}, samples.values);
    archive.arrays["training_time"] = NpyArray::fromDoubles({n}, samples.time);
    std::vector<double> indices(samples.sampleIndices.begin(), samples.sampleIndices.end());
    archive.arrays["training_sample_indices"] = NpyArray::fromDoubles({n}, indices);
    json::Value::Array representations;
    for (const auto& rep : set.representations()) {
        auto description = rep->describe();
        if (const auto* shared = dynamic_cast<const SharedRepresentation*>(rep.get())) {
            if (const auto* expansion = dynamic_cast<const FunctionExpansion*>(shared->function().get())) {
                // One dictionary; row j of the coefficients is eigenfunction j's own expansion.
                description.object()["dictionary"] = expansion->dictionary().describe();
                description.object()["coefficients"] = "expansion_coefficients";
                archive.arrays["expansion_coefficients"] = NpyArray::fromComplex(
                    {m, expansion->dictionary().size()}, expansion->coefficients());
                representations.emplace_back(std::move(description));
                continue;
            }
            const auto* interpolation = dynamic_cast<const FunctionInterpolation*>(shared->function().get());
            if (!interpolation) throw std::invalid_argument("this shared representation cannot be archived");
            // The base points belong to the interpolation; the values are eigenfunction_values.
            archive.arrays["interpolation_points"] = NpyArray::fromDoubles({n, shared->function()->coordinateDimension()},
                                                                          interpolation->points());
        } else if (const auto* window = dynamic_cast<const WindowMapRepresentation*>(rep.get())) {
            // Row j is eigenfunction j's own map from a raw window, so each one stands alone.
            std::vector<std::size_t> all(m);
            for (std::size_t j = 0; j < m; ++j) all[j] = j;
            archive.arrays["window_map"] = NpyArray::fromComplex({m, window->layout().rows()},
                window->materialized() ? window->map() : window->mapRows(all));
        }
        representations.emplace_back(std::move(description));
    }
    return json::Value::Object{{"schema", kEigenfunctionObjectsSchema}, {"origin", dataOriginJson(set.origin())},
                               {"input_channels", channelDescriptorsJson(inputChannels)},
                               {"representations", representations},
                               {"default_representation", set.representations().front()->kind()}};
}

bool hasEigenfunctionObjects(const NpzRecord& archive) {
    return metadataOf(archive).at("eigenfunction_objects").at("schema").string() == kEigenfunctionObjectsSchema;
}

EigenfunctionSet readEigenfunctionSet(const NpzRecord& archive, const std::vector<std::string>& kinds) {
    const auto block = metadataOf(archive).at("eigenfunction_objects");
    if (block.at("schema").string() != kEigenfunctionObjectsSchema)
        throw std::invalid_argument("archive has no eigenfunction objects");
    const auto rates = archive.at("continuous_eigenvalues").complexDoubles();
    const bool sampled = archive.contains("discrete_eigenvalues");
    const auto multipliers = sampled ? archive.at("discrete_eigenvalues").complexDoubles() : std::vector<std::complex<double>>{};
    const auto time = archive.at("training_time").doubles();
    const auto origin = dataOriginFromJson(block.at("origin"));
    auto samples = std::make_shared<EigenfunctionSamples>();
    samples->time = time;
    for (const double index : archive.at("training_sample_indices").doubles()) samples->sampleIndices.push_back(std::size_t(index));
    samples->values = archive.at("eigenfunction_values").complexDoubles();
    samples->modes = rates.size();
    const double dt = time.size() > 1 ? time[1] - time[0] : 0.0;
    const double interval = archive.contains("sample_interval") ? archive.at("sample_interval").doubles().front() : dt;
    std::vector<std::shared_ptr<const EigenfunctionRepresentation>> reps;
    for (const auto& description : block.at("representations").array()) {
        const auto kind = description.at("kind").string();
        const auto canonical = kind == "window_map" ? std::string("expansion") : kind;
        if (!kinds.empty() && std::find(kinds.begin(), kinds.end(), canonical) == kinds.end()) continue;
        if (kind == "interpolation") {
            if (!description.at("domains").isNull())
                throw std::invalid_argument("this eigenfunction archive predates shared interpolation "
                                            "(per-domain points); regenerate it");
            const auto& points = archive.at("interpolation_points");
            InterpolationData data;
            data.samples = samples->samples();
            data.stateDimension = points.shape.at(1);
            data.trainingStates = points.doubles();
            data.trainingValues = samples->values;
            data.method = interpolationMethodFromJson(description.at("method"), &data.neighbors);
            data.provenance = "eigenfunction interpolation";
            reps.push_back(std::make_shared<SharedRepresentation>(std::make_shared<const FunctionInterpolation>(
                std::move(data), functionDomainFromJson(description.at("domain")))));
        } else if (kind == "expansion" && !description.at("dictionary").isNull()) {
            auto dictionary = std::make_shared<const ExpansionDictionary>(ExpansionDictionary::fromJson(description.at("dictionary")));
            reps.push_back(std::make_shared<SharedRepresentation>(std::make_shared<const FunctionExpansion>(
                std::move(dictionary), archive.at("expansion_coefficients").complexDoubles(),
                functionDomainFromJson(description.at("domain")))));
        } else if (kind == "window_map" || (kind == "expansion" && description.at("family").string() == "linear-window")) {
            if (!archive.contains("window_map"))
                throw std::invalid_argument("this eigenfunction archive predates stand-alone window maps; regenerate it");
            reps.push_back(WindowMapRepresentation::fromMap(hankelWindowLayoutFromJson(description),
                                                            archive.at("window_map").complexDoubles(), rates.size()));
        } else {
            throw std::invalid_argument("unknown eigenfunction representation " + kind);
        }
    }
    auto spectrum = sampled ? Spectrum::fromDiscrete(multipliers, rates, interval, origin) : Spectrum::fromContinuousRates(rates);
    return EigenfunctionSet(std::move(spectrum), *samples, origin, std::move(reps));
}

Eigenfunction readEigenfunction(const NpzRecord& archive, std::size_t mode, const std::string& representation,
                               const ObservableCatalog* catalog) {
    const std::vector<std::string> kinds = representation.empty() ? std::vector<std::string>{} : std::vector<std::string>{representation};
    const auto set = readEigenfunctionSet(archive, kinds).extract({mode});
    auto eigenfunction = set.eigenfunction(0, {}, representation, catalog);
    return eigenfunction;
}

}  // namespace dynsyslab::core
