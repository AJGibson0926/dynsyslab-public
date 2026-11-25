#include "ExpansionArchive.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/EigenfunctionArchive.hpp>

#include <stdexcept>

namespace dynsyslab::workflow {
namespace c = dynsyslab::core;

c::FunctionDomain snapshotDomain(const std::string& sourceMetadataJson, std::size_t dimension,
                                      const std::optional<c::SystemConfiguration>& system) {
    namespace json = dynsyslab::core::json;
    const auto meta = sourceMetadataJson.empty() ? json::Value{} : json::parse(sourceMetadataJson);
    const auto method = meta.at("method").string();
    const bool solve = method.empty() || method == "solve";
    std::size_t stateDimension = 0;
    if (system) {
        try { stateDimension = c::createSystem(*system)->descriptor().stateDimension; } catch (const std::exception&) {}
    }
    if (solve && stateDimension == dimension) return c::FunctionDomain::physicalState(dimension);
    c::FunctionDomain data;
    data.dimension = dimension;
    return data;
}

bool attachExpansionContract(c::NpzRecord& archive, const c::EigenfunctionFitResult& bank,
                             const std::optional<c::SystemConfiguration>& system, const std::string& metadataJson,
                             std::string& reason) {
    try {
        const auto model = system ? c::KoopmanModel(bank, *system) : c::KoopmanModel(bank, std::size_t(0));
        auto arrays = archive.arrays;
        arrays.erase("metadata_json");
        archive = c::makeKoopmanArchive(model, metadataJson, std::move(arrays));
        return true;
    } catch (const std::exception& error) {
        reason = error.what();
        return false;
    }
}

bool writeExpansionEigenfunctions(c::NpzRecord& archive, const c::EigenfunctionFitResult& bank, const c::DataOrigin& origin,
                                  const std::vector<double>& time, const std::vector<std::complex<double>>& values,
                                  const std::vector<c::ChannelDescriptor>& channels, const std::string& metadataJson,
                                  std::string& reason) {
    namespace json = dynsyslab::core::json;
    try {
        c::EigenfunctionSamples samples;
        samples.modes = bank.spectrum.size();
        if (!values.empty()) {
            samples.values = values;
            samples.time = time;
            for (std::size_t k = 0; k < time.size(); ++k) samples.sampleIndices.push_back(k);
        }
        const c::EigenfunctionSet set(bank.spectrum, samples, origin,
            {std::make_shared<c::SharedRepresentation>(std::make_shared<const c::FunctionExpansion>(bank.functions))});
        auto objects = c::writeEigenfunctionSet(archive, set, channels);
        auto meta = metadataJson.empty() ? json::Value(json::Value::Object{}) : json::parse(metadataJson);
        meta.object()["eigenfunction_objects"] = std::move(objects);
        archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
        return true;
    } catch (const std::exception& error) {
        reason = error.what();
        return false;
    }
}

}  // namespace dynsyslab::workflow
