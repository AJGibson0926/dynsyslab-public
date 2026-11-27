#include "ExpansionArchive.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include "ObjectPublishing.hpp"

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

bool writeExpansionKrom(c::NpzRecord& archive, const c::EigenfunctionFitResult& bank, const c::DataOrigin& origin,
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
        const c::KoopmanModel model(bank.spectrum, samples, origin, {std::make_shared<const c::FunctionExpansion>(bank.functions)});
        auto extra = archive.arrays;  // the method's own arrays; the krom writes its spectrum and samples
        for (const char* key : {"metadata_json", "eigenfunction_values", "training_time", "training_sample_indices"}) extra.erase(key);
        const auto stamped = postprocess_detail::withDataObject(metadataJson, c::object_kind::krom, {});
        archive = c::makeKoopmanArchive(model, stamped, std::move(extra), channels);
        return true;
    } catch (const std::exception& error) {
        reason = error.what();
        return false;
    }
}

}  // namespace dynsyslab::workflow
