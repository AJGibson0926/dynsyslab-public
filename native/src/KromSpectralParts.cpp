// The spectral parts of a krom archive: the compressed operator of a spectral fit and the functions on
// its spectral sets (KoopmanModel::compression, KoopmanModel::spectralFunction). Bulk arrays are NPZ
// arrays; the krom block holds their descriptions.
#include "KromArchiveParts.hpp"

#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/FunctionDomain.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <cctype>

#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
std::string key(const std::string& value) {
    std::string out;
    for (const char c : value) out += (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') ? c : '_';
    return out;
}
}  // namespace

json::Value writeKromCompression(NpzRecord& archive, const KoopmanCompression& c) {
    c.validate();
    const auto n = c.dictionary.size();
    archive.arrays["compression_basis"] = NpyArray::fromDoubles({n, c.rank}, c.basis);
    archive.arrays["compression_action"] = NpyArray::fromDoubles({c.rank, c.rank}, c.action);
    json::Value::Object readouts;
    for (const auto& [name, readout] : c.readouts) {
        const auto array = "compression_readout__" + key(name);
        archive.arrays[array] = NpyArray::fromDoubles({n, readout.outputs}, readout.coefficients);
        readouts[name] = json::Value::Object{{"array", array}, {"outputs", double(readout.outputs)},
                                             {"relative_training_error", readout.relativeTrainingError}};
    }
    return json::Value::Object{{"schema", "dynsyslab.koopman-compression/v1"}, {"dictionary", spectralDictionaryJson(c.dictionary)},
        {"method", c.method}, {"assumptions", c.assumptions}, {"coordinate_contract", c.coordinateContract},
        {"sample_interval", c.sampleInterval}, {"rank", double(c.rank)}, {"samples", double(c.samples)},
        {"basis", "compression_basis"}, {"action", "compression_action"}, {"readouts", json::Value(std::move(readouts))}};
}

KoopmanCompression readKromCompression(const json::Value& block, const NpzRecord& archive) {
    if (block.at("schema").string() != "dynsyslab.koopman-compression/v1") throw std::invalid_argument("unknown krom compression schema");
    KoopmanCompression c;
    c.dictionary = spectralDictionaryFromJson(block.at("dictionary"));
    c.method = block.at("method").string(); c.assumptions = block.at("assumptions").string();
    c.coordinateContract = block.at("coordinate_contract").string("measured");
    c.sampleInterval = block.at("sample_interval").number(); c.rank = std::size_t(block.at("rank").number());
    c.samples = std::size_t(block.at("samples").number());
    c.basis = archive.at(block.at("basis").string()).doubles();
    c.action = archive.at(block.at("action").string()).doubles();
    for (const auto& [name, readout] : block.at("readouts").object())
        c.readouts[name] = FeatureReadout{std::size_t(readout.at("outputs").number()), archive.at(readout.at("array").string()).doubles(),
                                          readout.at("relative_training_error").number()};
    c.validate();
    return c;
}

json::Value::Array writeKromSpectralFunctions(NpzRecord& archive, const KoopmanModel& model, bool ownOnly) {
    json::Value::Array out;
    for (const auto& name : model.spectralFunctionKeys()) {
        if (ownOnly && !model.ownsSpectralFunction(name)) continue;
        const auto& function = *model.spectralFunction(name);
        const auto array = "spectral_function__" + key(name);
        archive.arrays[array] = NpyArray::fromComplex({function.outputCount(), function.dictionary().size()}, function.coefficients());
        out.emplace_back(json::Value::Object{{"key", name}, {"dictionary", function.dictionary().describe()},
            {"domain", functionDomainJson(function.coordinates())}, {"coefficients", array}});
    }
    return out;
}

void readKromSpectralFunctions(const json::Value& descriptions, const NpzRecord& archive, KoopmanModel& model) {
    for (const auto& description : descriptions.array()) {
        auto dictionary = std::make_shared<const ExpansionDictionary>(ExpansionDictionary::fromJson(description.at("dictionary")));
        model.setSpectralFunction(description.at("key").string(), FunctionExpansion(std::move(dictionary),
            archive.at(description.at("coefficients").string()).complexDoubles(), functionDomainFromJson(description.at("domain"))));
    }
}

}  // namespace dynsyslab::core::detail
