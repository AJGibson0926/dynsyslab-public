#include <dynsyslab/core/KoopmanArchive.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <stdexcept>

namespace dynsyslab::core {

ApproximateControlMap loadApproximateControlMap(
    const std::filesystem::path& path, const std::string& prefix,
    std::size_t modes, std::size_t inputs) {
    const auto archive = readNpz(path);
    const std::string base = prefix.empty() ? std::string{} : prefix + "_";
    const auto provenance = "loaded independently from " + path.u8string()
        + (prefix.empty() ? std::string{} : " [" + prefix + "]");
    const std::vector<std::string> matrixKeys{base + "matrix",
        base + "constant_control_matrix", "constant_control_matrix",
        "control_constant_matrix"};
    for (const auto& key : matrixKeys) if (archive.contains(key)) {
        const auto& value = archive.at(key);
        if (value.shape != std::vector<std::size_t>{modes,inputs})
            throw std::invalid_argument("constant control map shape does not match KROM");
        ConstantControlMap result{value.complexDoubles(), modes, inputs, provenance};
        const std::string keyPrefix = key.substr(
            0, key.size() - std::string("matrix").size());
        if (archive.contains(keyPrefix + "discrete_matrix"))
            result.discreteMatrix = archive.at(
                keyPrefix + "discrete_matrix").complexDoubles();
        if (archive.contains(keyPrefix + "sample_time"))
            result.sampleTime = archive.at(keyPrefix + "sample_time").doubles().front();
        return result;
    }
    const std::vector<std::string> operatorKeys{base + "operators",
        base + "bilinear_control_operators", "bilinear_control_operators",
        "control_bilinear_operators"};
    for (const auto& key : operatorKeys) if (archive.contains(key)) {
        const auto suffix = key.substr(0, key.size() - std::string("operators").size());
        const std::vector<std::string> offsetKeys{suffix + "offsets", base + "offsets",
            "bilinear_control_offsets", "control_bilinear_offsets"};
        const NpyArray* offsets = nullptr;
        for (const auto& candidate : offsetKeys)
            if (archive.contains(candidate)) { offsets = &archive.at(candidate); break; }
        if (!offsets)
            throw std::invalid_argument("bilinear control map offsets are missing");
        const auto& operators = archive.at(key);
        if (operators.shape != std::vector<std::size_t>{inputs,modes,modes}
            || offsets->shape != std::vector<std::size_t>{modes,inputs})
            throw std::invalid_argument("bilinear control map shape does not match KROM");
        BilinearControlMap result;
        result.operators = operators.complexDoubles();
        result.offsets = offsets->complexDoubles();
        result.modes = modes; result.inputs = inputs; result.provenance = provenance;
        const std::vector<std::string> conjugateKeys{suffix + "conjugate_operators",
            base + "conjugate_operators", "bilinear_control_conjugate_operators",
            "control_bilinear_conjugate_operators"};
        for (const auto& candidate : conjugateKeys)
            if (archive.contains(candidate)) {
                result.conjugateOperators = archive.at(candidate).complexDoubles();
                break;
            }
        return result;
    }
    throw std::invalid_argument(
        "control-map archive contains no supported representation");
}

ApproximateControlMap loadApproximateControlMap(
    const std::filesystem::path& path,const std::string& prefix,
    const KoopmanModel& model,bool externalMap) {
    const auto archive=readNpz(path);
    const auto key=prefix.empty()?"model_identity":prefix+"_model_identity";
    std::string identity=archive.contains(key)?archive.at(key).string():std::string{};
    if(identity.empty()&&archive.contains("input_map_manifest_json")) {
        const auto manifest=json::parse(archive.at("input_map_manifest_json").string());
        for(const auto& entry:manifest.array())
            if(prefix=="control_"+entry.at("name").string())identity=entry.at("model_identity").string();
    }
    if(identity.empty()&&!externalMap)
        throw std::invalid_argument("External control map records no KROM identity; refit it on this KROM, "
                                    "or set control_map_origin \"external\" to attach it on your word");
    if(!identity.empty()&&identity!=koopmanModelIdentity(model))
        throw std::invalid_argument("External control map does not match ordered eigenfunctions, normalization, physics or time convention");
    auto result=loadApproximateControlMap(path,prefix,model.eigenfunctionCount(),model.inputDimension());
    std::visit([&](auto& map) {
        map.modelIdentity=identity.empty()?koopmanModelIdentity(model):identity;
        if(identity.empty())map.provenance+="; external map attached without a recorded K-ROM identity";
    },result);
    return result;
}

}  // namespace dynsyslab::core
