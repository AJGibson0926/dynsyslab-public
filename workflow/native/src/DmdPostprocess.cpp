#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Extension.hpp>
#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/core/Text.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runDmd(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label) {
    dynsyslab::core::DMDOptions options;
    options.rank = static_cast<std::size_t>(settings.at("rank").number());
    options.singularValueTolerance = settings.at("svd_tolerance").number();
    options.exactModes = settings.at("exact").boolean(true);
    const auto result = dynsyslab::core::dynamicModeDecomposition(
        source.values, source.samples, source.dimension, uniformDeltaTime(source), options);
    dynsyslab::core::NpzRecord archive;
    archive.arrays["eigenvalues"] = dynsyslab::core::NpyArray::fromComplex(
        {result.rank}, result.discreteEigenvalues);
    archive.arrays["discrete_eigenvalues"] = archive.arrays["eigenvalues"];
    archive.arrays["continuous_eigenvalues"] = dynsyslab::core::NpyArray::fromComplex(
        {result.rank}, result.continuousEigenvalues);
    archive.arrays["modes"] = dynsyslab::core::NpyArray::fromComplex(
        {source.dimension, result.rank}, result.modes);
    archive.arrays["eigenfunction_values"] = dynsyslab::core::NpyArray::fromComplex(
        {source.samples, result.rank}, result.coordinates);
    archive.arrays["operator"] = dynsyslab::core::NpyArray::fromComplex(
        {result.rank, result.rank}, result.operatorMatrix);
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        metadata(source, "dmd", label,
            "\"eigenvalue_domain\":\"discrete-time\",\"delta_t\":"
            + std::to_string(result.deltaTime) + ",\"rank\":"
            + std::to_string(result.rank)));
    return publish(root, source, "koopman-dmd", label, std::move(archive));
}

std::filesystem::path runExtensionPostprocess(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label,
    const std::string& method, const std::string& extension) {
    json::Value::Object request;
    request["schema"] = "dynsyslab.extension-postprocess/v1";
    request["source_archive"] = source.path.generic_u8string();
    request["archive_root"] = root.generic_u8string();
    request["label"] = label;
    request["method"] = method;
    request["settings"] = settings;
    const auto response = json::parse(
        dynsyslab::core::ExtensionRegistry::instance().execute(
            extension, "postprocessor:" + method,
            json::dump(json::Value(std::move(request)))));
    auto output = std::filesystem::path(
        dynsyslab::core::wide(response.at("archive").string()));
    if (output.empty())
        throw std::runtime_error(
            "extension postprocessor did not return an archive path");
    if (output.is_relative()) output = root / output;
    if (!std::filesystem::is_regular_file(output))
        throw std::runtime_error("extension postprocessor archive does not exist");
    return output;
}

}  // namespace dynsyslab::workflow::postprocess_detail
