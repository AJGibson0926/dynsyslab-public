#include <optional>
#include <dynsyslab/core/SpectrumArchive.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include "ExpansionArchive.hpp"
#include "PostprocessInternal.hpp"
#include "TimeSeriesSource.hpp"

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
        source.values(), source.samples(), source.dimension(), uniformDeltaTime(source), options);
    dynsyslab::core::NpzRecord archive;
    dynsyslab::core::writeSpectrumArrays(archive, result.spectrum);
    archive.arrays["modes"] = dynsyslab::core::NpyArray::fromComplex(
        {source.dimension(), result.rank}, result.modes);
    archive.arrays["eigenfunction_values"] = dynsyslab::core::NpyArray::fromComplex(
        {source.samples(), result.rank}, result.coordinates);
    archive.arrays["operator"] = dynsyslab::core::NpyArray::fromComplex(
        {result.rank, result.rank}, result.operatorMatrix);
    const auto meta = metadata(source, "dmd", label,
            "\"eigenvalue_domain\":\"discrete-time\",\"delta_t\":"
            + std::to_string(result.deltaTime) + ",\"rank\":"
            + std::to_string(result.rank));
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(meta);
    {   // psi = Phi^+ x: the eigenfunctions as a linear expansion in the state coordinates.
        const std::size_t n = source.dimension(), r = result.rank;
        std::vector<std::complex<double>> identity(n * n);
        for (std::size_t i = 0; i < n; ++i) identity[i * n + i] = 1.0;
        numerics::LeastSquaresOptions solve;
        solve.method = numerics::LeastSquaresMethod::divideAndConquerSvd;
        const auto inverse = numerics::solveLeastSquares(result.modes, n, r, identity, n, solve);  // r x n
        dynsyslab::core::PolynomialDictionary linear;
        linear.coordinateDimension = n;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) linear.exponents.push_back(i == j ? 1 : 0);
        std::optional<dynsyslab::core::SystemConfiguration> system;
        try { system = systemConfiguration(source); } catch (const std::exception&) {}
        std::string reason;
        try {
            dynsyslab::core::EigenfunctionFitResult bank{dynsyslab::core::FunctionExpansion(std::move(linear), inverse, snapshotDomain(source.metadataJson, n, system)), result.spectrum, {}};
            // DMD (Schmid, J. Fluid Mech. 656, 2010; Tu et al., J. Comput. Dyn. 1, 2014) publishes its
            // eigenfunctions as a data object: psi = Phi^+ x as an expansion; modes and operator are properties.
            const dynsyslab::core::DataOrigin origin{"dmd", settings, {dataReference(source.path, "time_series")}};
            if (writeExpansionKrom(archive, bank, origin, source.time(), result.coordinates,
                                   source.record.channels, meta, reason))
                return publishObject(root, source, dynsyslab::core::object_kind::krom, "koopman-dmd", std::move(archive));
        } catch (const std::exception& error) {
            reason = error.what();
        }
        // No K-ROM (psi = Phi^+ x could not be formed): the spectrum, modes and operator as a report.
        archive.arrays["expansion_contract_skipped"] = dynsyslab::core::NpyArray::fromString(reason);
    }
    return publishObject(root, source, dynsyslab::core::object_kind::report, "koopman-dmd", std::move(archive));
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
