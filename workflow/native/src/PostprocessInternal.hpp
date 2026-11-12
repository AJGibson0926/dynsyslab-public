#pragma once

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/SINDy.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

namespace json = dynsyslab::core::json;

struct SourceData {
    std::filesystem::path path;
    std::vector<double> time;
    std::vector<double> values;
    std::size_t samples{};
    std::size_t dimension{};
    std::vector<double> inputs;
    std::size_t inputDimension{};
    std::string metadataJson;
};

std::string escaped(const std::string& value);
std::string timestamp();
std::string safeName(std::string value);
std::string safeArrayKey(std::string value);
void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned,
                    const std::string& message);
SourceData loadSource(const std::filesystem::path& path);
void normalizeSourceCoordinates(SourceData& source);
double uniformDeltaTime(const SourceData& source);
std::string metadata(const SourceData& source, const std::string& method,
                     const std::string& label, const std::string& details);
std::filesystem::path publish(const std::filesystem::path& root,
    const SourceData& source, const std::string& method, const std::string& label,
    dynsyslab::core::NpzRecord archive);
std::string jsonArray(const std::vector<double>& values);
std::string jsonStringArray(const std::vector<std::string>& values);
dynsyslab::core::SystemConfiguration systemConfiguration(const SourceData& source);
std::vector<double> exactVelocities(const SourceData& source, bool includeInputTerms = true);

std::filesystem::path runKoopmanSpectral(const SourceData&,const json::Value&,
    const std::filesystem::path&,const std::string&,const std::string&);
std::filesystem::path runLinearAnalysis(const std::filesystem::path& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runCoordinateTransform(const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runSindyDynamics(const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runDynamicsValidation(const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runObservableReconstruction(const std::filesystem::path& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runHdmd(const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runDmd(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runDmdControl(const std::filesystem::path& source,
    const json::Value& settings, const std::filesystem::path& root,
    const std::string& label, bool extended = false);
std::filesystem::path runExtensionPostprocess(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label,
    const std::string& method, const std::string& extension);

std::filesystem::path runPcaPostprocess(
    const std::string& method, const SourceData& source,
    const json::Value& settings, const std::filesystem::path& root,
    const std::string& label);

std::filesystem::path runControlRegression(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runKromExtraction(
    const std::filesystem::path& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runKromValidation(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runBoundaryEigenfunctions(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& precedingSpectrum,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runEigenvalueOptimization(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& precedingSpectrum,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runSindy(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& spectrumPath,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runCoordinateSindy(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& spectrumPath,const std::vector<std::complex<double>>& rates,
    const dynsyslab::core::KoopmanSINDyOptions& options,
    const std::filesystem::path& root,const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
