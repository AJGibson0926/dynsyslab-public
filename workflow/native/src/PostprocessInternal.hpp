#pragma once

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>
#include "ObjectPublishing.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

namespace json = dynsyslab::core::json;

// An opened archive: where it is, its metadata, and its samples as a time-series record (series,
// inputs, channels, origin). The accessors read the record's owned buffers.
struct SourceData {
    std::filesystem::path path;
    std::string metadataJson;
    dynsyslab::core::TimeSeriesRecord record;

    [[nodiscard]] const std::vector<double>& time() const { return record.series.timeVector(); }
    [[nodiscard]] const std::vector<double>& values() const { return record.series.valueVector(); }
    [[nodiscard]] std::size_t samples() const noexcept { return record.series.sampleCount(); }
    [[nodiscard]] std::size_t dimension() const noexcept { return record.series.componentCount(); }
    // Sample x input, or interval x input when record.inputsOnIntervals; empty without inputs.
    [[nodiscard]] const std::vector<double>& inputs() const {
        static const std::vector<double> none;
        return record.inputs ? record.inputs->valueVector() : none;
    }
    [[nodiscard]] std::size_t inputDimension() const noexcept {
        return record.inputs ? record.inputs->componentCount() : 0;
    }
};

// A copy of `like` (path, metadata, inputs, origin; channels when the width is unchanged) holding
// the given samples instead.
SourceData withSeries(const SourceData& like, std::vector<double> time, std::vector<double> values,
                      std::size_t dimension);
// The input-history series for a record with `samples` samples on `time`: one row per sample, or
// per interval (times of interval starts) when `onIntervals`.
dynsyslab::core::TimeSeries inputSeries(const std::vector<double>& time, std::vector<double> inputs,
                                        std::size_t inputDimension, bool onIntervals);

std::string escaped(const std::string& value);
std::string timestamp();
std::string safeName(std::string value);
std::string safeArrayKey(std::string value);
void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned,
                    const std::string& message);
// A record in the integration (dimensionless) coordinates; with `archiveCoordinates`, as stored
// (dimensional values kept, metadata "analysis_coordinates": "archive"; a step's setting
// "analysis_coordinates": "archive" asks for it).
SourceData loadSource(const std::filesystem::path& path, bool archiveCoordinates = false);
void normalizeSourceCoordinates(SourceData& source);
// True when the record is stored in dimensional units and loadSource rescales it.
bool storedInDimensionalUnits(const json::Value& metadata);
double uniformDeltaTime(const SourceData& source);
std::string metadata(const SourceData& source, const std::string& method,
                     const std::string& label, const std::string& details);
// Publishes `model` as a krom (core/KoopmanArchive.hpp) under <root>/objects/krom, named after its
// source and `method`, stamped with the step's provenance; a model loaded from a krom is published as
// its child (only what the step added). `metadataJson` and `extraArrays` hold the method's properties.
std::filesystem::path publishKrom(const std::filesystem::path& root, const SourceData& source,
                                  const dynsyslab::core::KoopmanModel& model, const std::string& method,
                                  const std::string& metadataJson,
                                  std::map<std::string, dynsyslab::core::NpyArray> extraArrays = {},
                                  const std::vector<dynsyslab::core::ChannelDescriptor>& inputChannels = {});
// Writes `archive` as a data object of `kind` (core::object_kind) under <root>/objects/<kind>, named
// after its source and `method`, with the data_object envelope.
std::filesystem::path publishObject(const std::filesystem::path& root, const SourceData& source,
                                    const std::string& kind, const std::string& method,
                                    dynsyslab::core::NpzRecord archive);
std::string jsonArray(const std::vector<double>& values);
std::string jsonStringArray(const std::vector<std::string>& values);
dynsyslab::core::SystemConfiguration systemConfiguration(const SourceData& source);
std::vector<double> exactVelocities(const SourceData& source, bool includeInputTerms = true);

std::filesystem::path runKoopmanSpectral(const SourceData&,const json::Value&,
    const std::filesystem::path&,const std::string&,const std::string&);
std::filesystem::path runLinearAnalysis(const std::filesystem::path& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runLinearSystems(const std::filesystem::path& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runEmpiricalObservability(const std::filesystem::path& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
// The input a coordinate change reads: a measured window integral reads its measurement array,
// every other change the source trajectory (CoordinateTransformPostprocess.cpp).
SourceData loadCoordinateInput(const std::filesystem::path& path, const dynsyslab::core::json::Value& settings);
std::filesystem::path runCoordinateTransform(const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runSindyDynamics(const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runDynamicsValidation(const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
std::filesystem::path runObservableReconstruction(const std::filesystem::path& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label);
// H-DMD pooled over the source and settings "additional_sources" (HdmdPooled.cpp).
std::filesystem::path runPooledHdmd(const SourceData& source, const json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);
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

std::filesystem::path runEigenfunctionErrors(const std::filesystem::path& sourcePath,
    const json::Value& settings, const std::filesystem::path& root, const std::string& label);
std::filesystem::path runEigenfunctionExpansion(const std::filesystem::path& sourcePath,
    const json::Value& settings, const std::filesystem::path& root, const std::string& label);
std::filesystem::path runDomainField(const std::filesystem::path& sourcePath,
    const json::Value& settings, const std::filesystem::path& root, const std::string& label);
std::filesystem::path runHankelSvd(const std::string& method, const SourceData& source,
    const json::Value& settings, const std::filesystem::path& root, const std::string& label);
std::filesystem::path runPcaPostprocess(
    const std::string& method, const SourceData& source,
    const json::Value& settings, const std::filesystem::path& root,
    const std::string& label);

std::filesystem::path runControlRegression(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runHistoryControlRegression(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
// B_psi for any KROM: exact / history_exact (analytic) or constant / bilinear regressed from forced records.
std::filesystem::path runInputMap(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
// Frozen-A certified bilinear design from a K-ROM map and held-input records. Writes a new
// derived/koopman-bilinear-certificate archive; the K-ROM is never modified.
std::filesystem::path runBilinearCertificate(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
// Numerical-differentiation benchmark of a measured channel under measurement noise.
std::filesystem::path runDerivativeBenchmark(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
// Polynomial K-ROM fitted to a K-ROM's scattered (H-DMD) eigenfunction values: same
// eigenvalues, smooth evaluation. Writes derived/koopman-polynomial-fit; fills report.
std::filesystem::path fitPolynomialEigenfunctions(const std::filesystem::path& kromPath, std::size_t degree,
    const std::filesystem::path& root, json::Value::Object& report);
// Eigenfunction subset selection for echo-only KLQE -> KLQR (docs/architecture/MODE_SELECTION.md).
std::filesystem::path runModeSelection(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
// Offline lifting of a state target into a KROM's coordinates (time, values archive).
std::filesystem::path runControlTarget(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);
// Generic archive-derivation methods (KromEditing, ReferenceSignal, ControlScore, FigureDataset).
#define DYNSYSLAB_POOLED_METHOD(name) std::filesystem::path name( \
    const std::vector<std::filesystem::path>& sources, const json::Value& settings, \
    const std::filesystem::path& root, const std::string& label)
DYNSYSLAB_POOLED_METHOD(runKromRecenter);
DYNSYSLAB_POOLED_METHOD(runAnimation);
DYNSYSLAB_POOLED_METHOD(runKromReadout);
DYNSYSLAB_POOLED_METHOD(runReferenceSignal);
DYNSYSLAB_POOLED_METHOD(runControlScore);
DYNSYSLAB_POOLED_METHOD(runIoDelayModel);
DYNSYSLAB_POOLED_METHOD(runFigureDataset);
DYNSYSLAB_POOLED_METHOD(runBifurcationPoints);  // BifurcationPoints.cpp
DYNSYSLAB_POOLED_METHOD(runObservableSpectrumComparison);  // ObservableSpectrumComparison.cpp
#undef DYNSYSLAB_POOLED_METHOD
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
