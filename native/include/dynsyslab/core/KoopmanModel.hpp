#pragma once

#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/DiscreteInputMap.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/ScatteredEigenfunctions.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/KoopmanDatasets.hpp>
#include <dynsyslab/core/SpectralKROM.hpp>

#include <complex>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace dynsyslab::core {

struct ConstantControlMap {
    std::vector<std::complex<double>> matrix; // mode x input
    std::size_t modes{};
    std::size_t inputs{};
    std::string provenance;
    std::vector<std::complex<double>> discreteMatrix;
    double sampleTime{};
    std::string modelIdentity;
};

struct BilinearControlMap {
    std::vector<std::complex<double>> operators; // input x output mode x predictor mode
    std::vector<std::complex<double>> offsets;   // output mode x input
    std::vector<std::complex<double>> conjugateOperators;
    std::size_t modes{};
    std::size_t inputs{};
    std::string provenance;
    std::string modelIdentity;
};

using ApproximateControlMap = std::variant<ConstantControlMap, BilinearControlMap>;

struct KoopmanEigenfunctionDescriptor {
    std::size_t index{};
    std::complex<double> continuousEigenvalue;
    std::string representation;
    std::string gradientBackend;
    bool continuousEigenvalueAvailable{true};
};

class KoopmanModel {
public:
    // Coordinate-only models need neither an ODE nor fictitious physical state.
    explicit KoopmanModel(PolynomialEigenfunctionSet eigenfunctions,
        std::size_t inputDimension = 0, EigenfunctionFitDiagnostics diagnostics = {});
    KoopmanModel(PolynomialEigenfunctionSet eigenfunctions, SystemConfiguration system,
        EigenfunctionFitDiagnostics diagnostics = {});
    [[nodiscard]] const PolynomialEigenfunctionSet& polynomialEigenfunctions() const;
    [[nodiscard]] const EigenfunctionSpectrum& spectrum() const noexcept;
    [[nodiscard]] const EigenfunctionCoordinates& coordinates() const noexcept;
    [[nodiscard]] const EigenfunctionFitDiagnostics& fitDiagnostics() const noexcept;
    [[nodiscard]] bool hasPhysicalSystem() const noexcept { return bool(system_); }
    [[nodiscard]] std::optional<std::size_t> physicalStateDimension() const noexcept;
    [[nodiscard]] std::size_t eigenfunctionCount() const noexcept { return liftedDimension(); }
    [[nodiscard]] std::size_t realLiftedDimension() const noexcept { return splitDimension(); }
    [[nodiscard]] std::vector<std::complex<double>> evaluateFromState(const double* x) const { return evaluate(x); }
    [[nodiscard]] std::vector<std::complex<double>> evaluateFromCoordinates(const double* q) const { return evaluateCoordinates(q); }
    [[nodiscard]] std::vector<double> continuousGeneratorReal() const { return generatorSplit(); }
    [[nodiscard]] std::vector<std::complex<double>> stateJacobianControlMap(const double* x) const { return exactControlMap(x); }
    // Historical names below remain source-compatible adapters.
    // Convenience entry point for manually specified polynomials in the
    // system's state variables. Use KoopmanSINDyResult::fromCoefficients and
    // the existing constructor for named analytic or causal-history coordinates.
    [[nodiscard]] static KoopmanModel fromPolynomialCoefficients(
        SystemConfiguration system, std::size_t coordinateDimension,
        std::vector<std::size_t> exponents,
        std::vector<std::complex<double>> coefficients,
        std::vector<std::complex<double>> continuousEigenvalues,
        std::vector<double> center = {}, std::vector<double> scale = {});

    KoopmanModel(KoopmanSINDyResult eigenfunctions, SystemConfiguration system);
    KoopmanModel(ScatteredEigenfunctionData eigenfunctions,
                 SystemConfiguration system);

    [[nodiscard]] std::size_t stateDimension() const;
    [[nodiscard]] std::size_t coordinateDimension() const noexcept;
    [[nodiscard]] std::size_t liftedDimension() const noexcept;
    [[nodiscard]] std::size_t splitDimension() const noexcept;
    [[nodiscard]] std::size_t inputDimension() const;
    [[nodiscard]] const KoopmanSINDyResult& eigenfunctions() const;
    [[nodiscard]] const std::vector<std::complex<double>>&
        continuousEigenvalues() const;
    [[nodiscard]] std::vector<std::complex<double>> discreteEigenvalues(double dt = 0) const;
    [[nodiscard]] std::string eigenfunctionRepresentation() const;
    [[nodiscard]] std::vector<KoopmanEigenfunctionDescriptor>
        eigenfunctionDescriptors() const;
    [[nodiscard]] const ScatteredEigenfunctionData*
        scatteredEigenfunctions() const noexcept;
    [[nodiscard]] const SystemConfiguration& systemConfiguration() const;
    [[nodiscard]] const std::vector<KoopmanDatasetReference>& datasetReferences() const noexcept;
    void addDatasetReference(KoopmanDatasetReference reference);
    // Audit trail from an archive binding; never changes the fitted coefficients.
    [[nodiscard]] const std::string& bindingEvidenceJson() const noexcept { return bindingEvidenceJson_; }
    void setBindingEvidenceJson(std::string evidence);

    [[nodiscard]] std::vector<std::complex<double>> evaluate(const double* state) const;
    [[nodiscard]] std::vector<double> evaluateSplit(const double* state) const;
    [[nodiscard]] const Observable& coordinateObservable() const;
    [[nodiscard]] bool requiresExternalCoordinates() const noexcept;
    [[nodiscard]] std::vector<std::complex<double>> evaluateCoordinates(const double* coordinates) const;
    [[nodiscard]] std::vector<double> evaluateCoordinatesSplit(const double* coordinates) const;
    [[nodiscard]] std::vector<double> splitValues(const std::vector<std::complex<double>>& values) const;
    // Inverse of the KROM's real split layout; never reevaluates eigenfunctions.
    [[nodiscard]] std::vector<std::complex<double>> valuesFromSplit(const std::vector<double>& values) const;
    [[nodiscard]] std::vector<std::complex<double>> controlMapFromLifted(
        const std::vector<std::complex<double>>& values, const std::string& representation) const;
    [[nodiscard]] std::vector<double> realControlMap(
        const std::vector<std::complex<double>>& matrix) const;
    [[nodiscard]] std::vector<double> generatorSplit() const;
    [[nodiscard]] std::vector<std::complex<double>> exactControlMap(const double* state) const;
    [[nodiscard]] std::vector<std::complex<double>> controlMap(
        const double* state, const std::string& representation = "exact") const;
    [[nodiscard]] std::vector<std::complex<double>> controlMapCoordinates(
        const double* coordinates, const std::string& representation) const;
    [[nodiscard]] std::vector<double> controlMapSplit(
        const double* state, const std::string& representation = "exact") const;

    void setControlMap(std::string name, ApproximateControlMap representation);
    void setDiscreteInputMap(std::string name, DiscreteInputMap map);
    [[nodiscard]] const std::map<std::string, DiscreteInputMap>& discreteInputMaps() const noexcept {return discreteInputMaps_;}
    [[nodiscard]] bool hasControlMap(const std::string& name) const;
    [[nodiscard]] const ApproximateControlMap& storedControlMap(const std::string& name) const;
    [[nodiscard]] const std::map<std::string, ApproximateControlMap>& controlMaps() const noexcept;
    void setReconstruction(std::string observable, ReconstructionResult reconstruction);
    [[nodiscard]] bool hasReconstruction(const std::string& observable) const;
    [[nodiscard]] const ReconstructionResult& reconstruction(
        const std::string& observable) const;
    [[nodiscard]] std::vector<std::complex<double>> reconstruct(
        const std::string& observable, const double* state) const;
    [[nodiscard]] const std::map<std::string, ReconstructionResult>&
        reconstructions() const noexcept;
    [[nodiscard]] KoopmanModel extract(const std::vector<std::size_t>& indices) const;
    [[nodiscard]] const SpectralData& spectralData() const noexcept { return spectralData_; }
    void setSpectralData(SpectralData data) { spectralData_=spectralDataFromJson(spectralDataJson(data)); }

private:
    void attachPhysicalSystem(SystemConfiguration system);
    bool isRealMode(std::size_t mode) const noexcept;
    std::optional<PolynomialEigenfunctionSet> polynomial_;
    std::optional<EigenfunctionSpectrum> scatteredSpectrum_;
    EigenfunctionCoordinates scatteredCoordinates_;
    EigenfunctionFitDiagnostics fitDiagnostics_;
    std::size_t suppliedInputDimension_{};
    // Read-only legacy snapshot; runtime evaluation uses polynomial_.
    KoopmanSINDyResult eigenfunctions_;
    SystemConfiguration systemConfiguration_;
    std::shared_ptr<ODESystem> system_;
    std::shared_ptr<const Observable> coordinateObservable_;
    std::shared_ptr<const ScatteredEigenfunctionSet> scattered_;
    std::map<std::string, ApproximateControlMap> controlMaps_;
    std::map<std::string, DiscreteInputMap> discreteInputMaps_;
    std::map<std::string, ReconstructionResult> reconstructions_;
    std::vector<KoopmanDatasetReference> datasetReferences_;
    std::string bindingEvidenceJson_{"{}"};
    SpectralData spectralData_;
};

}  // namespace dynsyslab::core
