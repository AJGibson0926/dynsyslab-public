#pragma once

#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/Eigenpair.hpp>
#include <dynsyslab/core/DiscreteInputMap.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/KoopmanDatasets.hpp>
#include <dynsyslab/core/KoopmanCompression.hpp>
#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <complex>
#include <filesystem>
#include <iterator>
#include <optional>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
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

// The name of a domain: "state" for the full state, "window" for series windows, else its coordinate
// names joined by commas, else "coordinates".
[[nodiscard]] std::string domainName(const FunctionDomain& domain);

struct EigenfunctionErrors {
    std::string representation, domain;
    std::size_t samples{};
    std::vector<double> valueRelative;    // ||R(p_k) - psi_k|| / ||psi_k|| over the samples, per mode
    std::vector<double> koopmanResidual;  // ||R(p_k+1) - mu^dt R(p_k)|| / ||R(p_k)||, per mode
};

// Empirical spectral weight of each eigenfunction over its training samples: w_kj = |psi_j(x_k)|^2 /
// sum_k |psi_j(x_k)|^2 (samples x modes; uniform for a mode that vanishes). An exploratory
// localization of where each mode lives along the data, not an identified spectral measure (for
// those, see KoopmanSpectralAnalysis.hpp; Colbrook & Townsend, Comm. Pure Appl. Math. 77, 2024).
[[nodiscard]] std::vector<double> empiricalSpectralWeight(const EigenfunctionSamples& samples);

// The krom a model's inherited parts come from: its archive (path and SHA-256) and which of its pairs,
// by stable mode id (Spectrum::modeId), in this model's order. A krom written from such a model
// references it as its parent and stores only what the model added (KoopmanArchive.hpp).
struct KromParent {
    std::filesystem::path path;
    std::string sha256;
    std::vector<std::size_t> modes;  // mode ids
};

struct KoopmanEigenfunctionDescriptor {
    std::size_t index{};
    std::complex<double> continuousEigenvalue;
    std::string representation;
    std::string gradientBackend;
    bool continuousEigenvalueAvailable{true};
};

// A Koopman reduced-order model (K-ROM): a plant for estimation and control. It consists of a
// Spectrum (Spectrum.hpp: the eigenpairs with their eigenfunctions, the continuous part, the
// observables' decomposition over both, the pseudospectrum and the operator), input maps B_psi (how
// inputs drive the eigenfunctions), and measurements: observables g read out through their Koopman
// modes (g ~ offset + C psi), which the spectrum stores. Without B or g it is still a model; they can be
// fitted later. The model adds what a plant needs beyond the spectrum: the evaluator of the first
// eigenfunction representation (useRepresentation), the physical system it is bound to, dataset
// references, the fit's diagnostics and its parent krom. The archive is the "krom" data object
// (KoopmanArchive.hpp).
class KoopmanModel {
public:
    using Representation = EigenfunctionRepresentation;  // every eigenfunction: one output per mode
    // A K-ROM from data: the spectrum, the eigenfunctions' sample values (may hold no samples), their
    // origin and representations (named by representationName; the first is evaluated). No physical
    // system: bindSystem attaches one.
    KoopmanModel(Spectrum spectrum, EigenfunctionSamples samples, DataOrigin origin,
                 std::vector<std::shared_ptr<const FunctionRepresentation>> representations,
                 std::size_t inputDimension = 0);
    // A K-ROM on a spectrum that already holds its eigenfunctions (at least one representation).
    explicit KoopmanModel(Spectrum spectrum, std::size_t inputDimension);
    // Coordinate-only models need neither an ODE nor fictitious physical state.
    // The eigenfunctions' shared expansion and, separately, their eigenvalues (one per output).
    KoopmanModel(FunctionExpansion functions, Spectrum spectrum,
        std::size_t inputDimension = 0, EigenfunctionFitDiagnostics diagnostics = {});
    KoopmanModel(FunctionExpansion functions, Spectrum spectrum, SystemConfiguration system,
        EigenfunctionFitDiagnostics diagnostics = {});
    explicit KoopmanModel(EigenfunctionFitResult fit, std::size_t inputDimension = 0);
    KoopmanModel(EigenfunctionFitResult fit, SystemConfiguration system);
    [[nodiscard]] const FunctionExpansion& expansion() const;
    // Mode `mode` as an Eigenfunction object on the physical state (or on the K-ROM coordinates
    // when they cannot be computed from a state): its representation is this model's expansion or
    // interpolation, its eigenvalue the model's.
    // Every representation of mode `mode` (the evaluated one first), each through its domain's state
    // map (`catalog` binds a catalog map; default: the system's); a function of q when a domain has none.
    [[nodiscard]] Eigenfunction eigenfunction(std::size_t mode, const ObservableCatalog* catalog = nullptr) const;
    // Eigenpair `mode`: its eigenvalue, eigenfunction (every representation, as eigenfunction()) and residual.
    [[nodiscard]] Eigenpair eigenpair(std::size_t mode, const ObservableCatalog* catalog = nullptr) const;
    // Every eigenpair in mode order, for (const auto& pair : krom.pairs()). Each is assembled on
    // dereference by eigenpair(j); the range stores and copies nothing (it must not outlive the model).
    class EigenpairRange;
    [[nodiscard]] EigenpairRange pairs(const ObservableCatalog* catalog = nullptr) const;

    // The eigenfunctions as data.
    [[nodiscard]] const EigenfunctionSamples& samples() const noexcept { return spectrum_->samples(); }
    [[nodiscard]] const DataOrigin& origin() const noexcept { return spectrum_->origin(); }
    void setOrigin(DataOrigin origin) { spectrum_->withOrigin(std::move(origin)); }
    // Sets the values at the training samples (a model fitted without them); an interpolation keeps its own.
    void setSamples(EigenfunctionSamples samples);
    [[nodiscard]] const std::vector<Representation>& representations() const noexcept { return spectrum_->representations(); }
    [[nodiscard]] std::vector<std::string> representationNames() const;
    // By name, else by kind when exactly one representation has it; empty selects the evaluated one.
    // An ambiguous kind is refused with the names to choose from.
    [[nodiscard]] const Representation& representation(const std::string& key = {}) const;
    // The name of `key`'s domain, checked against the observables of data handed in (empty: no check).
    [[nodiscard]] std::string recognizeDomain(const std::vector<std::string>& observables, const std::string& key = {}) const;
    // Every eigenfunction at a point of the representation's own domain.
    [[nodiscard]] std::vector<std::complex<double>> evaluateRepresentation(const double* point, const std::string& key = {}) const;
    // Errors at the training samples. An interpolation evaluates at its base points; a window map at
    // the windows of `input` (the record the eigenfunctions were fitted on).
    [[nodiscard]] EigenfunctionErrors errors(const std::string& key = {}, const TimeSeriesRecord* input = nullptr) const;
    // Replaces the representation of the same name, keeping its position.
    void replaceRepresentation(std::shared_ptr<const FunctionRepresentation> function);
    // Adds a representation (a new name); `first` makes it the evaluated one.
    void addRepresentation(std::shared_ptr<const FunctionRepresentation> function, bool first = false);
    // Evaluates through the named representation (moved to the front). An interpolation on catalog
    // observables binds them through the system; a window map evaluates only on windows.
    void useRepresentation(const std::string& key);
    // The interpolation, and switching its interpolant in place (sweeps). The domain and points stay.
    [[nodiscard]] const FunctionInterpolation& interpolation() const;
    void setInterpolationMethod(InterpolationMethod method, std::size_t neighbors);
    // Attaches a physical system to a model built from data (its catalog binds the domains).
    void bindSystem(SystemConfiguration system);
    [[nodiscard]] const Spectrum& spectrum() const noexcept;
    [[nodiscard]] const FunctionDomain& coordinates() const noexcept;
    [[nodiscard]] const EigenfunctionFitDiagnostics& fitDiagnostics() const noexcept;
    [[nodiscard]] bool hasPhysicalSystem() const noexcept { return bool(system_); }
    [[nodiscard]] std::optional<std::size_t> physicalStateDimension() const noexcept;
    [[nodiscard]] std::size_t eigenfunctionCount() const noexcept { return liftedDimension(); }
    [[nodiscard]] std::size_t realLiftedDimension() const noexcept { return splitDimension(); }
    [[nodiscard]] std::vector<std::complex<double>> evaluateFromState(const double* x) const { return evaluate(x); }
    [[nodiscard]] std::vector<std::complex<double>> evaluateFromCoordinates(const double* q) const { return evaluateCoordinates(q); }
    [[nodiscard]] std::vector<double> continuousGeneratorReal() const { return generatorSplit(); }
    [[nodiscard]] std::vector<std::complex<double>> stateJacobianControlMap(const double* x) const { return exactControlMap(x); }
    // Manually specified polynomials in the system's state variables (polynomialEigenfunctions in
    // EigenfunctionFit.hpp builds them on other domains).
    [[nodiscard]] static KoopmanModel fromPolynomialCoefficients(
        SystemConfiguration system, std::size_t coordinateDimension,
        std::vector<std::size_t> exponents,
        std::vector<std::complex<double>> coefficients,
        std::vector<std::complex<double>> continuousEigenvalues,
        std::vector<double> center = {}, std::vector<double> scale = {});

    KoopmanModel(InterpolationData eigenfunctions, Spectrum spectrum,
                 SystemConfiguration system);

    [[nodiscard]] std::size_t stateDimension() const;
    [[nodiscard]] std::size_t coordinateDimension() const noexcept;
    [[nodiscard]] std::size_t liftedDimension() const noexcept;
    [[nodiscard]] std::size_t splitDimension() const noexcept;
    [[nodiscard]] std::size_t inputDimension() const;
    // Additive-state actuation xdot = f(x) + B u (B row-major n x m): replaces the bound
    // system's input matrix in inputDimension() and the exact map grad psi(x) B, so KLQR and
    // KMPC can drive m inputs (Korda & Mezic 2018, Automatica 93:149-160). Empty restores it.
    void setAdditiveInputMatrix(std::vector<double> matrix);
    [[nodiscard]] const std::vector<double>& additiveInputMatrix() const noexcept { return additiveInput_; }
    [[nodiscard]] const std::vector<std::complex<double>>&
        continuousEigenvalues() const;
    [[nodiscard]] std::vector<std::complex<double>> discreteEigenvalues(double dt = 0) const;
    [[nodiscard]] std::string eigenfunctionRepresentation() const;
    [[nodiscard]] std::vector<KoopmanEigenfunctionDescriptor>
        eigenfunctionDescriptors() const;
    [[nodiscard]] const InterpolationData*
        scatteredEigenfunctions() const noexcept;
    [[nodiscard]] const SystemConfiguration& systemConfiguration() const;
    [[nodiscard]] const std::vector<KoopmanDatasetReference>& datasetReferences() const noexcept;
    void addDatasetReference(KoopmanDatasetReference reference);
    // Audit trail from an archive binding; never changes the fitted coefficients.
    [[nodiscard]] const std::string& bindingEvidenceJson() const noexcept { return bindingEvidenceJson_; }
    void setBindingEvidenceJson(std::string evidence);

    [[nodiscard]] std::vector<std::complex<double>> evaluate(const double* state) const;
    [[nodiscard]] std::vector<double> evaluateSplit(const double* state) const;
    // Row-major complex mode x physical-state Jacobian. History banks require
    // an explicit history-state provider and cannot use this instantaneous map.
    [[nodiscard]] std::vector<std::complex<double>> stateJacobian(const double* state) const;
    [[nodiscard]] std::vector<double> stateJacobianSplit(const double* state) const;
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
    // Chain rule for coordinate KROMs: B[mode,input] = sum_c dpsi_mode/dq_c(q) *
    // sensitivity[c*inputs+input], where sensitivity = d(qdot)/du supplied by the
    // caller (e.g. a history transform's newest-sample weight times L_g y).
    [[nodiscard]] std::vector<std::complex<double>> controlMapFromCoordinateSensitivity(
        const double* coordinates, const std::vector<double>& sensitivity) const;
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
    // Measurements: observables g read out through their Koopman modes C (g ~ offset + C psi), which the
    // spectrum stores (Spectrum::koopmanModes). Setting one makes it this model's own.
    void setReconstruction(std::string observable, ReconstructionResult reconstruction);
    [[nodiscard]] bool hasReconstruction(const std::string& observable) const;
    [[nodiscard]] const ReconstructionResult& reconstruction(
        const std::string& observable) const;
    [[nodiscard]] std::vector<std::complex<double>> reconstruct(
        const std::string& observable, const double* state) const;
    [[nodiscard]] const std::map<std::string, ReconstructionResult>&
        reconstructions() const noexcept;
    [[nodiscard]] KoopmanModel extract(const std::vector<std::size_t>& indices) const;
    // The spectrum's evidence (pair residuals, continuous part, pseudospectrum; Spectrum.hpp) and the
    // fit's diagnostics. Setting any of them makes it this model's own (a child krom stores it).
    void setResidual(std::size_t mode, std::optional<EigenpairResidual> residual);
    void setSpectralEvidence(const Spectrum& evidence);
    [[nodiscard]] const SpectralDiagnostics& spectralDiagnostics() const noexcept { return spectralDiagnostics_; }
    void setSpectralDiagnostics(SpectralDiagnostics diagnostics) { spectralDiagnostics_ = std::move(diagnostics); ownSpectralEvidence_ = true; }
    // The spectrum's compressed operator (Spectrum::compression), when the eigenpairs were derived from
    // one: what Rigged DMD, residuals on new data and operator rollout need. Inherited by child kroms
    // and kept by extraction (it is the operator, not the selection).
    [[nodiscard]] const KoopmanCompression* compression() const noexcept { return spectrum_->compression(); }
    void setCompression(KoopmanCompression compression);
    [[nodiscard]] bool ownsCompression() const noexcept { return ownCompression_; }
    // The spectrum's functions on spectral sets (Spectrum::spectralFunction): the pseudospectrum's
    // minimizers and the measures' wave packets. Setting one makes it this model's own.
    [[nodiscard]] const FunctionExpansion* spectralFunction(const std::string& key) const { return spectrum_->spectralFunction(key); }
    void setSpectralFunction(const std::string& key, FunctionExpansion function);
    [[nodiscard]] std::vector<std::string> spectralFunctionKeys() const { return spectrum_->spectralFunctionKeys(); }
    [[nodiscard]] bool ownsSpectralFunction(const std::string& key) const { return !inherited_.spectralFunctions.count(key); }

    // The krom this model was loaded from (set by the loader; extract keeps it with the chosen modes).
    // Everything the model has when it is set is inherited; parts added or replaced afterwards are its
    // own, and a krom written from the model stores only those, referencing the parent.
    [[nodiscard]] const std::optional<KromParent>& parent() const noexcept { return parent_; }
    void setParent(KromParent parent);
    // Forgets the parent: the model is written whole.
    void clearParent();
    [[nodiscard]] bool ownsRepresentation(const std::string& name) const { return !inherited_.representations.count(name); }
    [[nodiscard]] bool ownsReconstruction(const std::string& name) const { return !inherited_.readouts.count(name); }
    [[nodiscard]] bool ownsControlMap(const std::string& name) const { return !inherited_.controlMaps.count(name); }
    [[nodiscard]] bool ownsDiscreteInputMap(const std::string& name) const { return !inherited_.discreteMaps.count(name); }
    [[nodiscard]] bool ownsSpectralEvidence() const noexcept { return ownSpectralEvidence_; }

private:
    void attachPhysicalSystem(SystemConfiguration system);
    bool isRealMode(std::size_t mode) const noexcept;
    std::optional<FunctionExpansion> polynomial_;
    std::optional<Spectrum> spectrum_;
    EigenfunctionFitDiagnostics fitDiagnostics_;
    std::size_t suppliedInputDimension_{};
    std::vector<double> additiveInput_;
    // Read-only legacy snapshot; runtime evaluation uses polynomial_.
    SystemConfiguration systemConfiguration_;
    std::shared_ptr<ODESystem> system_;
    std::shared_ptr<const FunctionInterpolation> scattered_;
    std::map<std::string, ApproximateControlMap> controlMaps_;
    std::map<std::string, DiscreteInputMap> discreteInputMaps_;
    std::vector<KoopmanDatasetReference> datasetReferences_;
    std::string bindingEvidenceJson_{"{}"};
    struct Inherited { std::set<std::string> representations, readouts, controlMaps, discreteMaps, spectralFunctions; };
    void activate();  // builds the evaluator from representations_.front()
    void restampMaps();
    std::optional<KromParent> parent_;
    Inherited inherited_;
    bool ownSpectralEvidence_{false};
    SpectralDiagnostics spectralDiagnostics_;
    bool ownCompression_{false};
};

class KoopmanModel::EigenpairRange {
public:
    class iterator {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = Eigenpair;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = Eigenpair;
        iterator(const KoopmanModel* model, std::size_t mode, const ObservableCatalog* catalog)
            : model_(model), mode_(mode), catalog_(catalog) {}
        Eigenpair operator*() const { return model_->eigenpair(mode_, catalog_); }
        iterator& operator++() { ++mode_; return *this; }
        iterator operator++(int) { auto previous = *this; ++mode_; return previous; }
        bool operator==(const iterator& other) const { return mode_ == other.mode_ && model_ == other.model_; }
        bool operator!=(const iterator& other) const { return !(*this == other); }
    private:
        const KoopmanModel* model_;
        std::size_t mode_;
        const ObservableCatalog* catalog_;
    };
    EigenpairRange(const KoopmanModel& model, const ObservableCatalog* catalog) : model_(&model), catalog_(catalog) {}
    [[nodiscard]] iterator begin() const { return {model_, 0, catalog_}; }
    [[nodiscard]] iterator end() const { return {model_, model_->liftedDimension(), catalog_}; }
    [[nodiscard]] std::size_t size() const { return model_->liftedDimension(); }
private:
    const KoopmanModel* model_;
    const ObservableCatalog* catalog_;
};

inline KoopmanModel::EigenpairRange KoopmanModel::pairs(const ObservableCatalog* catalog) const {
    return {*this, catalog};
}


}  // namespace dynsyslab::core
