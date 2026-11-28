#pragma once
// Spectral Koopman fits and the analyses that need their compressed operator. Every fit is a krom
// (KoopmanModel): eigenvalues on its Spectrum, eigenfunctions as an expansion on the frozen
// dictionary, residuals on the eigenpairs, and the operator itself as compression() (Williams,
// Kevrekidis & Rowley 2015; Colbrook & Townsend 2024; Colbrook 2023).
#include <dynsyslab/core/KoopmanCompression.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <optional>
namespace dynsyslab::core {
struct KoopmanSnapshotPairs {
    std::vector<double> x,y; // samples x coordinates; never crosses a trajectory boundary
    std::size_t samples{},dimension{};
    std::vector<double> weights; // optional positive quadrature weights, normalized internally
    // Optional independent successor realizations, replicate x sample x coordinate.
    std::vector<double> successorReplicates;
    std::size_t replicateCount{};
};
struct SpectralDictionaryOptions {
    std::string kind{"linear"};
    std::size_t degree{2},centers{32},maximumFeatures{4096};
    bool constant{true},normalize{};
    double bandwidth{1};
    std::string indexSet{"total_degree"}, centerSelection{"farthest"};
    std::size_t angularOrder{4}, hermiteOrder{4}, centerIterations{100};
    double smoothness{2}, amplitude{1};
    std::vector<double> center, scale; // explicit frozen affine coordinates
};
SpectralDictionary buildSpectralDictionary(const std::vector<double>& trainingPoints,
    std::size_t dimension,const SpectralDictionaryOptions& = {});
struct KoopmanSpectralOptions {
    std::string method{"resdmd"}; // edmd, resdmd, mpedmd, hermitian_dmd, tlsdmd, multdmd, stochastic_resdmd
    std::size_t rank{};
    double relativeTolerance{1e-10},sampleInterval{1};
    bool acknowledgeMeasurePreserving{},acknowledgeSelfAdjoint{};
    std::vector<std::complex<double>> residualGrid;
    std::string weightPolicy{"uniform"}; // uniform or smooth-birkhoff
};
RkhsAdjointSpectrum rkhsAdjointSpectrum(const KoopmanSnapshotPairs&,double gaussianBandwidth,
    std::size_t rank=0,double tolerance=1e-10,const std::vector<std::complex<double>>& grid={});
// A spectral fit as a krom: its pairs numbered 0..rank-1 (discrete-time multipliers), its
// eigenfunctions an expansion on the frozen dictionary over `domain` (default: the data coordinates),
// one-step residuals on the pairs and, on options.residualGrid, the pseudospectrum with its minimizers
// as spectralFunction(pseudospectrumFunctionKey()); invariance diagnostics; the compression.
KoopmanModel fitKoopmanSpectral(const KoopmanSnapshotPairs&,const SpectralDictionary&,
                                const KoopmanSpectralOptions& = {},std::optional<FunctionDomain> domain = std::nullopt);
KoopmanModel fitStochasticKoopmanSpectral(const KoopmanSnapshotPairs&,const SpectralDictionary&,
                                          const KoopmanSpectralOptions&,std::optional<FunctionDomain> domain = std::nullopt);
StochasticSpectralResiduals stochasticKoopmanResiduals(const KoopmanModel&,const KoopmanSnapshotPairs&);
// The model's spectrum with each pair's one-step residual on `data` (Eigenpair.hpp) and, on `grid`,
// the pseudospectrum. The grid's minimizers (grid x features, on the frozen dictionary) go to
// *minimizers when given. Needs the model's compression.
Spectrum koopmanResiduals(const KoopmanModel&,const KoopmanSnapshotPairs&,
    const std::vector<std::complex<double>>& grid={},std::string scope="empirical-held-out",
    std::vector<std::complex<double>>* minimizers=nullptr);
SpectralInvariance koopmanInvariance(const std::vector<double>& x,const std::vector<double>& y,
    std::size_t samples,std::size_t features,double tolerance=1e-10);
// A least-squares readout of dictionary features named `name`, added to the model's compression.
void fitFeatureReadout(KoopmanModel&,const std::string& name,
    const std::vector<double>& points,const std::vector<double>& outputs,std::size_t outputDimension,
    const std::vector<double>& weights={});
// Rollout through the compressed operator: works for defective matrices and uses no future data.
std::vector<double> rolloutCompression(const KoopmanCompression&,const std::string& readout,
    const std::vector<double>& initialPoint,std::size_t steps);
// Independent initial points, sample-major input; output is trajectory x (steps+1) x readout
// output. The readout is factored once.
std::vector<double> rolloutCompressionBatch(const KoopmanCompression&,const std::string& readout,
    const std::vector<double>& initialPoints,std::size_t steps);
// The Fejer estimate of a stationary signal's spectral measure from its biased autocorrelation (Korda,
// Putinar & Mezic, Appl. Comput. Harmon. Anal. 48, 2020): band "correlation" (continuous, on the
// unit circle) and the measure of observable "signal" on it, with its moments.
SpectralMeasureEstimate correlationSpectralMeasure(const std::vector<std::complex<double>>& signal,
    std::size_t lags,const std::vector<double>& angles,bool acknowledgeStationarity);
// True when the compressed action is orthogonal to within `tolerance` (A^T A = I entrywise): a
// measure-preserving compression, on which Rigged DMD's unit-circle resolvent applies.
bool compressionIsUnitary(const KoopmanCompression&,double tolerance=1e-8);
// The model's eigenfunctions on its compression's frozen dictionary: the expansion representation
// with one term per dictionary feature (throws when there is none).
const FunctionExpansion& frozenDictionaryExpansion(const KoopmanModel&);
// The Koopman modes of an observable given on the frozen dictionary (dictionary x outputs, row-major):
// its exact coefficients on the krom's eigenfunctions (the observable must lie in their span), stored
// as the observable's modes on the krom's spectrum (KoopmanModel::setReconstruction).
void setDictionaryObservableModes(KoopmanModel&,const std::string& observable,
    const std::vector<double>& coefficients,std::size_t outputs);
// Rigged DMD (Colbrook, Drysdale & Horning 2024) of output `output` of observable `observable`, given as
// coefficients on the frozen dictionary, on spectral band `band` (unit circle, the angle
// grid, the resolvent kernel of `bandwidth` and `order`; added, or matched when it exists). Adds the
// observable's measure (density, total mass) with its generalized eigenfunctions (wave packets) and,
// when the observable has no Koopman modes yet, its modes (output 0 of a new observable), from which
// the measure's atoms follow (Spectrum::atoms). Requires a compression that compressionIsUnitary.
void riggedSpectralMeasure(KoopmanModel&,const std::vector<double>& observableCoefficients,
    const std::string& observable,std::size_t output,const std::string& band,
    const std::vector<double>& angles,double bandwidth,std::size_t order=1);
std::vector<double> koopmanPageRank(const std::vector<double>& transition,std::size_t size,
    const std::vector<double>& personalization,double damping=0.85,double tolerance=1e-12);
} // namespace dynsyslab::core
