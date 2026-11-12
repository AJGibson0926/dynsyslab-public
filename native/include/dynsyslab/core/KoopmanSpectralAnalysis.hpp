#pragma once
#include <dynsyslab/core/SpectralKROM.hpp>
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
};
SpectralDictionary buildSpectralDictionary(const std::vector<double>& trainingPoints,
    std::size_t dimension,const SpectralDictionaryOptions& = {});
struct KoopmanSpectralOptions {
    std::string method{"resdmd"}; // edmd, resdmd, mpedmd, hermitian_dmd, tlsdmd, multdmd
    std::size_t rank{};
    double relativeTolerance{1e-10},sampleInterval{1};
    bool acknowledgeMeasurePreserving{},acknowledgeSelfAdjoint{};
    std::vector<std::complex<double>> residualGrid;
    std::string weightPolicy{"uniform"}; // uniform or smooth-birkhoff
};
RkhsAdjointSpectrum rkhsAdjointSpectrum(const KoopmanSnapshotPairs&,double gaussianBandwidth,
    std::size_t rank=0,double tolerance=1e-10,const std::vector<std::complex<double>>& grid={});
SpectralKROM fitKoopmanSpectral(const KoopmanSnapshotPairs&,const SpectralDictionary&,
                               const KoopmanSpectralOptions& = {});
SpectralKROM fitStochasticKoopmanSpectral(const KoopmanSnapshotPairs&,const SpectralDictionary&,
                                         const KoopmanSpectralOptions&);
StochasticSpectralResiduals stochasticKoopmanResiduals(const SpectralKROM&,const KoopmanSnapshotPairs&);
SpectralResiduals koopmanResiduals(const SpectralKROM&,const KoopmanSnapshotPairs&,
    const std::vector<std::complex<double>>& grid={},std::string scope="empirical-held-out");
SpectralInvariance koopmanInvariance(const std::vector<double>& x,const std::vector<double>& y,
    std::size_t samples,std::size_t features,double tolerance=1e-10);
void fitSpectralReconstruction(SpectralKROM&,const std::string& name,
    const std::vector<double>& points,const std::vector<double>& outputs,std::size_t outputDimension,
    const std::vector<double>& weights={});
// Operator rollout works even for defective matrices; it uses no future measurements.
std::vector<double> rolloutSpectralKrom(const SpectralKROM&,const std::string& reconstruction,
    const std::vector<double>& initialPoint,std::size_t steps);
SpectralMeasure correlationSpectralMeasure(const std::vector<std::complex<double>>& signal,
    std::size_t lags,const std::vector<double>& angles,bool acknowledgeStationarity);
void riggedSpectralMeasure(SpectralKROM&,const std::vector<double>& observableCoefficients,
    const std::string& name,const std::vector<double>& angles,double bandwidth,std::size_t order=1);
std::vector<double> koopmanPageRank(const std::vector<double>& transition,std::size_t size,
    const std::vector<double>& personalization,double damping=0.85,double tolerance=1e-12);
} // namespace dynsyslab::core
