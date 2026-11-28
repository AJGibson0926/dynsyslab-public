#pragma once
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/math/PeriodicSnapshots.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>
namespace dynsyslab::core {
struct PeriodicSpectralReference {
    std::size_t channel{},validationOrbits{};
    double sampleInterval{},bandwidth{};std::size_t order{};
    std::vector<double> angles,centers,packetRms; // center-major/orbit
    math::SmoothedOrbitSpectrum spectrum;
};
struct PeriodicSpectralComparison {
    std::size_t rank{},features{};
    double fittedMass{},densityErrorNorm{},densityReferenceNorm{};
    std::vector<double> density,packetRms,packetErrorRms;
    std::vector<double> amplitudeErrorNorm,complexErrorNorm,packetReferenceNorm;
};
struct PeriodicDelayForecastAudit {
    std::size_t steps{},channels{};
    std::vector<std::size_t> stepIndices; // sorted native steps, starts at zero
    std::vector<double> referenceMeanSquare; // channel
    std::vector<double> predictionMeanSquare,errorMeanSquare; // time/channel, includes t=0
    std::vector<double> persistenceErrorMeanSquare; // hold the true current value fixed
};
struct PeriodicDelayPacketFeatures {
    std::vector<std::size_t> channels;
    std::vector<double> timeOffsets,centers,coefficients;
    std::size_t outputs{}; // real/imaginary pairs, coefficients input x output
};
// Resolve the newest-coordinate observable, compute frozen Rigged-DMD packets,
// and shift the history origin to its newest sample. No source archive changes.
PeriodicDelayPacketFeatures periodicDelayPacketFeatures(const KoopmanModel&,
    std::size_t channel,const std::vector<double>& centers,double bandwidth,std::size_t order);
// Newest-coordinate, past-only free forecasts; exact uniform-phase integrals
// for supplied finite Fourier series. Reads the frozen delay contract. No
// phase grid, refitting, temporal interpolation or Fourier-tail certificate.
std::vector<PeriodicDelayForecastAudit> auditPeriodicDelayForecast(const KoopmanModel&,
    const std::vector<math::PeriodicSnapshotOrbit>&validation,
    const std::vector<std::size_t>&stepsPerOrbit,
    const std::vector<std::vector<std::size_t>>&evaluationSteps={});
// Fitting owns no system equations: inputs are supplied physical-time Fourier
// orbits. Uniform phase products must be resolved unless explicitly opted out.
// Storage: auto, dense, streaming_qr. Streaming preserves the weighted joint
// source/successor least-squares problem; it does not subsample or truncate.
// General delay fitting accepts EDMD/ResDMD/mpEDMD; packet comparison below
// still requires mpEDMD. No spectrum is smoothed by this fitting operation.
KoopmanModel fitPeriodicDelayModel(const std::vector<math::PeriodicSnapshotOrbit>&,
    const std::vector<std::size_t>&channels,const std::vector<double>&timeOffsets,
    std::size_t phaseCount,const KoopmanSpectralOptions&,bool requireResolvedQuadrature=true,
    const std::string& storage="auto");
KoopmanModel fitPeriodicDelaySpectrum(const std::vector<math::PeriodicSnapshotOrbit>&,
    const std::vector<std::size_t>&channels,const std::vector<double>&timeOffsets,
    std::size_t phaseCount,const KoopmanSpectralOptions&,bool requireResolvedQuadrature=true,
    const std::string& storage="auto");
PeriodicSpectralReference periodicSpectralReference(
    const std::vector<math::PeriodicSnapshotOrbit>&training,
    const std::vector<math::PeriodicSnapshotOrbit>&validation,std::size_t channel,
    double sampleInterval,const std::vector<double>&angles,const std::vector<double>&centers,
    double bandwidth,std::size_t order);
// Frozen model, same signed offsets/channel order used at fit. Current-time
// coordinate readout is compared to the ORIGINAL observable, not its projection.
PeriodicSpectralComparison comparePeriodicDelaySpectrum(const KoopmanModel&,
    const std::vector<std::size_t>&channels,const std::vector<double>&timeOffsets,
    const std::vector<math::PeriodicSnapshotOrbit>&validation,const PeriodicSpectralReference&);
}
