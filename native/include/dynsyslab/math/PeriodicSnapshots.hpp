#pragma once
#include <dynsyslab/math/HarmonicAnalysis.hpp>
namespace dynsyslab::math {
struct PeriodicSnapshotOrbit {std::vector<FourierSeries> channels;double weight{1};};
struct PeriodicSnapshots {
    std::size_t samples{},dimension{};
    std::vector<double> x,y,weights;
};
// One orbit's delay dictionary on a uniform phase grid, row-major phase/feature.
// Feature order is offset/channel/(real,imag). Useful for bounded row-block fits.
std::vector<double> periodicDelaySamples(const std::vector<FourierSeries>&channels,
    std::size_t phaseCount,const std::vector<double>&timeOffsets,
    bool splitComplex=false,double phaseOrigin=0);
// Product quadrature: explicit orbit masses times uniform phase. Row order is
// orbit then phase; columns are delay then channel then real/imag (if split).
// A delay here is a SIGNED time offset: negative means past, positive future.
// No centering, normalization, interpolation or fitted dynamics are introduced.
PeriodicSnapshots periodicSnapshotPairs(const std::vector<PeriodicSnapshotOrbit>&orbits,
    std::size_t phaseCount,const std::vector<double>&timeOffsets,double sampleInterval,
    bool splitComplex=false,double phaseOrigin=0);

// Fourier coefficients of sum_{tap,channel} c[tap,channel] g_channel(t+offset).
// No phase sampling/aliasing; applies to complex readouts and complex channels.
FourierSeries periodicDelayReadout(const std::vector<FourierSeries>&channels,
    const std::vector<double>&timeOffsets,const std::vector<std::complex<double>>&coefficients);

struct PeriodicQuadratureAudit {
    std::size_t channels{},maximumHarmonic{},aliasFreePhaseCount{};
    double suppliedWeight{};
    std::vector<std::size_t> phaseCounts;
    std::vector<double> exactMass,sampledMass; // normalized weights; count-major/channel
};
// Checks second moments against Parseval. N>2*maximumHarmonic also resolves
// every pairwise product of these finite polynomials, including shifted copies.
// It does not certify omitted harmonics, orbit quadrature, or invariance.
PeriodicQuadratureAudit auditPeriodicQuadrature(const std::vector<PeriodicSnapshotOrbit>&,
    const std::vector<std::size_t>&phaseCounts,double phaseOrigin=0);
}
