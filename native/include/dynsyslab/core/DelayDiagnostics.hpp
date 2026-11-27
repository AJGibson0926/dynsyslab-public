#pragma once
// Delay-coordinate building blocks and checks (docs/architecture/IDENTIFICATION.md, "Delay
// diagnostics"). Delay coordinates: Takens, Lect. Notes Math. 898:366-381 (1981); Hankel/SVD
// coordinates: Broomhead & King, Physica D 20:217-236 (1986).
//
//   scalarDelayPairs     two taps per lag from data available at each end index: the older and
//                        newer samples p[e - l - 2h], p[e - l], their midpoint, the secant slope
//                        (newer - older) / (2 h dt) and the centred second difference. Every
//                        estimate refers to e - l - h, so it carries a latency of h dt.
//   delayConsistency     the relative error between each delayed tap (value, rate) of a forecast and
//                        the current tap shifted by its lag, read from the cubic Hermite interpolant
//                        of the current (value, rate) pair (Fritsch & Carlson, SIAM J. Numer. Anal.
//                        17:238-246, 1980); scaled by the reference tap norms.
//   extendLockedBasis    a frozen prefix of eigenpairs (U, s) of a window Gram matrix G stays
//                        bit-identical; the basis is extended with the leading eigenvectors of the
//                        complement-projected Gram P G P, P = I - U U^T (exact under repeated
//                        eigenvalues), signs fixed so each new column's largest entry is positive.
//                        Checks the prefix: U^T U = I and ||G u - s^2 u|| / s^2 <= 1e-6 per column.
//   historyEndpointDecoder  the exact linear map from history coordinates to the newest sample of a
//                        channel: y_end = d . xi, with d_i = U[end, i] s_i / gain (gain 1 for the
//                        native xi = U^T w / s; sqrt(windows) for unit training-RMS coordinates).
//                        The same map takes coordinate rates to the endpoint rate.
#include <dynsyslab/core/HankelOperator.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct DelayPairs {
    std::vector<double> raw;        // ends x (2 lags): (older, newer) per lag
    std::vector<double> state;      // ends x (2 lags): (midpoint, secant) per lag
    std::vector<double> curvature;  // ends x lags: centred second difference
    std::size_t ends{}, lags{};
    double latency{};               // half_width * dt
};
[[nodiscard]] DelayPairs scalarDelayPairs(const std::vector<double>& signal, const std::vector<std::size_t>& ends,
                                          const std::vector<std::size_t>& lags, double dt, std::size_t halfWidth = 1);

struct DelayConsistency {
    std::vector<double> lags;               // the delayed taps (lags[1:])
    std::vector<double> value, rate;        // relative errors per delayed tap; NaN when undefined
};
// time: increasing clock; prediction, reference: samples x (2 taps), (value, rate) per tap, tap 0 at lag 0.
[[nodiscard]] DelayConsistency delayConsistency(const std::vector<double>& time, const std::vector<double>& prediction,
                                                const std::vector<double>& reference, const std::vector<double>& lags,
                                                std::size_t stride = 1);

struct LockedBasis {
    std::vector<double> modes;           // dimension x rank, row-major
    std::vector<double> singularValues;  // rank
    double gramError{};                  // max |U^T G U / (s s^T) - I|
};
[[nodiscard]] LockedBasis extendLockedBasis(const std::vector<double>& gram, std::size_t dimension,
                                            const std::vector<double>& lockedModes,
                                            const std::vector<double>& lockedSingularValues, std::size_t rank);

// rows = taps * channels, oldest tap first, channels innermost.
[[nodiscard]] std::vector<double> historyEndpointDecoder(const HankelLeftBasis& basis, std::size_t rank,
                                                         std::size_t channels = 1, std::size_t channel = 0,
                                                         double gain = 1.0);

}  // namespace dynsyslab::core
