// Two-stage residual sector. Error structure from SafEDMD (Strasser, Worthmann et al.,
// arXiv:2402.03145); the sequential estimate is a DynSysLab construction. See the header.
#include <dynsyslab/core/SupportedResidualSector.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace dynsyslab::core {
namespace {

double rowNorm(const std::vector<double>& values, std::size_t row, std::size_t width) {
    double sum = 0.0;
    for (std::size_t j = 0; j < width; ++j) sum += values[row * width + j] * values[row * width + j];
    return std::sqrt(sum);
}

void validate(const SupportedSectorOptions& o) {
    if (!(o.inflation >= 1.0) || !std::isfinite(o.inflation)) throw std::invalid_argument("Sector inflation must be at least 1");
    if (!(o.zeroInputTolerance >= 0.0 && o.zeroInputTolerance < 1.0)) throw std::invalid_argument("zero-input tolerance must lie in [0, 1)");
    if (!(o.zeroResidualTolerance >= 0.0) || !std::isfinite(o.zeroResidualTolerance)) throw std::invalid_argument("zero-residual tolerance must be non-negative");
    if (!o.radialBins || !o.minimumPairsPerBin || !o.minimumZeroInputPairs || !o.minimumForcedPairs)
        throw std::invalid_argument("Sector support counts must be positive");
    if (!(o.stateRadius >= 0.0) || !(o.inputBound >= 0.0) || !(o.minimumCoverage > 0.0 && o.minimumCoverage <= 1.0))
        throw std::invalid_argument("Sector coverage settings are invalid");
    if (o.oneSided) {
        if (o.oneSided->channel != "state" && o.oneSided->channel != "input")
            throw std::invalid_argument("A one-sided sector declares channel state or input");
        if (o.oneSided->reason.empty()) throw std::invalid_argument("A one-sided sector needs a stated reason");
    }
}

}  // namespace

SupportedResidualSector estimateSupportedResidualSector(const BilinearKrom& model, const BilinearDataset& data,
    const SupportedSectorOptions& o) {
    validate(o);
    const auto n = model.stateDimension, m = model.inputDimension;
    const auto residual = bilinearResiduals(model, data);
    const auto N = data.samples;
    std::vector<double> zn(N), un(N), rn(N);
    double nextScale = 0.0, largestResidual = 0.0;
    SupportedResidualSector out;
    for (std::size_t k = 0; k < N; ++k) {
        zn[k] = rowNorm(data.z, k, n); un[k] = rowNorm(data.u, k, m); rn[k] = rowNorm(residual, k, n);
        nextScale = std::max(nextScale, rowNorm(data.zNext, k, n));
        largestResidual = std::max(largestResidual, rn[k]);
        out.largestStateNorm = std::max(out.largestStateNorm, zn[k]);
        out.largestInputNorm = std::max(out.largestInputNorm, un[k]);
    }
    const double zeroResidual = o.zeroResidualTolerance * (nextScale > 0.0 ? nextScale : 1.0);
    const double zeroInput = o.zeroInputTolerance * out.largestInputNorm;
    const double zeroState = 1e-13 * out.largestStateNorm;
    std::vector<char> hold(N);
    for (std::size_t k = 0; k < N; ++k) {
        hold[k] = un[k] <= zeroInput;
        if (hold[k]) {
            ++out.zeroInputPairs;
            out.largestZeroInputStateNorm = std::max(out.largestZeroInputStateNorm, zn[k]);
            out.largestZeroInputResidual = std::max(out.largestZeroInputResidual, rn[k]);
        } else {
            ++out.forcedPairs;
            out.largestForcedResidual = std::max(out.largestForcedResidual, rn[k]);
        }
    }
    out.bound.samples = N;
    out.bound.inflation = o.inflation;
    out.bound.maximumResidual = largestResidual;
    out.bound.guarantee = "two-stage data envelope on the supplied samples; no statement about unseen states";
    const auto refuse = [&](std::string why) {
        out.supported = false; out.bound.feasible = false; out.bound.cState = out.bound.cInput = 0.0;
        out.reason = std::move(why); out.bound.message = out.reason;
        return out;
    };
    const auto accept = [&](double cx, double cu, std::string why) {
        out.supported = true; out.bound.feasible = true;
        out.bound.cState = o.inflation * cx; out.bound.cInput = o.inflation * cu;
        out.reason = std::move(why); out.bound.message = out.reason;
        return out;
    };

    if (o.oneSided) {
        out.declarationReason = o.oneSided->reason;
        double c = 0.0;
        if (o.oneSided->channel == "state") {
            out.oneSided = "declared_state";
            if (out.forcedPairs < o.minimumForcedPairs) return refuse("declared c_x = 0 needs forced pairs to estimate c_u");
            for (std::size_t k = 0; k < N; ++k) {
                if (un[k] > zeroInput) c = std::max(c, rn[k] / un[k]);
                else if (rn[k] > zeroResidual) return refuse("declared c_x = 0 is contradicted: a zero-input hold has a nonzero residual");
            }
            return accept(0.0, c, "one-sided sector declared by the caller (c_x = 0)");
        }
        out.oneSided = "declared_input";
        for (std::size_t k = 0; k < N; ++k) {
            if (zn[k] > zeroState) c = std::max(c, rn[k] / zn[k]);
            else if (rn[k] > zeroResidual) return refuse("declared c_u = 0 is contradicted: a pair at z = 0 has a nonzero residual");
        }
        return accept(c, 0.0, "one-sided sector declared by the caller (c_u = 0)");
    }

    if (out.zeroInputPairs < o.minimumZeroInputPairs)
        return refuse("c_x is unsupported: " + std::to_string(out.zeroInputPairs) + " zero-input holds, need "
                      + std::to_string(o.minimumZeroInputPairs));
    if (out.forcedPairs < o.minimumForcedPairs)
        return refuse("c_u is unsupported: " + std::to_string(out.forcedPairs) + " forced pairs, need "
                      + std::to_string(o.minimumForcedPairs));
    if (o.stateRadius > 0.0 && out.largestZeroInputStateNorm < o.minimumCoverage * o.stateRadius)
        return refuse("c_x is unsupported: zero-input holds reach ||z|| = " + std::to_string(out.largestZeroInputStateNorm)
                      + ", below the required fraction of state_radius");
    if (o.inputBound > 0.0 && out.largestInputNorm < o.minimumCoverage * o.inputBound)
        return refuse("c_u is unsupported: forced pairs reach ||u|| = " + std::to_string(out.largestInputNorm)
                      + ", below the required fraction of the input bound");
    std::vector<std::size_t> all(o.radialBins, 0), holds(o.radialBins, 0);
    for (std::size_t k = 0; k < N; ++k) {
        if (!(out.largestStateNorm > 0.0)) break;
        const auto bin = std::min(o.radialBins - 1, static_cast<std::size_t>(zn[k] / out.largestStateNorm * double(o.radialBins)));
        ++all[bin];
        if (hold[k]) ++holds[bin];
    }
    for (std::size_t b = 0; b < o.radialBins; ++b)
        if (all[b] && holds[b] < o.minimumPairsPerBin)
            return refuse("c_x is unsupported: too few zero-input holds in ||z|| bin " + std::to_string(b + 1) + " of "
                          + std::to_string(o.radialBins));

    double cx = 0.0;
    if (out.largestZeroInputResidual > zeroResidual)
        for (std::size_t k = 0; k < N; ++k) {
            if (!hold[k]) continue;
            if (zn[k] > zeroState) cx = std::max(cx, rn[k] / zn[k]);
            else if (rn[k] > zeroResidual) {
                auto refused = refuse("a nonzero residual occurs at z = 0, u = 0: the target is not an equilibrium of the data");
                return refused;
            }
        }
    // Envelope term: every forced pair must satisfy the final bound. Identification term:
    // what a substantial input adds above the largest zero-input residual at the same radius.
    // The global c_x alone can hide the input channel, since it is set by the largest radius.
    double cu = 0.0;
    if (out.largestForcedResidual > zeroResidual) {
        std::vector<double> holdResidual(o.radialBins, 0.0);
        const auto binOf = [&](std::size_t k) {
            return std::min(o.radialBins - 1, static_cast<std::size_t>(zn[k] / out.largestStateNorm * double(o.radialBins)));
        };
        for (std::size_t k = 0; k < N; ++k) if (hold[k]) holdResidual[binOf(k)] = std::max(holdResidual[binOf(k)], rn[k]);
        for (std::size_t k = 0; k < N; ++k) {
            if (hold[k]) continue;
            cu = std::max(cu, std::max(rn[k] - cx * zn[k], 0.0) / un[k]);
            if (un[k] >= 0.1 * out.largestInputNorm) cu = std::max(cu, std::max(rn[k] - holdResidual[binOf(k)], 0.0) / un[k]);
        }
    }
    if (cu == 0.0 && out.largestForcedResidual > zeroResidual)
        return refuse("c_u is unsupported: the forced residual is nonzero but never exceeds the c_x envelope, "
                      "so the input channel is not identified");
    if (cx == 0.0 && cu == 0.0) { out.oneSided = "zero_residual"; return accept(0.0, 0.0, "residual is zero to tolerance on both channels"); }
    if (cx == 0.0 || cu == 0.0) out.oneSided = "zero_residual";
    return accept(cx, cu, cx == 0.0 ? "c_x = 0: the zero-input residual is zero to tolerance"
        : cu == 0.0 ? "c_u = 0: the forced residual is zero to tolerance"
                    : "both channels identified: c_x from zero-input holds, c_u from the remaining forced residual");
}

}  // namespace dynsyslab::core
