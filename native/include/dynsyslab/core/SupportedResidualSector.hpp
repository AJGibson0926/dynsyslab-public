#pragma once
// A residual sector ||r|| <= c_x ||z|| + c_u ||u|| whose two coefficients are each
// supported by data, for use in a certificate.
//
// The joint LP in ResidualSector.hpp returns a vertex, so it can place the whole residual
// on one channel. If the records never have u = 0, c_x = 0 is then feasible even when the
// autonomous residual is not zero. This estimator identifies the channels in sequence:
//
//   c_x = max over zero-input holds of ||r|| / ||z||
//   c_u = the larger of
//         max over forced pairs of max(||r|| - c_x ||z||, 0) / ||u||    (envelope), and
//         max over pairs with ||u|| >= 0.1 max ||u|| of
//             max(||r|| - largest hold residual in the same ||z|| bin, 0) / ||u||
//                                                                       (identification)
//
// It is refused (supported = false) when zero-input holds do not span the sampled ||z||
// range or the stated radius, when forced pairs do not reach the stated input bound, or
// when a coefficient is zero although that channel's residual is not. A one-sided sector
// needs an explicit declaration with a reason.
//
// The error structure is that of SafEDMD (Strasser, Worthmann et al., arXiv:2402.03145).
// The two-stage estimate is a DynSysLab construction, not a published method. It is still
// a data envelope: it gives no probability and no bound on unseen states.
#include <dynsyslab/core/ResidualSector.hpp>

#include <cstddef>
#include <optional>
#include <string>

namespace dynsyslab::core {

struct OneSidedSectorDeclaration {
    std::string channel;  // "state" declares c_x = 0; "input" declares c_u = 0
    std::string reason;   // required, recorded with the result
};

struct SupportedSectorOptions {
    double inflation{1.0};                // >= 1, multiplies both coefficients
    double zeroInputTolerance{1e-9};      // ||u|| at or below this * max ||u|| is a zero-input hold
    double zeroResidualTolerance{1e-10};  // relative to the largest ||z_next||
    std::size_t minimumZeroInputPairs{8}, minimumForcedPairs{8};
    std::size_t radialBins{4};            // bins of ||z|| over (0, largest sampled ||z||]
    std::size_t minimumPairsPerBin{2};    // zero-input pairs needed in every occupied bin
    double stateRadius{0.0};              // > 0: zero-input holds must reach coverage * radius
    double inputBound{0.0};               // > 0: forced pairs must reach coverage * bound
    double minimumCoverage{0.5};
    std::optional<OneSidedSectorDeclaration> oneSided;
};

struct SupportedResidualSector {
    ResidualSectorBound bound;     // bound.feasible is true only when supported
    bool supported{false};
    std::string reason;            // why it is, or is not, admissible
    std::string oneSided;          // "", "declared_state", "declared_input" or "zero_residual"
    std::string declarationReason;
    std::size_t zeroInputPairs{}, forcedPairs{};
    double largestZeroInputStateNorm{}, largestStateNorm{}, largestInputNorm{};
    double largestZeroInputResidual{}, largestForcedResidual{};
};

[[nodiscard]] SupportedResidualSector estimateSupportedResidualSector(const BilinearKrom& model,
    const BilinearDataset& data, const SupportedSectorOptions& options = {});

}  // namespace dynsyslab::core
