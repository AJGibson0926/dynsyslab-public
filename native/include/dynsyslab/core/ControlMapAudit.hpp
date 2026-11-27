#pragma once
// Audits of an approximate input map B_psi (docs/architecture/KROM_TOOLS.md, control_map_audit):
//   controlMapPathError  |B_approx(x) - B_exact(x)| along the states of a path (for example a
//                        tracking target), as RMS, relative Frobenius and pointwise relative
//                        median / 95th percentile / maximum;
//   controlDesignAudit   persistent-excitation diagnostics of controlled training data: input and
//                        state ranges, the rank of the constant design [u] and the rank and
//                        2-norm condition number of the bilinear design [u_j psi(x), u_j]
//                        (persistency of excitation: Willems et al. 2005, Syst. Control Lett.
//                        54:325-329; bilinear lifted input maps: Korda and Mezic 2018,
//                        Automatica 93:149-160).
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ControlMapPathError {
    std::size_t samples{};
    double rmse{}, relativeFrobenius{}, medianRelative{}, p95Relative{}, maximumRelative{};
};
// states: samples x stateDimension in the model's state coordinates; relative errors divide by
// max(|B_exact(x)|, relativeFloor).
[[nodiscard]] ControlMapPathError controlMapPathError(const KoopmanModel& model, const std::string& approximation,
                                                      const std::vector<double>& states, std::size_t samples,
                                                      double relativeFloor = 0.0);

struct ControlDesignAudit {
    std::size_t rows{}, inputs{}, bilinearColumns{}, constantRank{}, bilinearRank{};
    double bilinearCondition{};  // largest over smallest singular value (infinite when rank deficient)
    double inputMinimum{}, inputMaximum{}, inputMean{}, inputRms{};
    std::vector<double> stateMinimum, stateMaximum, stateP01, stateP99;
};
// states: samples x stateDimension (model coordinates), inputs: samples x inputs (the command held
// from each state); every stride-th row enters the designs, every row the ranges.
[[nodiscard]] ControlDesignAudit controlDesignAudit(const KoopmanModel& model, const std::vector<double>& states,
                                                    const std::vector<double>& inputs, std::size_t samples,
                                                    std::size_t inputCount, std::size_t stride = 10);

}  // namespace dynsyslab::core
