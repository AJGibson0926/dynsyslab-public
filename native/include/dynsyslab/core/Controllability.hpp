#pragma once
#include <dynsyslab/core/Observability.hpp>
namespace dynsyslab::core {
using ControllabilityMethod = ObservabilityMethod;
using ControllabilityOptions = ObservabilityOptions;
struct ControllabilityResult {
    std::vector<double> matrix; // state x (input * blockCount)
    std::vector<double> gramian;
    std::vector<double> singularValues;
    std::vector<double> gramianEigenvalues;
    std::vector<double> leftDirections; // state x state: most controllable first
    std::vector<double> rightDirections; // input history x min(rows,columns)
    ControllabilityMethod method{};
    std::size_t stateDimension{}, inputDimension{}, blockCount{}, rank{};
    double rankTolerance{}, horizon{}, conditionNumber{}, resolvedConditionNumber{}, effectiveRank{};
};
[[nodiscard]] ControllabilityResult analyzeControllability(
    const std::vector<double>& A, const std::vector<double>& B,
    std::size_t states, std::size_t inputs,
    const ControllabilityOptions& options = {});
// Constant B only. State-dependent maps require an explicit operating point and
// the matrix overload; they are not silently treated as LTI input matrices.
[[nodiscard]] ControllabilityResult analyzeControllability(
    const KoopmanModel& model, const std::string& controlMap,
    const ControllabilityOptions& options = {});
} // namespace dynsyslab::core
