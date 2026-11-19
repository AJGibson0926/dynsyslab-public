// PBH cross-check of the Kalman structure report (dynsyslab_linear_analysis).
// The Kalman decomposition itself is a dense kernel in dynsyslab_numerics.
#include <dynsyslab/core/SystemStructure.hpp>
#include <dynsyslab/core/Pbh.hpp>

#include <algorithm>

namespace dynsyslab::core {

StructureReport analyzeStructure(const StateSpaceModel& model, double tolerance) {
    auto out = kalmanStructureReport(model, tolerance);
    const bool discrete = model.domain == TimeDomain::discrete;
    const auto pc = pbhControllability(model.A, model.B, model.states, model.inputs, discrete, tolerance);
    const auto po = pbhObservability(model.A, model.C, model.states, model.outputs, discrete, tolerance);
    out.pbhControllable = pc.fullRank; out.pbhObservable = po.fullRank;
    out.pbhStabilizable = pc.stableHiddenModes; out.pbhDetectable = po.stableHiddenModes;
    out.pbhControllabilityMargin = pc.relativeMargins.empty() ? 0.0
        : *std::min_element(pc.relativeMargins.begin(), pc.relativeMargins.end());
    out.pbhObservabilityMargin = po.relativeMargins.empty() ? 0.0
        : *std::min_element(po.relativeMargins.begin(), po.relativeMargins.end());
    out.pbhAgrees = pc.fullRank == out.controllable && po.fullRank == out.observable
        && pc.stableHiddenModes == out.stabilizable && po.stableHiddenModes == out.detectable;
    return out;
}

}  // namespace dynsyslab::core
