#pragma once
// LTI realization of a Koopman reduced-order model in its real split
// coordinates: A = generator, B = a named control map, C = stacked named
// reconstructions (measurements). Reconstruction intercepts are affine output
// offsets; they do not change any linear-systems quantity and are omitted.
#include <dynsyslab/core/StateSpace.hpp>

#include <string>
#include <vector>

namespace dynsyslab::core {
class KoopmanModel;

struct KoopmanStateSpaceOptions {
    std::vector<std::string> measurements; // reconstruction names; stacked as C rows
    std::string controlMap;                // empty: no input
    // State-dependent maps ("exact", bilinear) are linearized at an explicit
    // operating point: a physical state, or the real split lifted values.
    std::vector<double> operatingState;
    std::vector<double> operatingLifted;
    TimeDomain domain{TimeDomain::continuous};
    // Discrete: a stored discrete input map of this name is used directly with
    // its own sample interval; otherwise the continuous realization is sampled
    // by zero-order hold at sampleTime.
    double sampleTime{};
};

[[nodiscard]] StateSpaceModel koopmanStateSpace(const KoopmanModel& model,
    const KoopmanStateSpaceOptions& options);

}  // namespace dynsyslab::core
