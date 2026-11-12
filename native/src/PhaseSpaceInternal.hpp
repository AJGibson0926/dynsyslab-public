#pragma once
#include <dynsyslab/core/PhaseSpace.hpp>
namespace dynsyslab::core::phase_detail {
void applyEquilibriumGuides(const ODESystem &, PhaseSpaceResult &);
void checkKnownInvariants(const ODESystem &, PhaseSpaceResult &);
struct Flow {
    PhaseState state, monodromy;
    PhaseCurve curve;
    bool complete{};
    std::string reason;
};
struct SectionGeometry {
    PhaseState scaledNormal, physicalNormal, reference, tangentBasis;
};
SectionGeometry sectionGeometry(const PhaseSpaceOptions &, std::size_t dimension);
double sectionValue(const PhaseState &, const SectionGeometry &);
PhasePoincareMap poincareFromTransition(const ODESystem &, const PhaseState &, const Flow &,
                                        const PhaseSpaceOptions &, const SectionGeometry &);
double norm(const PhaseState &x, const PhaseSpaceOptions &o);
bool inside(const PhaseState &x, const PhaseSpaceOptions &o);
PhaseState field(const ODESystem &system, const PhaseState &x);
PhaseState fieldAt(const ODESystem &, const PhaseState &, double time, bool periodic);
Flow flow(const ODESystem &system, const PhaseState &initial, double duration,
          const PhaseSpaceOptions &options, bool variational, bool stopAtDomain,
          const ExecutionContext *execution, double startTime = 0);
SectionGeometry sectionGeometry(const PhaseSpaceOptions &options, std::size_t dimension);
double sectionValue(const PhaseState &state, const SectionGeometry &section);
PhasePoincareMap poincareFromTransition(const ODESystem &system, const PhaseState &initial,
                                        const Flow &transition, const PhaseSpaceOptions &options,
                                        const SectionGeometry &section);
void findPoincareReturns(const ODESystem &, PhaseSpaceResult &, const ExecutionContext *);
void findBifurcationBranches(const ODESystem &, PhaseSpaceResult &, const ExecutionContext *);
void detectBifurcationEvents(const PhaseSystemFactory &, const PhaseSpaceOptions &,
                             const PhaseBifurcationOptions &, PhaseBifurcationBranch &);
void findForcedPeriodicOrbits(const ODESystem &, PhaseSpaceResult &, const ExecutionContext *);
void traceForcedManifolds(const ODESystem &, PhaseSpaceResult &, const ExecutionContext *);
std::vector<PhaseState> seeds(const PhaseSpaceOptions &o, std::size_t dimension);
void findFixedPoints(const ODESystem &system, PhaseSpaceResult &result,
                     const ExecutionContext *execution);
void traceSurface(const ODESystem &, PhaseSpaceResult &, std::size_t, bool, const ExecutionContext *);
void traceManifolds(const ODESystem &system, PhaseSpaceResult &result,
                    const ExecutionContext *execution);
void findPeriodicOrbits(const ODESystem &system, PhaseSpaceResult &result,
                        const ExecutionContext *execution);
} // namespace dynsyslab::core::phase_detail
