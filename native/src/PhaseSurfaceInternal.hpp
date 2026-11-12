#pragma once
#include "PhaseSpaceInternal.hpp"
namespace dynsyslab::core::phase_detail::surface {
using Ring = std::vector<PhaseState>;
struct Point {
    PhaseState state;
    double parent{}, residual{}, time{}, foot{};
};
struct Front {
    Ring points, previous;
};
struct Budget {
    std::size_t evaluations{}, maximum{200000};
};
double dot(const PhaseState &, const PhaseState &, const PhaseSpaceOptions &);
double length(const PhaseState &, const PhaseSpaceOptions &);
PhaseState subtract(PhaseState, const PhaseState &);
PhaseState interpolate(const Ring &, double);
PhaseState tangent(const Ring &, double, const PhaseSpaceOptions &);
PhaseState unit(PhaseState, const PhaseSpaceOptions &);
double angle(const PhaseState &, const PhaseState &, const PhaseSpaceOptions &);
Point advance(const ODESystem &, const Front &, double parent, double distance, double sign,
              const PhaseSpaceOptions &, Budget &, const ExecutionContext *);
PhaseSurface grow(const ODESystem &, const PhaseState &center, const Ring &basis, double sign,
                  const PhaseSpaceOptions &, const ExecutionContext *);
void validate(const ODESystem &, PhaseSurface &, const PhaseSurface &refined,
              const PhaseSpaceOptions &, const ExecutionContext *);
Ring tangentBasis(const PhaseFixedPoint &, bool stable, const PhaseSpaceOptions &);
} // namespace dynsyslab::core::phase_detail::surface
