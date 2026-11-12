#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
namespace dynsyslab::core {
void validatePhaseSpaceOptions(const ODESystem &system, const PhaseSpaceOptions &o) {
    const auto n = system.descriptor().stateDimension;
    if (n == 0 || n > 64 || o.lower.size() != n || o.upper.size() != n)
        throw std::invalid_argument(
            "Phase-space analysis requires a box matching 1..64 state components");
    for (std::size_t j = 0; j < n; ++j)
        if (!std::isfinite(o.lower[j]) || !std::isfinite(o.upper[j]) || o.lower[j] >= o.upper[j])
            throw std::invalid_argument("Phase-space bounds must be finite and strictly ordered");
    if (system.descriptor().identifier == "marmottant") {
        const auto &g = system.configuration().dimensionlessGroups;
        const bool buckled = o.lower[0] > 0 && o.upper[0] < std::min(g.at("buckling_radius_ratio"), g.at("rupture_radius_ratio"));
        const bool elastic =
            o.lower[0] > g.at("buckling_radius_ratio") && o.upper[0] < g.at("rupture_radius_ratio");
        const bool ruptured = o.lower[0] > std::max(g.at("breaking_radius_ratio"), g.at("rupture_radius_ratio"));
        if (!buckled && !elastic && !ruptured)
            throw std::invalid_argument("Marmottant sidecar requires a domain strictly within one "
                                        "buckled, elastic or fully ruptured shell regime; radius "
                                        "and velocity alone do not encode rupture history");
    }
    if (system.descriptor().identifier == "reddy_szeri") {
        const auto &c = system.configuration();
        const auto &g = c.dimensionlessGroups;
        if (o.dynamics != "periodic-forcing" && c.parameters.at("translation_enabled") != 0 && g.at("angular_frequency") != 0 &&
            (g.at("pressure_amplitude") != 0 || g.at("pressure_derivative_amplitude") != 0))
            throw std::invalid_argument("Reddy-Szeri phase analysis requires zero effective acoustic translation forcing");
        if (o.lower[0] <= 0 || (c.parameters.at("kappa") != 1 && o.lower[0] <= 1 && o.upper[0] >= 1))
            throw std::invalid_argument("Reddy-Szeri phase analysis requires positive radius and a domain strictly within one gas-exponent branch, or kappa=1");
    }
    if (o.dynamics != "autonomous-zero-input" && o.dynamics != "periodic-forcing")
        throw std::invalid_argument("Unknown phase-space dynamics mode");
    if (o.dynamics == "periodic-forcing") {
        const double period = phaseSpaceForcingPeriod(system, o.forcingPeriod);
        if (!std::isfinite(o.forcingSectionTime) || std::abs(o.forcingSectionTime / period) > 1e9 ||
            o.forcingPeriods < 1 || o.forcingPeriods > 64 || period * o.forcingPeriods > 2 * o.horizon ||
            period * o.forcingPeriods / o.outputStep > 100000 || o.maximumStep > period / 8)
            throw std::invalid_argument("Invalid forcing section, return-period count or sampling; maximum_step must resolve the forcing period with at least eight steps");
    }
    if (o.dynamics == "autonomous-zero-input" && system.configuration().driver.amplitude != 0)
        throw std::invalid_argument(
            "Phase-space analysis requires an autonomous configuration with zero driver");
    for (double value :
         {o.residualTolerance, o.relativeTolerance, o.absoluteTolerance, o.eigenvalueTolerance,
          o.seedOffset, o.horizon, o.outputStep, o.maximumStep, o.manifoldTolerance,
          o.transversalityTolerance})
        if (!std::isfinite(value) || value <= 0)
            throw std::invalid_argument(
                "Phase-space tolerances, offsets and times must be finite and positive");
    if (!std::isfinite(o.transientTime) || o.transientTime < 0 || o.seedOffset >= .1)
        throw std::invalid_argument("Invalid transient time or seed offset");
    if (!std::isfinite(o.sectionMinimumTime) || o.sectionMinimumTime < 0 ||
        o.sectionMinimumTime >= o.horizon || o.sectionDirection < -1 || o.sectionDirection > 1)
        throw std::invalid_argument("Invalid Poincare return time or crossing direction");
    if (o.manifoldIterations < 1 || o.manifoldIterations > 128 ||
        o.maximumManifoldPoints < 32 || o.maximumManifoldPoints > 10000 || o.manifoldTolerance >= .01)
        throw std::invalid_argument("Invalid section-manifold refinement limits");
    for(double value : {o.surface.seedRadius,o.surface.maximumDistance,o.surface.step,o.surface.meshSize,o.surface.chordTolerance,o.surface.tangentTolerance})
        if(!std::isfinite(value)||value<=0)throw std::invalid_argument("Surface geometry tolerances must be finite and positive");
    if(o.surface.initialPoints<8||o.surface.initialPoints>128||o.surface.maximumPoints<64||o.surface.maximumPoints>16384||
       o.surface.maximumPoints<=o.surface.initialPoints||o.surface.maximumFronts<3||o.surface.maximumFronts>240||o.surface.seedRadius>=.05||
       o.surface.maximumDistance<=o.surface.seedRadius||o.surface.maximumDistance>4||o.surface.step>.1||
       o.surface.meshSize>.25||o.surface.chordTolerance>=.02||o.surface.tangentTolerance>=.5)
        throw std::invalid_argument("Surface continuation limits are out of range");
    if (o.gridPerAxis < 2 || o.gridPerAxis > 32 || o.maximumSeeds < 1 || o.maximumSeeds > 4096 ||
        o.maximumNewtonIterations < 1 || o.maximumNewtonIterations > 200 ||
        o.maximumIntegrationSteps < 1 || o.maximumIntegrationSteps > 10000000)
        throw std::invalid_argument("Phase-space search limits are out of range");
    if (o.horizon / o.outputStep > 100000 || o.transientTime / o.outputStep > 100000)
        throw std::invalid_argument("Phase-space sampling exceeds 100000 points per flow");
    for (const auto &seed : o.equilibriumSeeds)
        if (seed.size() != n || !phase_detail::inside(seed, o))
            throw std::invalid_argument("Equilibrium seed must lie inside the declared box");
    if (o.equilibriumSeeds.size() > 4096 || o.periodicSeeds.size() > 128)
        throw std::invalid_argument("Too many explicit phase-space seeds");
    for (const auto &seed : o.periodicSeeds)
        if (seed.state.size() != n || !phase_detail::inside(seed.state, o) ||
            !std::isfinite(seed.period) || (o.dynamics == "periodic-forcing" ? seed.period < 0 : seed.period <= 2 * o.outputStep) ||
            seed.period > 2 * o.horizon)
            throw std::invalid_argument("Invalid periodic seed state or period");
    if (o.sectionAxis >= n || !std::isfinite(o.sectionValue))
        throw std::invalid_argument("Invalid Poincare section");
    if ((!o.sectionNormal.empty() || !o.sectionReference.empty()) &&
        (o.sectionNormal.size() != n || o.sectionReference.size() != n))
        throw std::invalid_argument(
            "Affine Poincare section normal and reference must match the state dimension");
    for (double value : o.sectionNormal)
        if (!std::isfinite(value))
            throw std::invalid_argument("Poincare section normal must be finite");
    for (double value : o.sectionReference)
        if (!std::isfinite(value))
            throw std::invalid_argument("Poincare section reference must be finite");
    (void)phase_detail::sectionGeometry(o, n);
    for (const auto &seed : o.sectionSeeds)
        if (seed.size() != n || !phase_detail::inside(seed, o))
            throw std::invalid_argument("Poincare seed must lie inside the declared box");
    if (o.sectionSeeds.size() > 256)
        throw std::invalid_argument("Too many explicit Poincare section seeds");
    if (o.bifurcation.enabled) {
        const auto &b = o.bifurcation;
        if (b.parameter.empty() || (b.branch != "equilibrium" && b.branch != "periodic-orbit") ||
            !(b.minimum < b.maximum) || b.initial < b.minimum || b.initial > b.maximum ||
            !std::isfinite(b.step) || b.step == 0 || b.minimumStep <= 0 ||
            b.maximumStep < b.minimumStep || b.maximumPoints < 2 || b.maximumPoints > 4096 ||
            b.maximumCorrectorIterations < 2 || b.maximumCorrectorIterations > 100)
            throw std::invalid_argument("Invalid bifurcation continuation options");
        if (!b.seedState.empty() && b.seedState.size() != n)
            throw std::invalid_argument("Bifurcation seed state must match the state dimension");
        if (b.branch == "periodic-orbit" && o.dynamics == "autonomous-zero-input" &&
            (!std::isfinite(b.seedPeriod) || b.seedPeriod <= 2 * o.outputStep))
            throw std::invalid_argument(
                "Autonomous periodic-orbit continuation requires a positive seed period");
    }
    if (o.vectorField) {
        if (o.vectorAxes.size() != 2 || o.vectorAxes[0] >= n || o.vectorAxes[1] >= n ||
            o.vectorAxes[0] == o.vectorAxes[1])
            throw std::invalid_argument("Vector field requires two distinct state axes");
        if (o.vectorGridPerAxis < 3 || o.vectorGridPerAxis > 25)
            throw std::invalid_argument("Vector grid must have 3..25 points per axis");
        if (n > 2 && (o.sliceState.size() != n || !phase_detail::inside(o.sliceState, o)))
            throw std::invalid_argument(
                "A projected vector field requires an explicit full-state slice");
    }
}
PhaseSpaceResult analyzePhaseSpace(const ODESystem &system, const PhaseSpaceOptions &o,
                                   const ExecutionContext *execution) {
    validatePhaseSpaceOptions(system, o);
    if (execution) {
        execution->checkpoint();
        execution->requireMemory(100000 * (system.descriptor().stateDimension + 1) *
                                     sizeof(double) * 4,
                                 "phase-space flow workspaces");
    }
    PhaseSpaceResult result;
    result.options = o;
    const bool periodic = o.dynamics == "periodic-forcing";
    if (periodic) result.options.forcingPeriod = phaseSpaceForcingPeriod(system, o.forcingPeriod);
    result.equationsFingerprint = equationsFingerprint(system.configuration());
    result.coordinateSystem = system.descriptor().preferredCoordinates;
    result.diagnostics.push_back(periodic
        ? "Periodically forced dynamics at a fixed forcing phase; finite shooting searches do not establish completeness"
        : "Autonomous zero-input dynamics; finite multistart searches do not establish completeness");
    if (periodic) {
        if (o.fixedPoints || o.periodicOrbits || o.manifolds)
            phase_detail::findForcedPeriodicOrbits(system, result, execution);
        if (o.manifolds) phase_detail::traceForcedManifolds(system, result, execution);
        result.diagnostics.push_back("Section manifolds are curves of the stroboscopic map, not whole flow manifolds; vector arrows show the instantaneous field at the recorded forcing section time");
    }
    if (!periodic && (o.fixedPoints || o.manifolds)) {
        phase_detail::applyEquilibriumGuides(system, result);
        phase_detail::findFixedPoints(system, result, execution);
    }
    if (!periodic && o.manifolds)
        phase_detail::traceManifolds(system, result, execution);
    if (!periodic && o.periodicOrbits && system.descriptor().stateDimension > 1)
        phase_detail::findPeriodicOrbits(system, result, execution);
    if (!periodic && o.poincareMap && !o.sectionSeeds.empty())
        phase_detail::findPoincareReturns(system, result, execution);
    if (o.vectorField) {
        auto &field = result.vectorField;
        field.axes = o.vectorAxes;
        field.sliceState = o.sliceState;
        const auto n = system.descriptor().stateDimension, a = o.vectorAxes[0], b = o.vectorAxes[1];
        for (std::size_t row = 0; row < o.vectorGridPerAxis; ++row)
            for (std::size_t col = 0; col < o.vectorGridPerAxis; ++col) {
                if (execution)
                    execution->checkpoint();
                auto x = o.sliceState;
                if (x.empty())
                    x.resize(n);
                x[a] = o.lower[a] + (o.upper[a] - o.lower[a]) * (col + .5) / o.vectorGridPerAxis;
                x[b] = o.lower[b] + (o.upper[b] - o.lower[b]) * (row + .5) / o.vectorGridPerAxis;
                try {
                    auto f = phase_detail::fieldAt(system, x, periodic ? o.forcingSectionTime : 0, periodic);
                    field.points.insert(field.points.end(), x.begin(), x.end());
                    field.vectors.insert(field.vectors.end(), f.begin(), f.end());
                } catch (const std::exception &) {
                } // Singular/invalid cells remain absent, never zero-filled.
            }
        if (n > 2)
            result.diagnostics.push_back("Vector arrows show the declared full-state slice, not a "
                                         "closed autonomous model of its projection");
    }
    if (!periodic) phase_detail::checkKnownInvariants(system, result);
    phase_detail::findBifurcationBranches(system, result, execution);
    result.options = o;
    if (periodic) result.options.forcingPeriod = phaseSpaceForcingPeriod(system, o.forcingPeriod); // Keep the user query distinct from verified system-generated seeds.
    if (execution)
        execution->report({"phase/complete", "Phase-space analysis complete", 1, 1, 1});
    return result;
}
} // namespace dynsyslab::core
