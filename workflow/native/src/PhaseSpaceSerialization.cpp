#include "PhaseSpaceSerialization.hpp"
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
json::Value array(const std::vector<double> &values) {
    json::Value::Array result;
    for (double q : values)
        result.emplace_back(std::isfinite(q) ? json::Value(q) : json::Value());
    return result;
}
json::Value complexArray(const std::vector<std::complex<double>> &values) {
    json::Value::Array result;
    for (auto q : values)
        result.push_back(array({q.real(), q.imag()}));
    return result;
}
std::vector<double> numbers(const json::Value &value) {
    std::vector<double> out;
    for (const auto &q : value.array()) {
        if (!std::holds_alternative<double>(q.storage()) || !std::isfinite(q.number()))
            throw std::invalid_argument("Phase-space arrays require finite numbers");
        out.push_back(q.number());
    }
    return out;
}
std::size_t count(const json::Value &v, const char *key, std::size_t fallback) {
    if (v.at(key).isNull())
        return fallback;
    double q = v.at(key).number(-1);
    if (q < 0 || q > 10000000 || q != std::floor(q))
        throw std::invalid_argument(std::string("Invalid count: ") + key);
    return std::size_t(q);
}
json::Value curveJson(const dynsyslab::core::PhaseCurve &curve) {
    return json::Value::Object{{"kind", curve.kind},
                               {"base_kind", curve.baseKind},
                               {"parameterization", curve.time.empty() ? "forcing-section curve" : "time"},
                               {"status", curve.status},
                               {"reason", curve.reason},
                               {"base_object", double(curve.baseObject)},
                               {"time", array(curve.time)},
                               {"values", array(curve.values)},
                               {"seed_offset", curve.seedOffset},
                               {"refinement_error", std::isfinite(curve.refinementError)
                                                        ? json::Value(curve.refinementError)
                                                        : json::Value()},
                               {"tangent_defect", curve.tangentDefect}};
}
json::Value poincareJson(const dynsyslab::core::PhasePoincareMap &map) {
    return json::Value::Object{
        {"kind", map.kind}, {"status", map.status}, {"reason", map.reason},
        {"initial_state", array(map.initialState)}, {"section_state", array(map.sectionState)},
        {"return_state", array(map.returnState)}, {"normal", array(map.normal)},
        {"tangent_basis", array(map.tangentBasis)},
        {"return_derivative", array(map.returnDerivative)},
        {"reduced_derivative", array(map.reducedDerivative)},
        {"return_time_derivative", array(map.returnTimeDerivative)},
        {"multipliers", complexArray(map.multipliers)}, {"return_time", map.returnTime},
        {"fixed_point_residual", map.fixedPointResidual},
        {"transversality", map.transversality},
        {"transversality_angle", map.transversalityAngle},
        {"endpoint_refinement_error", map.endpointRefinementError},
        {"derivative_refinement_error", map.derivativeRefinementError},
        {"direction", double(map.direction)}};
}
} // namespace
json::Value phaseSpaceOptionsJson(const dynsyslab::core::PhaseSpaceOptions &o) {
    json::Value::Array seeds, cycles, sectionSeeds, axes;
    for (const auto &x : o.equilibriumSeeds)
        seeds.push_back(array(x));
    for (const auto &s : o.periodicSeeds)
        cycles.emplace_back(json::Value::Object{{"state", array(s.state)}, {"period", s.period}});
    for (const auto &x : o.sectionSeeds)
        sectionSeeds.push_back(array(x));
    for (auto axis : o.vectorAxes)
        axes.emplace_back(double(axis));
    const auto &m=o.surface;
    json::Value surface=json::Value::Object{{"enabled",m.enabled},{"initial_points",double(m.initialPoints)},
        {"maximum_points",double(m.maximumPoints)},{"maximum_fronts",double(m.maximumFronts)},
        {"seed_radius",m.seedRadius},{"maximum_distance",m.maximumDistance},{"step",m.step},
        {"mesh_size",m.meshSize},{"chord_tolerance",m.chordTolerance},{"tangent_tolerance",m.tangentTolerance}};
    const auto &b=o.bifurcation;
    json::Value bifurcation=json::Value::Object{{"enabled",b.enabled},{"branch",b.branch},
        {"parameter",b.parameter},{"minimum",b.minimum},{"maximum",b.maximum},
        {"initial",b.initial},{"step",b.step},{"minimum_step",b.minimumStep},
        {"maximum_step",b.maximumStep},{"maximum_points",double(b.maximumPoints)},
        {"maximum_corrector_iterations",double(b.maximumCorrectorIterations)},
        {"seed_state",array(b.seedState)},{"seed_period",b.seedPeriod}};
    return json::Value::Object{{"surface", surface},{"bifurcation",bifurcation},
                               {"dynamics", o.dynamics},
                               {"forcing_period", o.forcingPeriod},
                               {"forcing_section_time", o.forcingSectionTime},
                               {"forcing_periods", double(o.forcingPeriods)},
                               {"manifold_iterations", double(o.manifoldIterations)},
                               {"maximum_manifold_points", double(o.maximumManifoldPoints)},
                               {"manifold_tolerance", o.manifoldTolerance},
                               {"lower", array(o.lower)},
                               {"upper", array(o.upper)},
                               {"equilibrium_seeds", seeds},
                               {"periodic_seeds", cycles},
                               {"section_seeds", sectionSeeds},
                               {"fixed_points", o.fixedPoints},
                               {"manifolds", o.manifolds},
                               {"periodic_orbits", o.periodicOrbits},
                               {"poincare_map", o.poincareMap},
                               {"vector_field", o.vectorField},
                               {"grid_per_axis", double(o.gridPerAxis)},
                               {"maximum_seeds", double(o.maximumSeeds)},
                               {"vector_grid_per_axis", double(o.vectorGridPerAxis)},
                               {"maximum_newton_iterations", double(o.maximumNewtonIterations)},
                               {"maximum_integration_steps", double(o.maximumIntegrationSteps)},
                               {"residual_tolerance", o.residualTolerance},
                               {"relative_tolerance", o.relativeTolerance},
                               {"absolute_tolerance", o.absoluteTolerance},
                               {"eigenvalue_tolerance", o.eigenvalueTolerance},
                               {"seed_offset", o.seedOffset},
                               {"horizon", o.horizon},
                               {"transient_time", o.transientTime},
                               {"output_step", o.outputStep},
                               {"maximum_step", o.maximumStep},
                               {"section_axis", double(o.sectionAxis)},
                               {"section_value", o.sectionValue},
                               {"section_normal", array(o.sectionNormal)},
                               {"section_reference", array(o.sectionReference)},
                               {"section_direction", double(o.sectionDirection)},
                               {"section_minimum_time", o.sectionMinimumTime},
                               {"transversality_tolerance", o.transversalityTolerance},
                               {"vector_axes", axes},
                               {"slice_state", array(o.sliceState)}};
}
dynsyslab::core::PhaseSpaceOptions phaseSpaceOptionsFromJson(const json::Value &v,
                                                             std::size_t dimension) {
    dynsyslab::core::PhaseSpaceOptions o;
    const auto &m=v.at("surface");
    o.surface.enabled=m.at("enabled").boolean(o.surface.enabled);
    o.surface.initialPoints=count(m,"initial_points",o.surface.initialPoints);
    o.surface.maximumPoints=count(m,"maximum_points",o.surface.maximumPoints);
    o.surface.maximumFronts=count(m,"maximum_fronts",o.surface.maximumFronts);
    o.surface.seedRadius=m.at("seed_radius").number(o.surface.seedRadius);
    o.surface.maximumDistance=m.at("maximum_distance").number(o.surface.maximumDistance);
    o.surface.step=m.at("step").number(o.surface.step);
    o.surface.meshSize=m.at("mesh_size").number(o.surface.meshSize);
    o.surface.chordTolerance=m.at("chord_tolerance").number(o.surface.chordTolerance);
    o.surface.tangentTolerance=m.at("tangent_tolerance").number(o.surface.tangentTolerance);
    const auto &b=v.at("bifurcation");
    o.bifurcation.enabled=b.at("enabled").boolean(o.bifurcation.enabled);
    o.bifurcation.branch=b.at("branch").string(o.bifurcation.branch);
    o.bifurcation.parameter=b.at("parameter").string(o.bifurcation.parameter);
    o.bifurcation.minimum=b.at("minimum").number(o.bifurcation.minimum);
    o.bifurcation.maximum=b.at("maximum").number(o.bifurcation.maximum);
    o.bifurcation.initial=b.at("initial").number(o.bifurcation.initial);
    o.bifurcation.step=b.at("step").number(o.bifurcation.step);
    o.bifurcation.minimumStep=b.at("minimum_step").number(o.bifurcation.minimumStep);
    o.bifurcation.maximumStep=b.at("maximum_step").number(o.bifurcation.maximumStep);
    o.bifurcation.maximumPoints=count(b,"maximum_points",o.bifurcation.maximumPoints);
    o.bifurcation.maximumCorrectorIterations=count(
        b,"maximum_corrector_iterations",o.bifurcation.maximumCorrectorIterations);
    o.bifurcation.seedState=numbers(b.at("seed_state"));
    o.bifurcation.seedPeriod=b.at("seed_period").number(o.bifurcation.seedPeriod);
    o.dynamics = v.at("dynamics").string(o.dynamics);
    o.forcingPeriod = v.at("forcing_period").number(o.forcingPeriod);
    o.forcingSectionTime = v.at("forcing_section_time").number(o.forcingSectionTime);
    o.forcingPeriods = count(v, "forcing_periods", o.forcingPeriods);
    o.manifoldIterations = count(v, "manifold_iterations", o.manifoldIterations);
    o.maximumManifoldPoints = count(v, "maximum_manifold_points", o.maximumManifoldPoints);
    o.manifoldTolerance = v.at("manifold_tolerance").number(o.manifoldTolerance);
    o.lower = numbers(v.at("lower"));
    o.upper = numbers(v.at("upper"));
    o.sliceState = numbers(v.at("slice_state"));
    for (const auto &q : v.at("equilibrium_seeds").array())
        o.equilibriumSeeds.push_back(numbers(q));
    for (const auto &q : v.at("periodic_seeds").array())
        o.periodicSeeds.push_back({numbers(q.at("state")), q.at("period").number()});
    for (const auto &q : v.at("section_seeds").array())
        o.sectionSeeds.push_back(numbers(q));
    o.fixedPoints = v.at("fixed_points").boolean(o.fixedPoints);
    o.manifolds = v.at("manifolds").boolean(o.manifolds);
    o.periodicOrbits = v.at("periodic_orbits").boolean(o.periodicOrbits);
    o.poincareMap = v.at("poincare_map").boolean(o.poincareMap);
    o.vectorField = v.at("vector_field").boolean(dimension > 1);
    o.gridPerAxis = count(v, "grid_per_axis", o.gridPerAxis);
    o.maximumSeeds = count(v, "maximum_seeds", o.maximumSeeds);
    o.vectorGridPerAxis = count(v, "vector_grid_per_axis", o.vectorGridPerAxis);
    o.maximumNewtonIterations = count(v, "maximum_newton_iterations", o.maximumNewtonIterations);
    o.maximumIntegrationSteps = count(v, "maximum_integration_steps", o.maximumIntegrationSteps);
    o.sectionAxis = count(v, "section_axis", o.sectionAxis);
    o.residualTolerance = v.at("residual_tolerance").number(o.residualTolerance);
    o.relativeTolerance = v.at("relative_tolerance").number(o.relativeTolerance);
    o.absoluteTolerance = v.at("absolute_tolerance").number(o.absoluteTolerance);
    o.eigenvalueTolerance = v.at("eigenvalue_tolerance").number(o.eigenvalueTolerance);
    o.seedOffset = v.at("seed_offset").number(o.seedOffset);
    o.horizon = v.at("horizon").number(o.horizon);
    o.transientTime = v.at("transient_time").number(o.transientTime);
    o.outputStep = v.at("output_step").number(o.outputStep);
    o.maximumStep = v.at("maximum_step").number(o.maximumStep);
    o.sectionValue = v.at("section_value")
                         .number(o.sectionAxis < o.lower.size() && o.sectionAxis < o.upper.size()
                                     ? (o.lower[o.sectionAxis] + o.upper[o.sectionAxis]) / 2
                                     : 0);
    o.sectionNormal = numbers(v.at("section_normal"));
    o.sectionReference = numbers(v.at("section_reference"));
    const double direction=v.at("section_direction").number(o.sectionDirection);
    if(direction < -1 || direction > 1 || direction != std::floor(direction))
        throw std::invalid_argument("Invalid Poincare section direction");
    o.sectionDirection=int(direction);
    o.sectionMinimumTime=v.at("section_minimum_time").number(o.sectionMinimumTime);
    o.transversalityTolerance=v.at("transversality_tolerance")
        .number(o.transversalityTolerance);
    if (!v.at("vector_axes").isNull()) {
        o.vectorAxes.clear();
        for (double q : numbers(v.at("vector_axes"))) {
            if (q < 0 || q >= dimension || q != std::floor(q))
                throw std::invalid_argument("Invalid vector axis");
            o.vectorAxes.push_back(std::size_t(q));
        }
    }
    return o;
}
namespace detail {
json::Value phaseSpaceResultJson(const dynsyslab::core::PhaseSpaceResult &r) {
    json::Value::Array points, manifolds, surfaces, orbits, returns, branches, diagnostics, axes, checks;
    for (const auto &c : r.knowledgeChecks)
        checks.emplace_back(json::Value::Object{{"fact_id", c.factId},
                                                {"kind", c.kind},
                                                {"status", c.status},
                                                {"reason", c.reason},
                                                {"state", array(c.state)},
                                                {"object_index", double(c.objectIndex)},
                                                {"residual", c.residual}});
    for (const auto &p : r.fixedPoints)
        points.emplace_back(json::Value::Object{{"state", array(p.state)},
                                                {"residual", p.residual},
                                                {"stability", p.stability},
                                                {"eigenvalues", complexArray(p.eigenvalues)},
                                                {"eigenvectors", complexArray(p.eigenvectors)},
                                                {"stable_dimension", double(p.stableDimension)},
                                                {"unstable_dimension", double(p.unstableDimension)},
                                                {"center_dimension", double(p.centerDimension)}});
    for (const auto &m : r.manifolds)
        manifolds.push_back(curveJson(m));
    for(const auto &m:r.manifoldSurfaces){
        json::Value::Array indices,offsets;
        for(auto q:m.triangles)indices.emplace_back(double(q));
        for(auto q:m.frontOffsets)offsets.emplace_back(double(q));
        surfaces.emplace_back(json::Value::Object{{"kind",m.kind},{"status",m.status},{"reason",m.reason},
            {"base_kind",m.baseKind},{"base_object",double(m.baseObject)},{"manifold_dimension",double(m.manifoldDimension)},
            {"vertices",array(m.vertices)},{"triangles",indices},{"front_offsets",offsets},{"front_distances",array(m.frontDistances)},
            {"shooting_residual",m.shootingResidual},{"refinement_error",m.refinementError},{"tangent_defect",m.tangentDefect},
            {"seed_radius",m.seedRadius},{"maximum_turn",m.maximumTurn},
            {"method","adaptive boundary-value fronts; full-state mesh; seed/mesh refinement and face-centroid field tangency"}});
    }
    for (const auto &p : r.periodicOrbits)
        orbits.emplace_back(
            json::Value::Object{{"curve", curveJson(p.curve)},
                                {"period", p.period},
                                {"section_state", array(p.sectionState)},
                                {"stable_dimension", double(p.stableDimension)},
                                {"unstable_dimension", double(p.unstableDimension)},
                                {"center_dimension", double(p.centerDimension)},
                                {"floquet_vectors", complexArray(p.floquetVectors)},
                                {"closure_residual", p.closureResidual},
                                {"phase_residual", p.phaseResidual},
                                {"neutral_multiplier_error", r.options.dynamics == "periodic-forcing" ? json::Value() : json::Value(p.neutralMultiplierError)},
                                {"tangent_identity_error", r.options.dynamics == "periodic-forcing" ? json::Value() : json::Value(p.tangentIdentityError)},
                                {"variational_refinement_error", p.variationalRefinementError},
                                {"floquet_multipliers", complexArray(p.floquetMultipliers)},
                                {"section_multipliers", complexArray(p.sectionMultipliers)},
                                {"return_time", p.returnTime},
                                {"return_residual", p.returnResidual},
                                {"transversality", r.options.dynamics == "periodic-forcing" ? json::Value() : json::Value(p.transversality)},
                                {"transversality_angle", r.options.dynamics == "periodic-forcing" ? json::Value() : json::Value(p.transversalityAngle)},
                                {"section_tangent_basis", array(p.sectionTangentBasis)},
                                {"return_derivative", array(p.returnDerivative)},
                                {"reduced_return_derivative", array(p.reducedReturnDerivative)},
                                {"stability", p.stability},
                                {"isolated", p.isolated}});
    for (const auto &map : r.poincareMaps)
        returns.push_back(poincareJson(map));
    for (const auto &branch : r.bifurcationBranches) {
        json::Value::Array branchPoints, events, branchDiagnostics;
        for (const auto &p : branch.points)
            branchPoints.emplace_back(json::Value::Object{
                {"parameter",p.parameter},{"state",array(p.state)},{"period",p.period},
                {"residual",p.residual},{"arclength_step",p.arclengthStep},
                {"correction_norm",p.correctionNorm},{"spectrum",complexArray(p.spectrum)},
                {"stability",p.stability}});
        for (const auto &event : branch.events)
            events.emplace_back(json::Value::Object{
                {"kind",event.kind},{"status",event.status},{"reason",event.reason},
                {"left_index",double(event.leftIndex)},{"right_index",double(event.rightIndex)},
                {"parameter",event.parameter},{"state",array(event.state)},
                {"indicator",event.indicator},{"refinement_error",event.refinementError}});
        for (const auto &message : branch.diagnostics)
            branchDiagnostics.emplace_back(message);
        branches.emplace_back(json::Value::Object{
            {"kind",branch.kind},{"parameter",branch.parameter},{"status",branch.status},
            {"reason",branch.reason},{"points",branchPoints},{"events",events},
            {"diagnostics",branchDiagnostics}});
    }
    for (const auto &d : r.diagnostics)
        diagnostics.emplace_back(d);
    for (auto axis : r.vectorField.axes)
        axes.emplace_back(double(axis));
    return json::Value::Object{
        {"schema", r.schema},
        {"method_version", r.methodVersion},
        {"equations_fingerprint", r.equationsFingerprint},
        {"coordinate_system", r.coordinateSystem},
        {"state_dimension", double(r.options.lower.size())},
        {"options", phaseSpaceOptionsJson(r.options)},
        {"fixed_points", points},
        {"manifolds", manifolds},
        {"manifold_surfaces", surfaces},
        {"periodic_orbits", orbits},
        {"poincare_maps", returns},
        {"bifurcation_branches", branches},
        {"knowledge_checks", checks},
        {"vector_field", json::Value::Object{{"points", array(r.vectorField.points)},
                                             {"vectors", array(r.vectorField.vectors)},
                                             {"axes", axes},
                                             {"slice_state", array(r.vectorField.sliceState)}}},
        {"diagnostics", diagnostics},
        {"equilibrium_searches", double(r.equilibriumSearches)},
        {"equilibrium_failures", double(r.equilibriumFailures)},
        {"periodic_searches", double(r.periodicSearches)},
        {"periodic_failures", double(r.periodicFailures)}};
}
} // namespace detail
} // namespace dynsyslab::workflow
