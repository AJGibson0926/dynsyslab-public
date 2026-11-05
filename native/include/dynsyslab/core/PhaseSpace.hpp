#pragma once
#include <complex>
#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>
#include <functional>
#include <memory>
namespace dynsyslab::core {
using PhaseState = std::vector<double>;
struct PeriodicOrbitSeed {
    PhaseState state;
    double period{};
};
struct PhaseSurfaceOptions {
    bool enabled{true};
    std::size_t initialPoints{16}, maximumPoints{4096}, maximumFronts{64};
    double seedRadius{.002}, maximumDistance{.5}, step{.02}, meshSize{.04};
    double chordTolerance{.001}, tangentTolerance{.08};
};
struct PhaseBifurcationOptions {
    bool enabled{};
    std::string branch{"equilibrium"}; // equilibrium or periodic-orbit
    std::string parameter;
    double minimum{}, maximum{}, initial{}, step{.02}, minimumStep{1e-5}, maximumStep{.1};
    std::size_t maximumPoints{128}, maximumCorrectorIterations{16};
    PhaseState seedState;
    double seedPeriod{};
};
struct PhaseSpaceOptions {
    PhaseSurfaceOptions surface;
    PhaseBifurcationOptions bifurcation;
    PhaseState lower, upper;
    std::string dynamics{"autonomous-zero-input"};
    double forcingPeriod{}, forcingSectionTime{}; // native model time; zero period infers the clock
    std::size_t forcingPeriods{1}, manifoldIterations{16}, maximumManifoldPoints{1024};
    double manifoldTolerance{2e-4}; // domain-scaled chord and interpolation tolerance
    std::vector<PhaseState> equilibriumSeeds;
    std::vector<PeriodicOrbitSeed> periodicSeeds;
    std::vector<PhaseState> sectionSeeds;
    bool fixedPoints{true}, manifolds{true}, periodicOrbits{true}, poincareMap{true},
        vectorField{true};
    std::size_t gridPerAxis{7}, maximumSeeds{256}, vectorGridPerAxis{11};
    std::size_t maximumNewtonIterations{30}, maximumIntegrationSteps{200000};
    double residualTolerance{1e-8}, relativeTolerance{1e-9}, absoluteTolerance{1e-11};
    double eigenvalueTolerance{1e-7}, seedOffset{1e-4}; // offset in domain-scaled coordinates
    double horizon{40}, transientTime{30}, outputStep{0.05}, maximumStep{0.05};
    std::size_t sectionAxis{0};
    double sectionValue{0};
    PhaseState sectionNormal, sectionReference; // optional affine section n.(x-reference)=0
    int sectionDirection{1}; // -1, 0 (either), or +1
    double sectionMinimumTime{}, transversalityTolerance{1e-7};
    std::vector<std::size_t> vectorAxes{0, 1};
    PhaseState sliceState; // explicit values for unplotted dimensions
};
struct PhaseFixedPoint {
    PhaseState state;
    std::vector<std::complex<double>> eigenvalues, eigenvectors; // columns, row-major matrix
    double residual{};
    std::string stability;
    std::size_t stableDimension{}, unstableDimension{}, centerDimension{};
};
struct PhaseCurve {
    std::string kind, status, reason;
    std::string baseKind{"equilibrium"};
    std::size_t baseObject{};
    std::vector<double> time, values; // model coordinates, sample-major
    double seedOffset{}, refinementError{}, tangentDefect{};
};
struct PhasePeriodicOrbit {
    PhaseCurve curve;
    double period{}, closureResidual{}, phaseResidual{}, neutralMultiplierError{};
    std::vector<std::complex<double>> floquetMultipliers;
    std::string stability;
    std::size_t stableDimension{}, unstableDimension{}, centerDimension{};
    PhaseState sectionState;
    std::vector<std::complex<double>> floquetVectors;
    std::vector<std::complex<double>> sectionMultipliers;
    std::vector<double> returnDerivative, reducedReturnDerivative, sectionTangentBasis;
    double returnTime{}, returnResidual{}, transversality{}, transversalityAngle{};
    double tangentIdentityError{}, variationalRefinementError{};
    bool isolated{};
};
struct PhasePoincareMap {
    std::string kind{"first-positive-return"}, status, reason;
    PhaseState initialState, sectionState, returnState, normal;
    std::vector<double> tangentBasis, returnDerivative, reducedDerivative, returnTimeDerivative;
    std::vector<std::complex<double>> multipliers;
    double returnTime{}, fixedPointResidual{}, transversality{}, transversalityAngle{};
    double endpointRefinementError{}, derivativeRefinementError{};
    int direction{};
};
struct PhaseBifurcationPoint {
    double parameter{}, period{}, residual{}, arclengthStep{}, correctionNorm{};
    PhaseState state;
    std::vector<std::complex<double>> spectrum;
    std::string stability;
};
struct PhaseBifurcationEvent {
    std::string kind, status{"candidate"}, reason;
    std::size_t leftIndex{}, rightIndex{};
    double parameter{}, indicator{}, refinementError{};
    PhaseState state;
};
struct PhaseBifurcationBranch {
    std::string kind, parameter, status, reason;
    std::vector<PhaseBifurcationPoint> points;
    std::vector<PhaseBifurcationEvent> events;
    std::vector<std::string> diagnostics;
};
struct PhaseSurface {
    std::string kind, status, reason, baseKind{"equilibrium"};
    std::size_t baseObject{}, manifoldDimension{2};
    PhaseState vertices; // full native coordinates, vertex-major
    std::vector<std::size_t> triangles, frontOffsets;
    std::vector<double> frontDistances; // accumulated chord-distance estimates
    double shootingResidual{}, refinementError{}, tangentDefect{};
    double seedRadius{}, maximumTurn{};
};
struct PhaseVectorField {
    std::vector<double> points, vectors;
    std::vector<std::size_t> axes;
    PhaseState sliceState;
};
struct PhaseKnowledgeCheck {
    std::string factId, kind, status, reason;
    PhaseState state;
    std::size_t objectIndex{};
    double residual{};
};
struct PhaseSpaceResult {
    std::string schema{"dynsyslab.phase-space/v1"}, methodVersion{"4"};
    std::string equationsFingerprint, coordinateSystem;
    PhaseSpaceOptions options;
    std::vector<PhaseFixedPoint> fixedPoints;
    std::vector<PhaseCurve> manifolds;
    std::vector<PhaseSurface> manifoldSurfaces;
    std::vector<PhasePeriodicOrbit> periodicOrbits;
    std::vector<PhasePoincareMap> poincareMaps;
    std::vector<PhaseBifurcationBranch> bifurcationBranches;
    PhaseVectorField vectorField;
    std::vector<std::string> diagnostics;
    std::vector<PhaseKnowledgeCheck> knowledgeChecks;
    std::size_t equilibriumSearches{}, equilibriumFailures{}, periodicSearches{},
        periodicFailures{};
};
using PhaseSystemFactory = std::function<std::unique_ptr<ODESystem>(double)>;
// Resolve a common native-time period from the configured sinusoidal clocks.
// Explicit periods must be integer multiples of each active clock.
double phaseSpaceForcingPeriod(const ODESystem &system, double requested = 0);
// Explicit call only. Autonomous zero input or fixed-phase periodic forcing on a box.
// Returned numerical candidates do not certify global completeness or basin boundaries.
void validatePhaseSpaceOptions(const ODESystem &system, const PhaseSpaceOptions &options);
PhaseSpaceResult analyzePhaseSpace(const ODESystem &system, const PhaseSpaceOptions &options,
                                   const ExecutionContext *execution = nullptr);
// First oriented positive return to the configured affine section. The seed must lie on it.
PhasePoincareMap firstPoincareReturn(const ODESystem &system, const PhaseState &sectionSeed,
                                     const PhaseSpaceOptions &options,
                                     const ExecutionContext *execution = nullptr);
// Pseudo-arclength continuation over a caller-supplied one-parameter system family.
PhaseBifurcationBranch continuePhaseBifurcation(const PhaseSystemFactory &factory,
                                                const PhaseSpaceOptions &options,
                                                const PhaseBifurcationOptions &bifurcation,
                                                const ExecutionContext *execution = nullptr);
} // namespace dynsyslab::core
