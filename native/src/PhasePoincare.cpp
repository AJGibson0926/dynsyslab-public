#include "PhaseSpaceInternal.hpp"

#include <algorithm>
#include <cmath>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core::phase_detail {
namespace {
double dot(const PhaseState &a, const PhaseState &b) {
    double value = 0;
    for (std::size_t j = 0; j < a.size(); ++j)
        value += a[j] * b[j];
    return value;
}

PhaseState sample(const PhaseCurve &curve, std::size_t row, std::size_t dimension) {
    return {curve.values.begin() + row * dimension,
            curve.values.begin() + (row + 1) * dimension};
}

bool crosses(double left, double right, int direction, double tolerance) {
    if (direction > 0)
        return left < -tolerance && right >= -tolerance;
    if (direction < 0)
        return left > tolerance && right <= tolerance;
    return (left < -tolerance && right >= -tolerance) ||
           (left > tolerance && right <= tolerance);
}

double matrixDifference(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size())
        return std::numeric_limits<double>::infinity();
    double result = 0;
    for (std::size_t j = 0; j < a.size(); ++j)
        result = std::max(result, std::abs(a[j] - b[j]) / (1 + std::abs(b[j])));
    return result;
}

struct Bracket {
    double left{}, right{}, leftValue{}, rightValue{};
};

Bracket returnBracket(const ODESystem &system, const PhaseState &initial,
                      const PhaseSpaceOptions &options, const SectionGeometry &section,
                      const ExecutionContext *execution) {
    auto scanOptions = options;
    scanOptions.outputStep = std::min(options.outputStep, options.maximumStep);
    const auto scan = flow(system, initial, options.horizon, scanOptions, false, true, execution);
    const auto n = initial.size();
    const double eventTolerance = std::max(options.absoluteTolerance,
                                           options.residualTolerance * .01);
    const double minimumTime = std::max(options.sectionMinimumTime,
                                        32 * std::numeric_limits<double>::epsilon() *
                                            (1 + options.horizon));
    double previous = sectionValue(initial, section);
    for (std::size_t row = 1; row < scan.curve.time.size(); ++row) {
        const auto currentState = sample(scan.curve, row, n);
        const double current = sectionValue(currentState, section);
        if (scan.curve.time[row] > minimumTime &&
            crosses(previous, current, options.sectionDirection, eventTolerance))
            return {scan.curve.time[row - 1], scan.curve.time[row], previous, current};
        previous = current;
    }
    throw std::runtime_error(scan.complete
                                 ? "No oriented first positive section return within horizon"
                                 : "Section-return scan truncated: " + scan.reason);
}

Flow refineReturn(const ODESystem &system, const PhaseState &initial,
                  const PhaseSpaceOptions &options, const SectionGeometry &section,
                  Bracket bracket, const ExecutionContext *execution) {
    Flow best;
    const double timeTolerance = 16 * std::numeric_limits<double>::epsilon() *
                                 (1 + std::abs(bracket.right));
    for (int iteration = 0; iteration < 64; ++iteration) {
        if (execution)
            execution->checkpoint();
        const double candidateTime = (bracket.left + bracket.right) / 2;
        auto candidate = flow(system, initial, candidateTime, options, true, false, execution);
        if (!candidate.complete)
            throw std::runtime_error("Section-return refinement failed: " + candidate.reason);
        const double value = sectionValue(candidate.state, section);
        best = std::move(candidate);
        if (std::abs(value) <= options.absoluteTolerance ||
            bracket.right - bracket.left <= timeTolerance)
            break;
        if ((bracket.leftValue <= 0 && value >= 0) ||
            (bracket.leftValue >= 0 && value <= 0)) {
            bracket.right = candidateTime;
            bracket.rightValue = value;
        } else {
            bracket.left = candidateTime;
            bracket.leftValue = value;
        }
    }
    return best;
}
} // namespace

SectionGeometry sectionGeometry(const PhaseSpaceOptions &o, std::size_t n) {
    SectionGeometry section;
    section.reference.assign(n, 0);
    PhaseState normal(n, 0);
    if (!o.sectionNormal.empty()) {
        if (o.sectionNormal.size() != n || o.sectionReference.size() != n)
            throw std::invalid_argument(
                "Affine Poincare section requires normal and reference matching the state dimension");
        normal = o.sectionNormal;
        section.reference = o.sectionReference;
    } else {
        if (o.sectionAxis >= n)
            throw std::invalid_argument("Poincare section axis is outside the state dimension");
        normal[o.sectionAxis] = 1;
        section.reference[o.sectionAxis] = o.sectionValue;
    }
    double normalNorm = std::sqrt(dot(normal, normal));
    if (!std::isfinite(normalNorm) || normalNorm <= 0)
        throw std::invalid_argument("Poincare section normal must be finite and nonzero");
    section.physicalNormal = normal;
    for (auto &value : section.physicalNormal)
        value /= normalNorm;
    section.scaledNormal.resize(n);
    for (std::size_t j = 0; j < n; ++j)
        section.scaledNormal[j] = section.physicalNormal[j] * (o.upper[j] - o.lower[j]);

    const std::size_t tangentDimension = n > 0 ? n - 1 : 0;
    section.tangentBasis.assign(n * tangentDimension, 0);
    std::size_t column = 0;
    for (std::size_t candidate = 0; candidate < n && column < tangentDimension; ++candidate) {
        PhaseState q(n, 0);
        q[candidate] = 1;
        const double normalComponent = dot(q, section.physicalNormal);
        for (std::size_t row = 0; row < n; ++row)
            q[row] -= normalComponent * section.physicalNormal[row];
        for (std::size_t previous = 0; previous < column; ++previous) {
            double projection = 0;
            for (std::size_t row = 0; row < n; ++row)
                projection += q[row] * section.tangentBasis[row * tangentDimension + previous];
            for (std::size_t row = 0; row < n; ++row)
                q[row] -= projection * section.tangentBasis[row * tangentDimension + previous];
        }
        const double qNorm = std::sqrt(dot(q, q));
        if (qNorm <= 64 * std::numeric_limits<double>::epsilon())
            continue;
        for (std::size_t row = 0; row < n; ++row)
            section.tangentBasis[row * tangentDimension + column] = q[row] / qNorm;
        ++column;
    }
    if (column != tangentDimension)
        throw std::runtime_error("Unable to construct Poincare-section tangent coordinates");
    return section;
}

double sectionValue(const PhaseState &state, const SectionGeometry &section) {
    double value = 0;
    for (std::size_t j = 0; j < state.size(); ++j)
        value += section.physicalNormal[j] * (state[j] - section.reference[j]);
    return value;
}

PhasePoincareMap poincareFromTransition(const ODESystem &system, const PhaseState &initial,
                                        const Flow &transition, const PhaseSpaceOptions &o,
                                        const SectionGeometry &section) {
    const auto n = initial.size(), k = n - 1;
    if (!transition.complete || transition.monodromy.size() != n * n)
        throw std::invalid_argument("Poincare derivative requires a complete variational transition");
    PhasePoincareMap result;
    result.initialState = initial;
    result.sectionState = initial;
    result.returnState = transition.state;
    result.normal = section.physicalNormal;
    result.tangentBasis = section.tangentBasis;
    result.returnTime = transition.curve.time.back() - transition.curve.time.front();
    result.fixedPointResidual = norm([&] {
        auto difference = transition.state;
        for (std::size_t j = 0; j < n; ++j)
            difference[j] -= initial[j];
        return difference;
    }(), o);
    const auto endField = field(system, transition.state);
    result.transversality = dot(section.physicalNormal, endField);
    double normalScale = 0, fieldScale = 0;
    for (std::size_t j = 0; j < n; ++j) {
        normalScale += section.scaledNormal[j] * section.scaledNormal[j];
        const double scaledField = endField[j] / (o.upper[j] - o.lower[j]);
        fieldScale += scaledField * scaledField;
    }
    const double denominatorScale = std::sqrt(normalScale * fieldScale);
    result.transversalityAngle = denominatorScale > 0
                                    ? std::abs(result.transversality) / denominatorScale
                                    : 0;
    result.direction = result.transversality > 0 ? 1 : (result.transversality < 0 ? -1 : 0);
    if (!std::isfinite(result.transversalityAngle) ||
        result.transversalityAngle <= o.transversalityTolerance) {
        result.status = "grazing";
        result.reason = "Section-return derivative is singular or ill-conditioned at a tangency";
        return result;
    }

    result.returnTimeDerivative.resize(n);
    result.returnDerivative.assign(n * n, 0);
    for (std::size_t column = 0; column < n; ++column) {
        double normalMonodromy = 0;
        for (std::size_t row = 0; row < n; ++row)
            normalMonodromy +=
                section.physicalNormal[row] * transition.monodromy[row * n + column];
        result.returnTimeDerivative[column] = -normalMonodromy / result.transversality;
        for (std::size_t row = 0; row < n; ++row)
            result.returnDerivative[row * n + column] =
                transition.monodromy[row * n + column] +
                endField[row] * result.returnTimeDerivative[column];
    }
    result.reducedDerivative.assign(k * k, 0);
    for (std::size_t row = 0; row < k; ++row)
        for (std::size_t column = 0; column < k; ++column)
            for (std::size_t a = 0; a < n; ++a)
                for (std::size_t b = 0; b < n; ++b)
                    result.reducedDerivative[row * k + column] +=
                        section.tangentBasis[a * k + row] *
                        result.returnDerivative[a * n + b] *
                        section.tangentBasis[b * k + column];
    if (k)
        result.multipliers =
            dynsyslab::numerics::eigensystem(result.reducedDerivative, k, false).eigenvalues;
    result.status = "complete";
    return result;
}

void findPoincareReturns(const ODESystem &system, PhaseSpaceResult &result,
                         const ExecutionContext *execution) {
    for (const auto &seed : result.options.sectionSeeds) {
        try {
            result.poincareMaps.push_back(
                firstPoincareReturn(system, seed, result.options, execution));
        } catch (const ExecutionCancelled &) {
            throw;
        } catch (const std::exception &error) {
            PhasePoincareMap failed;
            failed.initialState = seed;
            failed.sectionState = seed;
            failed.status = "failed";
            failed.reason = error.what();
            result.poincareMaps.push_back(std::move(failed));
            result.diagnostics.push_back(std::string("Poincare seed: ") + error.what());
        }
    }
}
} // namespace dynsyslab::core::phase_detail

namespace dynsyslab::core {
PhasePoincareMap firstPoincareReturn(const ODESystem &system, const PhaseState &seed,
                                     const PhaseSpaceOptions &options,
                                     const ExecutionContext *execution) {
    if (options.dynamics != "autonomous-zero-input")
        throw std::invalid_argument(
            "State-defined Poincare returns require autonomous dynamics; use the forcing section for periodic forcing");
    if (seed.size() != system.descriptor().stateDimension || !phase_detail::inside(seed, options))
        throw std::invalid_argument("Poincare seed must lie inside the declared phase-space box");
    const auto section = phase_detail::sectionGeometry(options, seed.size());
    const double sectionResidual = std::abs(phase_detail::sectionValue(seed, section));
    if (sectionResidual > std::max(options.residualTolerance,
                                   10 * options.absoluteTolerance))
        throw std::invalid_argument("Poincare seed does not lie on the configured section");
    const auto bracket =
        phase_detail::returnBracket(system, seed, options, section, execution);
    auto coarse =
        phase_detail::refineReturn(system, seed, options, section, bracket, execution);
    auto refinedOptions = options;
    refinedOptions.maximumStep *= .5;
    refinedOptions.relativeTolerance *= .125;
    refinedOptions.absoluteTolerance *= .125;
    const auto refinedBracket =
        phase_detail::returnBracket(system, seed, refinedOptions, section, execution);
    auto refined = phase_detail::refineReturn(system, seed, refinedOptions, section,
                                               refinedBracket, execution);
    auto result = phase_detail::poincareFromTransition(system, seed, refined, options, section);
    result.endpointRefinementError = phase_detail::norm([&] {
        auto difference = refined.state;
        for (std::size_t j = 0; j < difference.size(); ++j)
            difference[j] -= coarse.state[j];
        return difference;
    }(), options);
    const auto coarseMap =
        phase_detail::poincareFromTransition(system, seed, coarse, options, section);
    result.derivativeRefinementError = phase_detail::matrixDifference(
        coarseMap.reducedDerivative, result.reducedDerivative);
    if (result.status == "complete" &&
        (result.endpointRefinementError > 20 * options.residualTolerance ||
         result.derivativeRefinementError > 1e-5)) {
        result.status = "unresolved";
        result.reason = "Section return failed tighter-integration endpoint or derivative validation";
    }
    return result;
}
} // namespace dynsyslab::core
