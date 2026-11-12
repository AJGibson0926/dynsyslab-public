#include "PhaseSpaceInternal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace dynsyslab::core::phase_detail {
namespace {
struct Indicator {
    bool available{};
    double value{};
};

bool bracketed(double left, double right, double tolerance) {
    if (!std::isfinite(left) || !std::isfinite(right))
        return false;
    if (std::abs(left) <= tolerance && std::abs(right) <= tolerance)
        return false;
    return left * right <= 0;
}

Indicator closestReal(const std::vector<std::complex<double>> &spectrum, double target,
                      double tolerance) {
    Indicator result;
    double distance = std::numeric_limits<double>::infinity();
    for (const auto value : spectrum) {
        if (std::abs(value.imag()) > 100 * tolerance)
            continue;
        const double candidate = value.real() - target;
        if (std::abs(candidate) < distance) {
            result = {true, candidate};
            distance = std::abs(candidate);
        }
    }
    return result;
}

Indicator closestComplexRealPart(const std::vector<std::complex<double>> &spectrum,
                                 double tolerance) {
    Indicator result;
    double distance = std::numeric_limits<double>::infinity();
    for (const auto value : spectrum) {
        if (std::abs(value.imag()) <= 100 * tolerance)
            continue;
        if (std::abs(value.real()) < distance) {
            result = {true, value.real()};
            distance = std::abs(value.real());
        }
    }
    return result;
}

Indicator closestComplexRadius(const std::vector<std::complex<double>> &spectrum,
                               double tolerance) {
    Indicator result;
    double distance = std::numeric_limits<double>::infinity();
    for (const auto value : spectrum) {
        if (std::abs(value.imag()) <= 100 * tolerance)
            continue;
        const double candidate = std::abs(value) - 1;
        if (std::abs(candidate) < distance) {
            result = {true, candidate};
            distance = std::abs(candidate);
        }
    }
    return result;
}

double scalarField(const PhaseSystemFactory &factory, double parameter, double state) {
    auto system = factory(parameter);
    PhaseState x{state};
    return field(*system, x).front();
}

std::pair<std::string, std::string>
classifyScalar(const PhaseSystemFactory &factory, const PhaseSpaceOptions &options,
               const PhaseBifurcationOptions &bifurcation, double state, double parameter) {
    const double width = options.upper.front() - options.lower.front();
    const double range = bifurcation.maximum - bifurcation.minimum;
    const double hx = std::max(1e-4 * width,
                               std::cbrt(std::numeric_limits<double>::epsilon()) *
                                   std::max(1., std::abs(state)));
    const double hp = std::max(1e-4 * range,
                               std::cbrt(std::numeric_limits<double>::epsilon()) *
                                   std::max(1., std::abs(parameter)));
    if (state - 2 * hx < options.lower.front() || state + 2 * hx > options.upper.front() ||
        parameter - hp < bifurcation.minimum || parameter + hp > bifurcation.maximum)
        return {"scalar_zero_eigenvalue_candidate",
                "The local derivative stencil reached the declared continuation boundary"};
    const double f0 = scalarField(factory, parameter, state);
    const double fxp = scalarField(factory, parameter, state + hx);
    const double fxm = scalarField(factory, parameter, state - hx);
    const double fxpp = scalarField(factory, parameter, state + 2 * hx);
    const double fxmm = scalarField(factory, parameter, state - 2 * hx);
    const double fpp = scalarField(factory, parameter + hp, state);
    const double fpm = scalarField(factory, parameter - hp, state);
    const double fMu = (fpp - fpm) / (2 * hp);
    const double fXx = (fxp - 2 * f0 + fxm) / (hx * hx);
    const double fXxx = (fxpp - 2 * fxp + 2 * fxm - fxmm) / (2 * hx * hx * hx);
    const double jacobianPlus = factory(parameter + hp)->jacobian(0, &state).front();
    const double jacobianMinus = factory(parameter - hp)->jacobian(0, &state).front();
    const double fXMu = (jacobianPlus - jacobianMinus) / (2 * hp);
    const double residualScale = std::max({1., std::abs(f0), std::abs(fxp), std::abs(fxm)});
    const double structuralTolerance =
        std::max(1e-6 * residualScale, 100 * options.residualTolerance);
    const bool preserved = std::max(std::abs(fpp), std::abs(fpm)) <=
                           std::max(1e-7 * residualScale, 1000 * options.residualTolerance);
    const double symmetryStep = std::min(.01 * width,
        .45 * std::min(state - options.lower.front(), options.upper.front() - state));
    bool reflectionOdd = false;
    if (symmetryStep > hx) {
        const double positive = scalarField(factory, parameter, state + symmetryStep);
        const double negative = scalarField(factory, parameter, state - symmetryStep);
        reflectionOdd = std::abs(positive + negative - 2 * f0) <=
                        std::max(1e-7 * (1 + std::abs(positive) + std::abs(negative)),
                                 1000 * options.residualTolerance);
    }
    std::string kind = "scalar_zero_eigenvalue_candidate";
    if (std::abs(fMu) > structuralTolerance && std::abs(fXx) > structuralTolerance)
        kind = "fold";
    else if (preserved && reflectionOdd && std::abs(fXMu) > structuralTolerance &&
             std::abs(fXxx) > structuralTolerance)
        kind = "pitchfork";
    else if (preserved && std::abs(fXMu) > structuralTolerance &&
             std::abs(fXx) > structuralTolerance)
        kind = "transcritical";
    std::ostringstream evidence;
    evidence << "finite-difference jet: f_mu=" << fMu << ", f_xx=" << fXx
             << ", f_xmu=" << fXMu << ", f_xxx=" << fXxx
             << "; preserved_branch=" << (preserved ? "true" : "false")
             << "; reflection_odd=" << (reflectionOdd ? "true" : "false");
    return {kind, evidence.str()};
}

void addEvent(const PhaseSystemFactory &factory, const PhaseSpaceOptions &options,
              const PhaseBifurcationOptions &bifurcation, PhaseBifurcationBranch &branch,
              std::size_t rightIndex, const std::string &kind, double leftIndicator,
              double rightIndicator, bool scalarClassification) {
    const auto &left = branch.points[rightIndex - 1];
    const auto &right = branch.points[rightIndex];
    const double denominator = std::abs(leftIndicator) + std::abs(rightIndicator);
    const double fraction = denominator > 0 ? std::abs(leftIndicator) / denominator : .5;
    PhaseBifurcationEvent event;
    event.kind = kind;
    event.leftIndex = rightIndex - 1;
    event.rightIndex = rightIndex;
    event.parameter = left.parameter + fraction * (right.parameter - left.parameter);
    event.state.resize(left.state.size());
    for (std::size_t j = 0; j < event.state.size(); ++j)
        event.state[j] = left.state[j] + fraction * (right.state[j] - left.state[j]);
    event.indicator = leftIndicator + fraction * (rightIndicator - leftIndicator);
    event.refinementError = std::abs(right.parameter - left.parameter);
    event.reason = "Spectral sign change bracketed by adjacent continuation points";
    if (scalarClassification) {
        const auto classified =
            classifyScalar(factory, options, bifurcation, event.state.front(), event.parameter);
        event.kind = classified.first;
        event.reason = classified.second;
        event.status = event.kind == "scalar_zero_eigenvalue_candidate"
                           ? "candidate"
                           : "structurally-classified-candidate";
    }
    branch.events.push_back(std::move(event));
}
} // namespace

void detectBifurcationEvents(const PhaseSystemFactory &factory,
                             const PhaseSpaceOptions &options,
                             const PhaseBifurcationOptions &bifurcation,
                             PhaseBifurcationBranch &branch) {
    const double tolerance = options.eigenvalueTolerance;
    for (std::size_t right = 1; right < branch.points.size(); ++right) {
        const auto &leftPoint = branch.points[right - 1];
        const auto &rightPoint = branch.points[right];
        if (bifurcation.branch == "equilibrium") {
            const auto leftZero = closestReal(leftPoint.spectrum, 0, tolerance);
            const auto rightZero = closestReal(rightPoint.spectrum, 0, tolerance);
            if (leftZero.available && rightZero.available &&
                bracketed(leftZero.value, rightZero.value, tolerance))
                addEvent(factory, options, bifurcation, branch, right,
                         "equilibrium_zero_eigenvalue_candidate", leftZero.value,
                         rightZero.value, options.lower.size() == 1);
            const auto leftHopf = closestComplexRealPart(leftPoint.spectrum, tolerance);
            const auto rightHopf = closestComplexRealPart(rightPoint.spectrum, tolerance);
            if (leftHopf.available && rightHopf.available &&
                bracketed(leftHopf.value, rightHopf.value, tolerance))
                addEvent(factory, options, bifurcation, branch, right, "hopf_candidate",
                         leftHopf.value, rightHopf.value, false);
            continue;
        }
        for (const auto target : {1., -1.}) {
            const auto left = closestReal(leftPoint.spectrum, target, tolerance);
            const auto rightValue = closestReal(rightPoint.spectrum, target, tolerance);
            if (left.available && rightValue.available &&
                bracketed(left.value, rightValue.value, tolerance))
                addEvent(factory, options, bifurcation, branch, right,
                         target > 0 ? "periodic_plus_one_multiplier_candidate"
                                    : "period_doubling_candidate",
                         left.value, rightValue.value, false);
        }
        const auto leftCircle = closestComplexRadius(leftPoint.spectrum, tolerance);
        const auto rightCircle = closestComplexRadius(rightPoint.spectrum, tolerance);
        if (leftCircle.available && rightCircle.available &&
            bracketed(leftCircle.value, rightCircle.value, tolerance))
            addEvent(factory, options, bifurcation, branch, right,
                     "neimark_sacker_candidate", leftCircle.value, rightCircle.value, false);
    }
}
} // namespace dynsyslab::core::phase_detail
