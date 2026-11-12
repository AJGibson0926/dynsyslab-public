#include "PhaseSpaceInternal.hpp"

#include <algorithm>
#include <cmath>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
struct Evaluation {
    PhaseState residual;
    std::vector<double> jacobian; // residual rows by physical variables, including parameter
    std::vector<std::complex<double>> spectrum;
    std::string stability;
    double period{};
};

double maximumAbs(const PhaseState &values) {
    double result = 0;
    for (double value : values)
        result = std::max(result, std::abs(value));
    return result;
}

PhaseState subtract(PhaseState a, const PhaseState &b) {
    for (std::size_t j = 0; j < a.size(); ++j)
        a[j] -= b[j];
    return a;
}

std::string classifySpectrum(const std::vector<std::complex<double>> &spectrum,
                             bool periodic, double tolerance) {
    std::size_t stable = 0, unstable = 0, center = 0;
    for (const auto value : spectrum) {
        const double indicator = periodic ? std::abs(value) - 1 : value.real();
        if (indicator < -10 * tolerance)
            ++stable;
        else if (indicator > 10 * tolerance)
            ++unstable;
        else
            ++center;
    }
    if (periodic) {
        if (center)
            return "nonhyperbolic_periodic_orbit";
        if (stable == spectrum.size())
            return "attracting_periodic_orbit";
        if (unstable == spectrum.size())
            return "repelling_periodic_orbit";
        return "saddle_periodic_orbit";
    }
    if (center)
        return "nonhyperbolic";
    if (stable == spectrum.size())
        return "sink";
    if (unstable == spectrum.size())
        return "source";
    return "saddle";
}

class BranchProblem {
  public:
    BranchProblem(const PhaseSystemFactory &factory, const PhaseSpaceOptions &options,
                  const PhaseBifurcationOptions &bifurcation)
        : factory_(factory), options_(options), bifurcation_(bifurcation),
          n_(options.lower.size()), periodic_(bifurcation.branch == "periodic-orbit"),
          autonomous_(periodic_ && options.dynamics == "autonomous-zero-input"),
          section_(phase_detail::sectionGeometry(options, n_)) {
        variables_ = n_ + (autonomous_ ? 1 : 0);
        scales_.reserve(variables_ + 1);
        for (std::size_t j = 0; j < n_; ++j)
            scales_.push_back(options.upper[j] - options.lower[j]);
        if (autonomous_)
            scales_.push_back(std::max(bifurcation.seedPeriod, 10 * options.outputStep));
        scales_.push_back(bifurcation.maximum - bifurcation.minimum);
        sectionScale_ = 0;
        for (double value : section_.scaledNormal)
            sectionScale_ += value * value;
        sectionScale_ = std::sqrt(sectionScale_);
    }

    std::size_t variables() const { return variables_; }
    std::size_t dimension() const { return n_; }
    bool autonomousPeriodic() const { return autonomous_; }
    const PhaseState &scales() const { return scales_; }
    double tolerance() const { return std::max(1e-10, 10 * options_.residualTolerance); }

    Evaluation evaluate(const PhaseState &y, const ExecutionContext *execution,
                        bool derivatives = true) const {
        const double parameter = y.back();
        auto system = factory_(parameter);
        if (!system || system->descriptor().stateDimension != n_)
            throw std::runtime_error("Bifurcation family changed state dimension");
        PhaseState state(y.begin(), y.begin() + n_);
        Evaluation result;
        if (!periodic_) {
            result.residual = phase_detail::field(*system, state);
            result.spectrum = dynsyslab::numerics::eigensystem(
                                  system->jacobian(0, state.data()), n_, false)
                                  .eigenvalues;
            result.stability = classifySpectrum(result.spectrum, false,
                                                options_.eigenvalueTolerance);
            if (derivatives) {
                result.jacobian = system->jacobian(0, state.data());
                const auto parameterColumn = parameterDerivative(y, execution, false);
                appendColumn(result.jacobian, parameterColumn, n_, n_);
            }
            return result;
        }

        double duration = autonomous_ ? y[n_]
                                      : phaseSpaceForcingPeriod(*system, options_.forcingPeriod) *
                                            options_.forcingPeriods;
        if (!std::isfinite(duration) || duration <= 2 * options_.outputStep)
            throw std::runtime_error("Periodic branch proposed an invalid return time");
        auto run = phase_detail::flow(*system, state, duration, options_, true, false, execution,
                                      autonomous_ ? 0 : options_.forcingSectionTime);
        if (!run.complete)
            throw std::runtime_error("Periodic branch return failed: " + run.reason);
        result.residual = subtract(run.state, state);
        result.period = duration;
        if (autonomous_)
            result.residual.push_back(phase_detail::sectionValue(state, section_) / sectionScale_);
        if (autonomous_) {
            const auto map = phase_detail::poincareFromTransition(*system, state, run, options_,
                                                                  section_);
            if (map.status != "complete")
                throw std::runtime_error(map.reason);
            result.spectrum = map.multipliers;
        } else {
            result.spectrum = dynsyslab::numerics::eigensystem(run.monodromy, n_, false).eigenvalues;
        }
        result.stability = classifySpectrum(result.spectrum, true,
                                            options_.eigenvalueTolerance);
        if (!derivatives)
            return result;

        result.jacobian.assign(variables_ * (variables_ + 1), 0);
        for (std::size_t row = 0; row < n_; ++row) {
            for (std::size_t column = 0; column < n_; ++column)
                result.jacobian[row * (variables_ + 1) + column] =
                    run.monodromy[row * n_ + column] - (row == column ? 1. : 0.);
            if (autonomous_)
                result.jacobian[row * (variables_ + 1) + n_] =
                    phase_detail::field(*system, run.state)[row];
        }
        if (autonomous_)
            for (std::size_t column = 0; column < n_; ++column)
                result.jacobian[n_ * (variables_ + 1) + column] =
                    section_.physicalNormal[column] / sectionScale_;
        const auto parameterColumn = parameterDerivative(y, execution, true);
        for (std::size_t row = 0; row < variables_; ++row)
            result.jacobian[row * (variables_ + 1) + variables_] = parameterColumn[row];
        return result;
    }

    PhaseState scaledResidual(const Evaluation &evaluation) const {
        auto result = evaluation.residual;
        for (std::size_t row = 0; row < n_; ++row)
            result[row] /= scales_[row];
        return result;
    }

    std::vector<double> scaledJacobian(const Evaluation &evaluation) const {
        std::vector<double> result(variables_ * (variables_ + 1));
        for (std::size_t row = 0; row < variables_; ++row) {
            const double rowScale = row < n_ ? scales_[row] : 1.;
            for (std::size_t column = 0; column < variables_ + 1; ++column)
                result[row * (variables_ + 1) + column] =
                    evaluation.jacobian[row * (variables_ + 1) + column] *
                    scales_[column] / rowScale;
        }
        return result;
    }

    bool admissible(const PhaseState &y) const {
        if (y.size() != variables_ + 1 || y.back() < bifurcation_.minimum ||
            y.back() > bifurcation_.maximum)
            return false;
        PhaseState state(y.begin(), y.begin() + n_);
        if (!phase_detail::inside(state, options_))
            return false;
        return !autonomous_ ||
               (y[n_] > 2 * options_.outputStep && y[n_] <= 2 * options_.horizon);
    }

  private:
    PhaseState parameterDerivative(const PhaseState &y, const ExecutionContext *execution,
                                   bool periodic) const {
        const double parameter = y.back();
        const double range = bifurcation_.maximum - bifurcation_.minimum;
        const double h = std::max(std::cbrt(std::numeric_limits<double>::epsilon()) *
                                      std::max({1., std::abs(parameter), range}),
                                  1e-7 * range);
        const double minus = std::max(bifurcation_.minimum, parameter - h);
        const double plus = std::min(bifurcation_.maximum, parameter + h);
        if (!(plus > minus))
            throw std::runtime_error("Bifurcation parameter derivative has no finite interval");
        PhaseState state(y.begin(), y.begin() + n_);
        auto value = [&](double mu) {
            auto system = factory_(mu);
            if (!periodic)
                return phase_detail::field(*system, state);
            const double duration = autonomous_
                                        ? y[n_]
                                        : phaseSpaceForcingPeriod(*system, options_.forcingPeriod) *
                                              options_.forcingPeriods;
            auto run = phase_detail::flow(*system, state, duration, options_, false, false,
                                          execution,
                                          autonomous_ ? 0 : options_.forcingSectionTime);
            if (!run.complete)
                throw std::runtime_error("Parameter-difference return failed: " + run.reason);
            auto endpoint = run.state;
            if (autonomous_)
                endpoint.push_back(phase_detail::sectionValue(state, section_) / sectionScale_);
            return endpoint;
        };
        auto left = value(minus), right = value(plus);
        for (std::size_t row = 0; row < left.size(); ++row)
            left[row] = (right[row] - left[row]) / (plus - minus);
        return left;
    }

    static void appendColumn(std::vector<double> &matrix, const PhaseState &column,
                             std::size_t rows, std::size_t existingColumns) {
        std::vector<double> expanded(rows * (existingColumns + 1));
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t col = 0; col < existingColumns; ++col)
                expanded[row * (existingColumns + 1) + col] =
                    matrix[row * existingColumns + col];
            expanded[row * (existingColumns + 1) + existingColumns] = column[row];
        }
        matrix = std::move(expanded);
    }

    const PhaseSystemFactory &factory_;
    const PhaseSpaceOptions &options_;
    const PhaseBifurcationOptions &bifurcation_;
    std::size_t n_{}, variables_{};
    bool periodic_{}, autonomous_{};
    phase_detail::SectionGeometry section_;
    PhaseState scales_;
    double sectionScale_{};
};

PhaseState normalizedSecant(const PhaseState &left, const PhaseState &right,
                            const PhaseState &scales, const PhaseState *orientation = nullptr) {
    PhaseState tangent(left.size());
    double length = 0;
    for (std::size_t j = 0; j < left.size(); ++j) {
        tangent[j] = (right[j] - left[j]) / scales[j];
        length += tangent[j] * tangent[j];
    }
    length = std::sqrt(length);
    if (!std::isfinite(length) || length <= 1e-14)
        throw std::runtime_error("Bifurcation continuation produced a zero secant");
    for (auto &value : tangent)
        value /= length;
    if (orientation) {
        double agreement = 0;
        for (std::size_t j = 0; j < tangent.size(); ++j)
            agreement += tangent[j] * (*orientation)[j];
        if (agreement < 0)
            for (auto &value : tangent)
                value = -value;
    }
    return tangent;
}

PhaseState correct(const BranchProblem &problem, PhaseState y, const PhaseState &prediction,
                   const PhaseState &tangent, const PhaseBifurcationOptions &bifurcation,
                   const ExecutionContext *execution, double &correctionNorm,
                   std::size_t &iterations) {
    const auto m = problem.variables();
    const auto &scales = problem.scales();
    correctionNorm = 0;
    for (iterations = 0; iterations < bifurcation.maximumCorrectorIterations; ++iterations) {
        if (execution)
            execution->checkpoint();
        const auto evaluation = problem.evaluate(y, execution);
        auto residual = problem.scaledResidual(evaluation);
        double arc = 0;
        for (std::size_t j = 0; j < y.size(); ++j)
            arc += tangent[j] * (y[j] - prediction[j]) / scales[j];
        const double cost = std::max(maximumAbs(residual), std::abs(arc));
        if (cost <= problem.tolerance())
            return y;
        auto matrix = problem.scaledJacobian(evaluation);
        matrix.resize((m + 1) * (m + 1));
        for (std::size_t column = 0; column < m + 1; ++column)
            matrix[m * (m + 1) + column] = tangent[column];
        residual.push_back(arc);
        for (auto &value : residual)
            value = -value;
        const auto delta =
            dynsyslab::numerics::solveLeastSquares(matrix, m + 1, m + 1, residual, 1);
        bool accepted = false;
        for (double alpha = 1; alpha >= 1. / 1024; alpha *= .5) {
            auto candidate = y;
            for (std::size_t j = 0; j < candidate.size(); ++j)
                candidate[j] += alpha * delta[j] * scales[j];
            if (!problem.admissible(candidate))
                continue;
            try {
                const auto test = problem.scaledResidual(problem.evaluate(candidate, execution, false));
                double testArc = 0;
                for (std::size_t j = 0; j < candidate.size(); ++j)
                    testArc += tangent[j] * (candidate[j] - prediction[j]) / scales[j];
                if (std::max(maximumAbs(test), std::abs(testArc)) < cost) {
                    y = std::move(candidate);
                    correctionNorm = std::max(correctionNorm, alpha * maximumAbs(delta));
                    accepted = true;
                    break;
                }
            } catch (const std::exception &) {
            }
        }
        if (!accepted)
            break;
    }
    throw std::runtime_error("Pseudo-arclength corrector did not converge");
}

PhaseBifurcationPoint point(const BranchProblem &problem, const PhaseState &y,
                            const ExecutionContext *execution) {
    const auto evaluation = problem.evaluate(y, execution, false);
    PhaseBifurcationPoint result;
    result.parameter = y.back();
    result.state.assign(y.begin(), y.begin() + problem.dimension());
    result.period = evaluation.period;
    result.residual = maximumAbs(problem.scaledResidual(evaluation));
    result.spectrum = evaluation.spectrum;
    result.stability = evaluation.stability;
    return result;
}

} // namespace

PhaseBifurcationBranch continuePhaseBifurcation(const PhaseSystemFactory &factory,
                                                const PhaseSpaceOptions &options,
                                                const PhaseBifurcationOptions &bifurcation,
                                                const ExecutionContext *execution) {
    if (!factory || (bifurcation.branch != "equilibrium" &&
                     bifurcation.branch != "periodic-orbit"))
        throw std::invalid_argument("Bifurcation continuation requires a system family and branch kind");
    if (!(bifurcation.minimum < bifurcation.maximum) ||
        bifurcation.initial < bifurcation.minimum || bifurcation.initial > bifurcation.maximum ||
        !std::isfinite(bifurcation.step) || bifurcation.step == 0 ||
        bifurcation.minimumStep <= 0 || bifurcation.maximumStep < bifurcation.minimumStep ||
        bifurcation.maximumPoints < 2)
        throw std::invalid_argument("Invalid bifurcation interval, step, or point budget");
    BranchProblem problem(factory, options, bifurcation);
    PhaseBifurcationBranch branch;
    branch.kind = bifurcation.branch;
    branch.parameter = bifurcation.parameter;
    branch.status = "truncated";
    branch.diagnostics.push_back(
        "Pseudo-arclength continuation is finite numerical evidence; detected stability changes are candidates, not proofs of bifurcation type");

    PhaseState first = bifurcation.seedState;
    if (first.empty())
        first = factory(bifurcation.initial)->initialState();
    if (first.size() != options.lower.size())
        throw std::invalid_argument("Bifurcation seed state has the wrong dimension");
    if (bifurcation.branch == "periodic-orbit" && options.dynamics == "autonomous-zero-input")
        first.push_back(bifurcation.seedPeriod);
    first.push_back(bifurcation.initial);
    PhaseState initialTangent(first.size(), 0);
    initialTangent.back() = bifurcation.step > 0 ? 1 : -1;
    auto prediction = first;
    prediction.back() += bifurcation.step;
    if (!problem.admissible(prediction))
        throw std::invalid_argument("The first bifurcation continuation step leaves the declared interval");
    double correctionNorm = 0;
    std::size_t iterations = 0;
    first = correct(problem, first, first, initialTangent, bifurcation, execution,
                    correctionNorm, iterations);
    auto second = correct(problem, prediction, prediction, initialTangent, bifurcation, execution,
                          correctionNorm, iterations);
    branch.points.push_back(point(problem, first, execution));
    branch.points.push_back(point(problem, second, execution));
    auto tangent = normalizedSecant(first, second, problem.scales());
    double step = std::clamp(std::abs(bifurcation.step), bifurcation.minimumStep,
                             bifurcation.maximumStep);
    while (branch.points.size() < bifurcation.maximumPoints) {
        prediction = second;
        for (std::size_t j = 0; j < prediction.size(); ++j)
            prediction[j] += step * tangent[j] * problem.scales()[j];
        if (!problem.admissible(prediction)) {
            branch.status = "complete";
            branch.reason = "Continuation reached the declared parameter or state boundary";
            break;
        }
        try {
            auto next = correct(problem, prediction, prediction, tangent, bifurcation, execution,
                                correctionNorm, iterations);
            auto nextPoint = point(problem, next, execution);
            nextPoint.arclengthStep = step;
            nextPoint.correctionNorm = correctionNorm;
            branch.points.push_back(std::move(nextPoint));
            auto nextTangent = normalizedSecant(second, next, problem.scales(), &tangent);
            first = std::move(second);
            second = std::move(next);
            tangent = std::move(nextTangent);
            if (iterations <= 4)
                step = std::min(bifurcation.maximumStep, step * 1.25);
        } catch (const ExecutionCancelled &) {
            throw;
        } catch (const std::exception &error) {
            step *= .5;
            if (step < bifurcation.minimumStep) {
                branch.reason = error.what();
                break;
            }
        }
    }
    if (branch.points.size() == bifurcation.maximumPoints && branch.reason.empty())
        branch.reason = "Maximum continuation point budget reached";
    phase_detail::detectBifurcationEvents(factory, options, bifurcation, branch);
    return branch;
}

namespace phase_detail {
void findBifurcationBranches(const ODESystem &system, PhaseSpaceResult &result,
                             const ExecutionContext *execution) {
    const auto bifurcation = result.options.bifurcation;
    if (!bifurcation.enabled)
        return;
    auto base = system.configuration();
    const auto parameter = bifurcation.parameter;
    const bool ordinary = base.parameters.count(parameter) != 0;
    const bool driver = parameter == "driver.amplitude" || parameter == "driver.frequency_hz" ||
                        parameter == "driver.phase_radians";
    if (!ordinary && !driver)
        throw std::invalid_argument("Unknown bifurcation parameter: " + parameter);
    PhaseSystemFactory factory = [base, parameter, ordinary](double value) mutable {
        auto configuration = base;
        if (ordinary)
            configuration.parameters[parameter] = value;
        else if (parameter == "driver.amplitude")
            configuration.driver.amplitude = value;
        else if (parameter == "driver.frequency_hz")
            configuration.driver.frequencyHz = value;
        else
            configuration.driver.phaseRadians = value;
        return createSystem(configuration);
    };
    result.bifurcationBranches.push_back(
        continuePhaseBifurcation(factory, result.options, bifurcation, execution));
}
} // namespace phase_detail
} // namespace dynsyslab::core
