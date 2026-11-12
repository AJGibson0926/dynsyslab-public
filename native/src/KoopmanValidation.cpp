#include <dynsyslab/core/KoopmanValidation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

KoopmanErrorStatistics statistics(
    const std::vector<std::complex<double>>& actual,
    const std::vector<std::complex<double>>& predicted) {
    if (actual.empty() || actual.size() != predicted.size())
        throw std::invalid_argument("Koopman error statistics require equal nonempty samples");
    double errorSquared = 0.0, actualSquared = 0.0, maximum = 0.0;
    double actualMaximum = 0.0;
    std::vector<double> magnitudes(actual.size()), relative(actual.size());
    for (std::size_t index = 0; index < actual.size(); ++index) {
        magnitudes[index] = std::abs(actual[index] - predicted[index]);
        errorSquared += magnitudes[index] * magnitudes[index];
        actualSquared += std::norm(actual[index]);
        maximum = std::max(maximum, magnitudes[index]);
        actualMaximum = std::max(actualMaximum, std::abs(actual[index]));
    }
    const double tiny = std::numeric_limits<double>::min();
    const double floor = std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(actualMaximum, tiny);
    for (std::size_t index = 0; index < actual.size(); ++index)
        relative[index] = magnitudes[index] / std::max(std::abs(actual[index]), floor);
    std::sort(relative.begin(), relative.end());
    const auto quantile = [&](double fraction) {
        const double position = fraction * static_cast<double>(relative.size() - 1);
        const auto lower = static_cast<std::size_t>(std::floor(position));
        const auto upper = static_cast<std::size_t>(std::ceil(position));
        const double weight = position - static_cast<double>(lower);
        return relative[lower] * (1.0 - weight) + relative[upper] * weight;
    };
    return {actual.size(), std::sqrt(errorSquared / actual.size()), maximum,
        std::sqrt(errorSquared) / std::max(std::sqrt(actualSquared), tiny),
        quantile(0.5), quantile(0.95)};
}

std::complex<double> interpolate(const std::vector<double>& time,
    const std::vector<std::complex<double>>& values, double target) {
    const auto upper = std::lower_bound(time.begin(), time.end(), target);
    if (upper == time.begin()) return values.front();
    if (upper == time.end()) return values.back();
    const auto right = static_cast<std::size_t>(upper - time.begin());
    const auto left = right - 1;
    const double weight = (target - time[left]) / (time[right] - time[left]);
    return values[left] * (1.0 - weight) + values[right] * weight;
}

std::vector<double> defaultHorizons(const std::vector<double>& time) {
    double minimumStep = std::numeric_limits<double>::infinity();
    for (std::size_t index = 1; index < time.size(); ++index)
        minimumStep = std::min(minimumStep, time[index] - time[index - 1]);
    const double duration = time.back() - time.front();
    std::vector<double> candidates{10.0 * minimumStep, 100.0 * minimumStep,
        0.01 * duration, 0.05 * duration, 0.1 * duration, 0.25 * duration,
        0.5 * duration};
    candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
        [&](double value) { return !(value > minimumStep && value <= duration); }),
        candidates.end());
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end(),
        [](double left, double right) {
            return std::abs(left - right) <= 1.0e-12 * std::max({1.0, left, right});
        }), candidates.end());
    return candidates;
}

std::vector<std::complex<double>> derivative(
    const std::vector<double>& time,
    const std::vector<std::complex<double>>& values) {
    std::vector<std::complex<double>> result(values.size());
    result.front() = (values[1] - values[0]) / (time[1] - time[0]);
    result.back() = (values.back() - values[values.size() - 2])
        / (time.back() - time[time.size() - 2]);
    for (std::size_t index = 1; index + 1 < values.size(); ++index)
        result[index] = (values[index + 1] - values[index - 1])
            / (time[index + 1] - time[index - 1]);
    return result;
}

}  // namespace

EigenfunctionValidationResult validateEigenfunctions(
    const KoopmanModel& model, const std::vector<double>& time,
    const std::vector<double>& states, std::size_t samples,
    std::size_t stateDimension, EigenfunctionValidationInput input, std::vector<double> horizons,
    const std::vector<double>& trainingMinimum,
    const std::vector<double>& trainingMaximum) {
    const bool coordinates=input==EigenfunctionValidationInput::coordinates;
    const std::size_t expectedDimension = coordinates
        ? model.coordinateDimension() : model.stateDimension();
    if (samples < 3 || time.size() != samples
        || states.size() != samples * stateDimension
        || stateDimension != expectedDimension)
        throw std::invalid_argument("Koopman validation trajectory dimensions are incompatible");
    for(double value:time)if(!std::isfinite(value))throw std::invalid_argument("Validation time must be finite");
    for(double value:states)if(!std::isfinite(value))throw std::invalid_argument("Validation values must be finite");
    for (std::size_t index = 1; index < samples; ++index)
        if (!(time[index] > time[index - 1]))
            throw std::invalid_argument("Koopman validation time must increase strictly");
    if (horizons.empty()) {
        horizons = defaultHorizons(time);
        if(!model.spectrum().hasContinuousRates()) {
            const double dt=model.spectrum().sampleInterval();
            for(auto& h:horizons)h=std::max(1.,std::round(h/dt))*dt;
        }
    }
    for (const double value : horizons)
        if (!(value > 0.0) || !std::isfinite(value))
            throw std::invalid_argument("Koopman validation horizons must be finite and positive");
    std::sort(horizons.begin(), horizons.end());
    horizons.erase(std::unique(horizons.begin(), horizons.end()), horizons.end());

    const std::size_t modes = model.liftedDimension();
    std::vector<std::vector<std::complex<double>>> values(
        modes, std::vector<std::complex<double>>(samples));
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const auto lifted = coordinates
            ? model.evaluateCoordinates(states.data() + sample * stateDimension)
            : model.evaluate(states.data() + sample * stateDimension);
        for (std::size_t mode = 0; mode < modes; ++mode)
            values[mode][sample] = lifted[mode];
    }

    KoopmanValidationResult result;
    result.samples = samples; result.stateDimension = stateDimension;
    result.inputKind=input;
    result.trainingDomainAvailable = trainingMinimum.size() == stateDimension
        && trainingMaximum.size() == stateDimension;
    if (result.trainingDomainAvailable) {
        for (std::size_t sample = 0; sample < samples; ++sample) {
            bool outside = false;
            for (std::size_t axis = 0; axis < stateDimension; ++axis) {
                const double width = std::max(trainingMaximum[axis]
                    - trainingMinimum[axis], std::numeric_limits<double>::epsilon());
                const double value = states[sample * stateDimension + axis];
                const double excursion = std::max({
                    (trainingMinimum[axis] - value) / width,
                    (value - trainingMaximum[axis]) / width, 0.0});
                outside = outside || excursion > 0.0;
                result.maximumRelativeExcursion =
                    std::max(result.maximumRelativeExcursion, excursion);
            }
            if (outside) ++result.outsideTrainingSamples;
        }
        result.outsideTrainingFraction = static_cast<double>(
            result.outsideTrainingSamples) / samples;
    }

    for (std::size_t mode = 0; mode < modes; ++mode) {
        const auto rate=model.spectrum().continuousRate(mode);
        const auto lambda = rate.value_or(std::complex<double>{});
        std::vector<std::complex<double>> stepActual, stepPredicted;
        std::vector<std::complex<double>> anchoredActual, anchoredPredicted;
        stepActual.reserve(samples - 1); stepPredicted.reserve(samples - 1);
        anchoredActual.reserve(samples - 1); anchoredPredicted.reserve(samples - 1);
        for (std::size_t index = 1; index < samples; ++index) {
            stepActual.push_back(values[mode][index]);
            stepPredicted.push_back(model.spectrum().transition(mode,time[index] - time[index - 1])
                * values[mode][index - 1]);
            anchoredActual.push_back(values[mode][index]);
            anchoredPredicted.push_back(model.spectrum().transition(mode,time[index] - time.front())
                * values[mode].front());
        }
        const auto observedDerivative = derivative(time, values[mode]);
        std::vector<std::complex<double>> generatorPrediction(samples);
        for (std::size_t index = 0; index < samples; ++index)
            generatorPrediction[index] = lambda * values[mode][index];
        KoopmanModeValidation item{lambda, statistics(stepActual, stepPredicted),
            statistics(anchoredActual, anchoredPredicted),
            rate?statistics(observedDerivative, generatorPrediction):KoopmanErrorStatistics{}, {}};
        item.continuousRateAvailable=rate.has_value();item.generatorAvailable=rate.has_value();
        long double incrementSquared=0, residualSquared=0;
        for(std::size_t index=1;index<samples;++index) {
            incrementSquared+=std::norm(values[mode][index]-values[mode][index-1]);
            residualSquared+=std::norm(stepActual[index-1]-stepPredicted[index-1]);
        }
        item.oneStepChangeDefined=incrementSquared>0;
        if(item.oneStepChangeDefined)
            item.oneStepChangeRelativeL2=static_cast<double>(std::sqrt(residualSquared/incrementSquared));
        for (const double horizon : horizons) {
            std::vector<std::complex<double>> actual, predicted;
            for (std::size_t index = 0; index < samples; ++index) {
                const double target = time[index] + horizon;
                if (target > time.back()) break;
                actual.push_back(interpolate(time, values[mode], target));
                predicted.push_back(model.spectrum().transition(mode,horizon) * values[mode][index]);
            }
            if (!actual.empty()) item.horizons.push_back(
                {horizon, statistics(actual, predicted)});
        }
        result.modes.push_back(std::move(item));
    }
    return result;
}

KoopmanValidationResult validateKoopmanSemigroup(const KoopmanModel& model,
    const std::vector<double>& time,const std::vector<double>& states,std::size_t samples,
    std::size_t dimension,std::vector<double> horizons,const std::vector<double>& minimum,
    const std::vector<double>& maximum) {
    return validateEigenfunctions(model,time,states,samples,dimension,
        model.requiresExternalCoordinates()?EigenfunctionValidationInput::coordinates:
        EigenfunctionValidationInput::physicalState,std::move(horizons),minimum,maximum);
}
EigenfunctionValidationResult validateEigenfunctions(const PolynomialEigenfunctionSet& bank,
    const std::vector<double>& time,const std::vector<double>& q,std::vector<double> horizons,
    const std::vector<double>& minimum,const std::vector<double>& maximum) {
    return validateEigenfunctions(KoopmanModel(bank),time,q,time.size(),bank.coordinateDimension(),
        EigenfunctionValidationInput::coordinates,std::move(horizons),minimum,maximum);
}

}  // namespace dynsyslab::core
