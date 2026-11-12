#include <dynsyslab/core/TimeSeries.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

bool FourierOptions::operator==(const FourierOptions& other) const noexcept {
    return oneSided == other.oneSided && removeMean == other.removeMean
        && window == other.window;
}

TimeSeries::TimeSeries(std::vector<double> time, std::vector<double> values,
                       std::size_t componentCount)
    : ownedTime_(std::make_shared<const std::vector<double>>(std::move(time))),
      ownedValues_(std::make_shared<const std::vector<double>>(std::move(values))) {
    initialize(ownedTime_->data(), ownedValues_->data(), ownedTime_->size(), componentCount);
    if (ownedValues_->size() != sampleCount_ * componentCount_)
        throw std::invalid_argument("TimeSeries value shape does not match its time axis");
}

TimeSeries TimeSeries::view(const std::vector<double>& time, const double* values,
                            std::size_t sampleCount, std::size_t componentCount) {
    if (time.size() != sampleCount)
        throw std::invalid_argument("TimeSeries view time and sample counts differ");
    TimeSeries result;
    result.initialize(time.data(), values, sampleCount, componentCount);
    return result;
}

void TimeSeries::initialize(const double* time, const double* values,
                            std::size_t samples, std::size_t components) {
    if (!samples || !components || !time || !values)
        throw std::invalid_argument("TimeSeries requires non-empty time and value arrays");
    timeData_ = time; valuesData_ = values;
    sampleCount_ = samples; componentCount_ = components;
    uniformTimeSteps_ = true; deltaTime_ = 0.0;
    if (!std::isfinite(timeData_[0]))
        throw std::invalid_argument("TimeSeries times must be finite");
    if (samples == 1) return;
    deltaTime_ = timeData_[1] - timeData_[0];
    if (!std::isfinite(timeData_[1]) || !(deltaTime_ > 0.0))
        throw std::invalid_argument("TimeSeries times must be finite and strictly increasing");
    for (std::size_t sample = 2; sample < samples; ++sample) {
        const double step = timeData_[sample] - timeData_[sample - 1];
        if (!std::isfinite(timeData_[sample]) || !(step > 0.0))
            throw std::invalid_argument("TimeSeries times must be finite and strictly increasing");
        const double scale = std::max({1.0, std::abs(deltaTime_), std::abs(step)});
        if (std::abs(step - deltaTime_) > 1.0e-12 * scale) uniformTimeSteps_ = false;
    }
}

double TimeSeries::deltaTime() const {
    if (!uniformTimeSteps_ || sampleCount_ < 2)
        throw std::logic_error("TimeSeries does not have a uniform nonzero time step");
    return deltaTime_;
}

double TimeSeries::time(std::size_t sample) const {
    if (sample >= sampleCount_) throw std::out_of_range("TimeSeries sample index");
    return timeData_[sample];
}

double TimeSeries::value(std::size_t sample, std::size_t component) const {
    if (sample >= sampleCount_ || component >= componentCount_)
        throw std::out_of_range("TimeSeries value index");
    return valuesData_[sample * componentCount_ + component];
}

FourierSidecar TimeSeries::fourier() const { return FourierSidecar(*this); }

}  // namespace dynsyslab::core
