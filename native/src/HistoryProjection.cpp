#include <dynsyslab/core/HistoryObservation.hpp>
#include "HistoryCorrelation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dynsyslab::core {
using detail::CorrelationPlan;
using Vector=Eigen::VectorXd;
HistoryProjection projectHistory(const std::vector<double>& observations,
                                 const HankelLeftBasis& basis) {
    // Validate once before the FFT; an invalid scalar must not contaminate every window.
    const FrozenHistoryTransform checked(basis,0);
    for(double value:observations)if(!std::isfinite(value))
        throw std::invalid_argument("History observations must be finite");
    if (basis.rows < 2 || !basis.rank
        || basis.leftSingularVectors.size() != basis.rows * basis.rank
        || basis.singularValues.size() != basis.rank
        || observations.size() < basis.rows)
        throw std::invalid_argument("invalid history projection inputs");
    CorrelationPlan plan(observations, basis.rows);
    HistoryProjection result;
    result.windows = observations.size() - basis.rows + 1;
    result.rank = basis.rank;
    result.coordinates.resize(result.windows * result.rank);
    for (std::size_t mode = 0; mode < basis.rank; ++mode) {
        if (!(basis.singularValues[mode] > std::numeric_limits<double>::epsilon()))
            throw std::invalid_argument("history basis contains a zero singular value");
        Vector vector(static_cast<Eigen::Index>(basis.rows));
        for (std::size_t delay = 0; delay < basis.rows; ++delay)
            vector[static_cast<Eigen::Index>(delay)]
                = basis.leftSingularVectors[delay * basis.rank + mode];
        const Vector projected = plan.correlate(vector) / basis.singularValues[mode];
        for (std::size_t window = 0; window < result.windows; ++window)
            result.coordinates[window * result.rank + mode]
                = projected[static_cast<Eigen::Index>(window)];
    }
    return result;
}

CausalHistoryProjector::CausalHistoryProjector(HankelLeftBasis basis)
    :CausalHistoryProjector(FrozenHistoryTransform(std::move(basis),0)) {}
CausalHistoryProjector::CausalHistoryProjector(FrozenHistoryTransform transform)
    :transform_(std::move(transform)),ring_(basis().rows) {
    const auto& b=basis();const auto& selected=transform_.retainedCoordinates();
    weights_.reserve(selected.size()*b.rows);
    for(const auto mode:selected)for(std::size_t delay=0;delay<b.rows;++delay)
        weights_.push_back(b.leftSingularVectors[delay*b.rank+mode]);
}

std::optional<std::vector<double>> CausalHistoryProjector::push(double observation) {
    append(observation);
    return projectLatest();
}

void CausalHistoryProjector::append(double observation) {
    if(transform_.sampleInterval()>0)throw std::invalid_argument("Clocked history requires a timestamp");
    appendValue(observation);
}
std::optional<std::vector<double>> CausalHistoryProjector::push(double time,double observation) {
    append(time,observation);return projectLatest();
}
void CausalHistoryProjector::append(double time,double observation) {
    if(!std::isfinite(time))throw std::invalid_argument("History timestamp must be finite");
    const double dt=transform_.sampleInterval();
    if(lastTime_ && (time<=*lastTime_ || (dt>0 && std::abs(time-*lastTime_-dt)>
        1e-8*dt+8*std::numeric_limits<double>::epsilon()*std::max(std::abs(time),std::abs(*lastTime_)))))
        throw std::invalid_argument("History measurement clock changed; reset or explicitly resample the stream");
    appendValue(observation);lastTime_=time;
}
void CausalHistoryProjector::appendValue(double observation) {
    if (!std::isfinite(observation))
        throw std::invalid_argument("causal history observation must be finite");
    ring_[writeIndex_] = observation;
    writeIndex_ = (writeIndex_ + 1) % ring_.size();
    ++samplesSeen_;
}

std::optional<std::vector<double>> CausalHistoryProjector::projectLatest() const {
    if (samplesSeen_ < ring_.size()) return std::nullopt;
    const auto& b=basis();const auto& selected=transform_.retainedCoordinates();
    std::vector<double> result(selected.size());
    // Oldest to newest: the ring's two contiguous segments, in the original summation order.
    const std::size_t tail = b.rows - writeIndex_;
    for (std::size_t axis = 0; axis < selected.size(); ++axis) {
        const double* weight = weights_.data() + axis * b.rows;
        const double* older = ring_.data() + writeIndex_;
        const double* newer = ring_.data();
        double sum = 0.0;
        for (std::size_t j = 0; j < tail; ++j) sum += older[j] * weight[j];
        for (std::size_t j = 0; j < writeIndex_; ++j) sum += newer[j] * weight[tail + j];
        result[axis] = sum / b.singularValues[selected[axis]];
    }
    return result;
}

void CausalHistoryProjector::reset() noexcept {
    std::fill(ring_.begin(), ring_.end(), 0.0);
    writeIndex_ = 0;
    samplesSeen_ = 0;
    lastTime_.reset();
}

} // namespace dynsyslab::core
