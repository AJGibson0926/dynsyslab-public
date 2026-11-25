#include <dynsyslab/core/FunctionRepresentation.hpp>
#include <dynsyslab/core/Observable.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {

Observable::Observable(ObservableDescriptor descriptor, std::size_t stateDimension,
                       std::size_t outputDimension, Evaluator evaluator,
                       Jacobian analyticJacobian)
    : descriptor_(std::move(descriptor)), stateDimension_(stateDimension),
      outputDimension_(outputDimension), evaluator_(std::move(evaluator)),
      analyticJacobian_(std::move(analyticJacobian)) {
    if (descriptor_.name.empty() || !stateDimension_ || !outputDimension_ || !evaluator_)
        throw std::invalid_argument("observable requires a name, dimensions, and evaluator");
    if (descriptor_.label.empty()) descriptor_.label = descriptor_.name;
    if (!descriptor_.outputNames.empty()
        && descriptor_.outputNames.size() != outputDimension_)
        throw std::invalid_argument("observable output names do not match its dimension");
}

const ObservableDescriptor& Observable::descriptor() const noexcept { return descriptor_; }
Observable Observable::withDecomposition(ObservableDecompositionReference reference) const {
    if(reference.artifact.empty()||reference.modelIdentity.empty()||
       (reference.kind!="reconstruction"&&reference.kind!="eigenfunction"&&reference.kind!="spectral-measure"&&reference.kind!="wave-packet"))
        throw std::invalid_argument("Observable decomposition reference needs artifact, identity and a supported kind");
    auto result=*this;result.descriptor_.decompositions.push_back(std::move(reference));return result;
}
std::size_t Observable::stateDimension() const noexcept { return stateDimension_; }
std::size_t Observable::outputDimension() const noexcept { return outputDimension_; }
std::optional<std::size_t> Observable::stateComponentIndex(std::size_t output) const {
    if (output >= outputDimension_) throw std::out_of_range("observable output component");
    return output < stateComponents_.size() ? stateComponents_[output] : std::nullopt;
}
bool Observable::hasAnalyticJacobian() const noexcept { return bool(analyticJacobian_); }

TimeSeries Observable::evaluate(const TimeSeries& series) const {
    const auto* window = representation_ ? representation_->coordinates().window.get() : nullptr;
    if (!window && series.componentCount() != stateDimension_)
        throw std::invalid_argument(descriptor_.name + ": the series has " + std::to_string(series.componentCount())
                                    + " components, the observable takes " + std::to_string(stateDimension_));
    const bool complex = descriptor_.complexValued;
    std::vector<double> time, values;
    std::size_t width = 0;
    for (std::size_t k = 0; k < series.sampleCount(); ++k) {
        Values g;
        if (window) {
            if (!window->complete(series.sampleCount(), k)) continue;
            const auto q = window->window_at(series, k);
            g = evaluate(q.data());
        } else {
            std::vector<double> x(stateDimension_);
            for (std::size_t i = 0; i < stateDimension_; ++i) x[i] = series.value(k, i);
            g = evaluate(x.data());
        }
        width = complex ? 2 * g.size() : g.size();
        time.push_back(series.time(k));
        for (const auto& value : g) {
            values.push_back(value.real());
            if (complex) values.push_back(value.imag());
        }
    }
    if (time.empty()) throw std::invalid_argument(descriptor_.name + ": the series is shorter than one window");
    return TimeSeries(std::move(time), std::move(values), width);
}

Observable::Values Observable::evaluate(const double* state) const {
    auto result = evaluator_(state);
    if (result.size() != outputDimension_)
        throw std::runtime_error("observable evaluator returned the wrong dimension");
    for (const auto value : result)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::runtime_error("observable evaluation is non-finite");
    return result;
}

Observable::Values Observable::jacobian(const double* state,
                                        ObservableJacobianBackend backend,
                                        double relativeStep) const {
    if (backend == ObservableJacobianBackend::automatic)
        backend = analyticJacobian_ ? ObservableJacobianBackend::analytic
                                    : ObservableJacobianBackend::finiteDifference;
    if (backend == ObservableJacobianBackend::analytic) {
        if (!analyticJacobian_)
            throw std::logic_error("observable has no analytic Jacobian");
        auto result = analyticJacobian_(state);
        if (result.size() != outputDimension_ * stateDimension_)
            throw std::runtime_error("observable Jacobian returned the wrong shape");
        return result;
    }
    if (!(relativeStep > 0.0) || !std::isfinite(relativeStep))
        throw std::invalid_argument("observable finite-difference step must be positive");
    Values result(outputDimension_ * stateDimension_);
    std::vector<double> plus(state, state + stateDimension_), minus(plus);
    for (std::size_t axis = 0; axis < stateDimension_; ++axis) {
        const double step = relativeStep * std::max(1.0, std::abs(state[axis]));
        plus[axis] += step; minus[axis] -= step;
        const auto high = evaluate(plus.data()), low = evaluate(minus.data());
        for (std::size_t output = 0; output < outputDimension_; ++output)
            result[output * stateDimension_ + axis] = (high[output] - low[output]) / (2.0 * step);
        plus[axis] = state[axis]; minus[axis] = state[axis];
    }
    return result;
}

std::vector<double> Observable::evaluateReal(const double* state) const {
    const auto values = evaluate(state); std::vector<double> result(values.size());
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (std::abs(values[index].imag()) > 1.0e-12)
            throw std::invalid_argument("a complex observable cannot enter a real-valued method");
        result[index] = values[index].real();
    }
    return result;
}

std::vector<double> Observable::evaluateRealBatch(const std::vector<double>& states,
                                                  std::size_t samples) const {
    if (states.size() != samples * stateDimension_)
        throw std::invalid_argument("observable batch shape mismatch");
    std::vector<double> result(samples * outputDimension_);
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const auto value = evaluateReal(states.data() + sample * stateDimension_);
        std::copy(value.begin(), value.end(), result.begin() + sample * outputDimension_);
    }
    return result;
}

Observable Observable::stateComponent(std::size_t dimension, std::size_t index,
                                      std::string name, std::string label,
                                      std::string coordinates) {
    if (index >= dimension) throw std::out_of_range("observable state component");
    ObservableDescriptor descriptor{std::move(name), std::move(label), "",
                                    std::move(coordinates), {}, false};
    Observable result{std::move(descriptor), dimension, 1,
        [index](const double* state) { return Values{state[index]}; },
        [dimension, index](const double*) {
            Values result(dimension); result[index] = 1.0; return result;
        }};
    result.stateComponents_ = {index};
    return result;
}

Observable Observable::identity(std::size_t dimension, std::vector<std::string> names,
                                std::vector<std::string> labels, std::string coordinates) {
    if (names.size() != dimension || labels.size() != dimension)
        throw std::invalid_argument("identity observable names do not match state dimension");
    ObservableDescriptor descriptor{"state", "Full state", "", std::move(coordinates),
                                    std::move(names), false};
    Observable result{std::move(descriptor), dimension, dimension,
        [dimension](const double* state) {
            Values result(dimension); for (std::size_t i=0;i<dimension;++i) result[i]=state[i];
            return result;
        },
        [dimension](const double*) {
            Values result(dimension*dimension);
            for (std::size_t i=0;i<dimension;++i) result[i*dimension+i]=1.0;
            return result;
        }};
    for (std::size_t i = 0; i < dimension; ++i) {
        result.stateComponents_.push_back(i);
    }
    return result;
}

Observable Observable::stack(std::string name, std::string label,
                             const std::vector<Observable>& observables) {
    if (observables.empty()) throw std::invalid_argument("observable stack is empty");
    const std::size_t stateDimension = observables.front().stateDimension();
    std::size_t outputDimension = 0; bool complex = false, analytic = true;
    std::vector<std::string> outputNames;
    for (const auto& item : observables) {
        if (item.stateDimension() != stateDimension)
            throw std::invalid_argument("stacked observables have different state dimensions");
        outputDimension += item.outputDimension(); complex |= item.descriptor().complexValued;
        analytic &= item.hasAnalyticJacobian();
        if (item.outputDimension() == 1) outputNames.push_back(item.descriptor().name);
        else outputNames.insert(outputNames.end(), item.descriptor().outputNames.begin(),
                                item.descriptor().outputNames.end());
    }
    ObservableDescriptor descriptor{std::move(name), std::move(label), "",
        observables.front().descriptor().stateCoordinateSystem, std::move(outputNames), complex};
    Evaluator evaluate = [observables](const double* state) {
        Values result;
        for (const auto& item : observables) {
            const auto value = item.evaluate(state); result.insert(result.end(), value.begin(), value.end());
        }
        return result;
    };
    Jacobian jacobian;
    if (analytic) jacobian = [observables](const double* state) {
        Values result;
        for (const auto& item : observables) {
            const auto value = item.jacobian(state, ObservableJacobianBackend::analytic);
            result.insert(result.end(), value.begin(), value.end());
        }
        return result;
    };
    std::vector<AnalyticExpression> formulas;
    bool closed = !complex;
    for (const auto& item : observables) {
        const auto own = item.expressions();
        closed = closed && !own.empty();
        formulas.insert(formulas.end(), own.begin(), own.end());
    }
    // A stack of formulas is a formula observable (an expansion); otherwise it stacks evaluators.
    Observable result = closed ? fromExpressions(descriptor, std::move(formulas))
                               : Observable{descriptor, stateDimension, outputDimension, std::move(evaluate), std::move(jacobian)};
    for (const auto& item : observables)
        for (std::size_t i = 0; i < item.outputDimension(); ++i)
            result.stateComponents_.push_back(item.descriptor().stateCoordinateSystem == result.descriptor().stateCoordinateSystem
                ? item.stateComponentIndex(i) : std::nullopt);
    return result;
}

std::vector<double> Observable::pushForwardRealBatch(
    const std::vector<double>& states, const std::vector<double>& velocities,
    std::size_t samples, ObservableJacobianBackend backend) const {
    if (states.size()!=samples*stateDimension_ || velocities.size()!=states.size())
        throw std::invalid_argument("observable push-forward batch shape mismatch");
    std::vector<double> result(samples*outputDimension_);
    for (std::size_t row=0;row<samples;++row) {
        const auto gradient=jacobian(states.data()+row*stateDimension_,backend);
        for (std::size_t output=0;output<outputDimension_;++output) {
            std::complex<double> value{};
            for (std::size_t axis=0;axis<stateDimension_;++axis)
                value+=gradient[output*stateDimension_+axis]*velocities[row*stateDimension_+axis];
            if (std::abs(value.imag())>1.0e-12)
                throw std::invalid_argument("real observable push-forward has a complex component");
            result[row*outputDimension_+output]=value.real();
        }
    }
    return result;
}

ObservableCatalog::ObservableCatalog(std::vector<Observable> observables,
    std::map<std::string, std::string> aliases,
    std::map<std::string, std::vector<std::string>> coordinatePresets, ObservableDerivativeLinks derivatives)
    : derivatives_(std::move(derivatives)), aliases_(std::move(aliases)), coordinatePresets_(std::move(coordinatePresets)) {
    for (auto& item : observables) {
        const auto name = item.descriptor().name;
        if (!observables_.emplace(name, std::move(item)).second)
            throw std::invalid_argument("duplicate observable name: " + name);
    }
    for (const auto& entry:derivatives_) {
        const auto from=observables_.find(entry.first), to=observables_.find(entry.second.target);
        if (from==observables_.end() || to==observables_.end() ||
            from->second.outputDimension()!=to->second.outputDimension() ||
            from->second.stateDimension()!=to->second.stateDimension() ||
            from->second.descriptor().stateCoordinateSystem!=to->second.descriptor().stateCoordinateSystem)
            throw std::invalid_argument("invalid observable time-derivative link: "+entry.first);
    }
    if (observables_.empty()) throw std::invalid_argument("observable catalog is empty");
    for (const auto& alias : aliases_)
        if (observables_.find(alias.second) == observables_.end())
            throw std::invalid_argument("observable alias refers to an unknown name");
    for (const auto& preset : coordinatePresets_)
        for (const auto& name : preset.second)
            if (!contains(name)) throw std::invalid_argument("observable coordinate preset is invalid");
}

bool ObservableCatalog::contains(const std::string& name) const {
    const auto alias = aliases_.find(name);
    return observables_.find(alias == aliases_.end() ? name : alias->second) != observables_.end();
}
const Observable& ObservableCatalog::at(const std::string& name) const {
    const auto alias = aliases_.find(name);
    return observables_.at(alias == aliases_.end() ? name : alias->second);
}
const ObservableTimeDerivative* ObservableCatalog::timeDerivative(const std::string& name) const {
    const auto it=derivatives_.find(at(name).descriptor().name);
    return it==derivatives_.end()?nullptr:&it->second;
}
Observable ObservableCatalog::coordinates(const std::vector<std::string>& names,
                                          std::string name) const {
    std::vector<Observable> members;
    for (const auto& item : names) members.push_back(at(item));
    if (name.empty()) {
        name = "coordinates[";
        for (std::size_t i=0;i<names.size();++i) { if (i) name += ','; name += names[i]; }
        name += ']';
    }
    return Observable::stack(name, name, members);
}
Observable ObservableCatalog::coordinatePreset(const std::string& name) const {
    return coordinates(coordinatePresets_.at(name), name);
}
std::vector<std::string> ObservableCatalog::names() const {
    std::vector<std::string> result; result.reserve(observables_.size());
    for (const auto& item : observables_) result.push_back(item.first); return result;
}
const std::map<std::string, std::string>& ObservableCatalog::aliases() const noexcept {
    return aliases_;
}
const std::map<std::string, std::vector<std::string>>&
ObservableCatalog::coordinatePresets() const noexcept { return coordinatePresets_; }

}  // namespace dynsyslab::core
