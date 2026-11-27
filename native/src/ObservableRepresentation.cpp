// Observables backed by a FunctionRepresentation (expansion or interpolation), and Eigenfunction.
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/FunctionRepresentation.hpp>

#include <algorithm>
#include <stdexcept>

namespace dynsyslab::core {
Observable Observable::fromRepresentation(ObservableDescriptor descriptor, std::shared_ptr<const FunctionRepresentation> rep,
                                          const ObservableCatalog* catalog) {
    if (!rep) throw std::invalid_argument("an observable representation is required");
    if (rep->coordinates().needsCatalog()) {
        if (!catalog) throw std::invalid_argument("the representation domain is catalog observables; pass the system catalog");
        rep = rep->onDomain(rep->coordinates().bound(*catalog));
    }
    // No map (q is supplied) or the identity: the observable is a function of q directly.
    const auto& domain = rep->coordinates();
    const std::shared_ptr<const Observable> map = domain.isPhysicalState() ? nullptr : domain.stateMap;
    const std::size_t outputs = rep->outputCount(), qDim = rep->coordinateDimension();
    const std::size_t xDim = map ? map->stateDimension() : qDim;
    if (map && map->outputDimension() != qDim) throw std::invalid_argument("coordinate map width differs from the representation domain");
    descriptor.complexValued = true;
    if (descriptor.outputNames.empty())
        for (std::size_t r = 0; r < outputs; ++r)
            descriptor.outputNames.push_back(outputs == 1 ? descriptor.name : descriptor.name + "[" + std::to_string(r) + "]");
    Evaluator evaluator = [rep, map](const double* x) {
        if (!map) return rep->evaluate(x);
        const auto q = map->evaluateReal(x);
        return rep->evaluate(q.data());
    };
    Jacobian jacobian;
    if (rep->exactJacobian() && (!map || map->hasAnalyticJacobian()))
        jacobian = [rep, map, outputs, qDim, xDim](const double* x) {
            if (!map) return rep->jacobian(x);
            const auto q = map->evaluateReal(x);
            const auto inner = rep->jacobian(q.data());                               // outputs x qDim
            const auto outer = map->jacobian(x, ObservableJacobianBackend::analytic);  // qDim x xDim
            Values out(outputs * xDim);
            for (std::size_t r = 0; r < outputs; ++r)
                for (std::size_t a = 0; a < qDim; ++a)
                    for (std::size_t i = 0; i < xDim; ++i) out[r * xDim + i] += inner[r * qDim + a] * outer[a * xDim + i];
            return out;
        };
    Observable out(std::move(descriptor), xDim, outputs, std::move(evaluator), std::move(jacobian));
    out.representationName_ = ::dynsyslab::core::representationName(*rep);
    out.representation_ = std::move(rep);
    return out;
}

std::string representationName(const FunctionRepresentation& representation) {
    if (representation.kind() == "interpolation") return "interpolation";
    if (representation.coordinates().window) return "time-delay";
    if (const auto* expansion = dynamic_cast<const FunctionExpansion*>(&representation)) {
        const auto kinds = expansion->dictionary().kinds();
        if (kinds.size() == 1) return kinds.front() == "monomial" ? "polynomial" : kinds.front();
    }
    return representation.kind();
}

std::string Observable::representationName() const {
    return representation_ ? representationName_ : std::string("evaluator");
}

std::vector<std::string> Observable::representationNames() const {
    std::vector<std::string> names{representationName()};
    for (const auto& alternative : alternatives_) names.push_back(alternative.name);
    return names;
}

Observable::Alternative Observable::current() const {
    return {representationName(), stateDimension_, evaluator_, analyticJacobian_, stateComponents_, representation_};
}

void Observable::use(Alternative alternative) {
    representationName_ = alternative.representation ? alternative.name : std::string();
    stateDimension_ = alternative.stateDimension;
    evaluator_ = std::move(alternative.evaluator);
    analyticJacobian_ = std::move(alternative.jacobian);
    stateComponents_ = std::move(alternative.stateComponents);
    representation_ = std::move(alternative.representation);
}

Observable Observable::withRepresentation(std::shared_ptr<const FunctionRepresentation> representation,
                                          const ObservableCatalog* catalog, bool asDefault, std::string name) const {
    auto descriptor = descriptor_;
    descriptor.outputNames.clear();
    const auto built = fromRepresentation(std::move(descriptor), std::move(representation), catalog);
    if (built.outputDimension_ != outputDimension_)
        throw std::invalid_argument(descriptor_.name + ": another representation must have the same outputs");
    auto alternative = built.current();
    if (!name.empty()) alternative.name = std::move(name);
    const auto names = representationNames();
    if (std::find(names.begin(), names.end(), alternative.name) != names.end())
        throw std::invalid_argument(descriptor_.name + " already has a " + alternative.name + " representation");
    auto out = *this;
    if (asDefault) {
        out.alternatives_.insert(out.alternatives_.begin(), out.current());
        out.use(std::move(alternative));
    } else {
        out.alternatives_.push_back(std::move(alternative));
    }
    return out;
}

Observable Observable::as(const std::string& name) const {
    if (name == representationName()) return *this;
    for (std::size_t k = 0; k < alternatives_.size(); ++k)
        if (alternatives_[k].name == name) {
            auto out = *this;
            auto chosen = std::move(out.alternatives_[k]);
            out.alternatives_[k] = out.current();
            out.use(std::move(chosen));
            return out;
        }
    std::string known;
    for (const auto& n : representationNames()) known += (known.empty() ? "" : ", ") + n;
    throw std::invalid_argument(descriptor_.name + " has no " + name + " representation; it has: " + known);
}

std::vector<AnalyticExpression> Observable::expressions() const {
    const auto* expansion = dynamic_cast<const FunctionExpansion*>(representation_.get());
    if (expansion && !descriptor_.complexValued && expansion->coordinates().isPhysicalState()) {
        const auto parts = expansion->realExpressions();  // real, imaginary per output
        if (!parts) return {};
        std::vector<AnalyticExpression> out;
        for (std::size_t r = 0; r < parts->size(); r += 2) {
            const auto imaginary = (*parts)[r + 1].toJson();
            if (imaginary.at("operation").string() != "constant" || imaginary.at("scalar").number() != 0) return {};
            out.push_back((*parts)[r]);
        }
        return out;
    }
    if (representation_ || stateComponents_.size() != outputDimension_) return {};
    std::vector<AnalyticExpression> out;
    for (const auto& index : stateComponents_) {
        if (!index) return {};
        out.push_back(AnalyticExpression::state(stateDimension_, *index));
    }
    return out;
}

std::string Observable::representationKind() const {
    return representation_ ? representation_->kind() : "evaluator";
}

Eigenfunction::Eigenfunction(Observable function, Eigenvalue eigenvalue)
    : Observable(std::move(function)), eigenvalue_(eigenvalue) {
    if (outputDimension() != 1) throw std::invalid_argument("an eigenfunction is one scalar function");
}

Eigenfunction Eigenfunction::fromRepresentation(ObservableDescriptor descriptor, std::shared_ptr<const FunctionRepresentation> rep,
                                                std::size_t output, Eigenvalue eigenvalue,
                                                const ObservableCatalog* catalog) {
    if (!rep || output >= rep->outputCount()) throw std::out_of_range("eigenfunction output");
    auto single = rep->outputCount() == 1 ? rep : rep->selectOutputs({output});
    return Eigenfunction(Observable::fromRepresentation(std::move(descriptor), std::move(single), catalog), eigenvalue);
}

}  // namespace dynsyslab::core
