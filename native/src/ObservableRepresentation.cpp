// Observables backed by a FunctionRepresentation (expansion or interpolation), and Eigenfunction.
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/FunctionRepresentation.hpp>

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
    out.representation_ = std::move(rep);
    return out;
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
