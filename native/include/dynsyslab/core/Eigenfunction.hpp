#pragma once
// An Eigenfunction is an Observable with an Eigenvalue. Everything an Observable does (evaluation,
// Jacobians, partial and flow derivatives, its representation property: expansion or
// interpolation) it does; in addition its value evolves as exp(lambda t) along solutions.
#include <dynsyslab/core/Eigenvalue.hpp>
#include <dynsyslab/core/Observable.hpp>

namespace dynsyslab::core {

class Eigenfunction : public Observable {
public:
    // `function` must have one output (complex valued).
    Eigenfunction(Observable function, Eigenvalue eigenvalue);
    // Output `output` of a representation, on the physical state through the domain's state map
    // (Observable::fromRepresentation).
    [[nodiscard]] static Eigenfunction fromRepresentation(ObservableDescriptor descriptor,
        std::shared_ptr<const FunctionRepresentation> representation, std::size_t output, Eigenvalue eigenvalue,
        const ObservableCatalog* catalog = nullptr);

    [[nodiscard]] const Eigenvalue& eigenvalue() const noexcept { return eigenvalue_; }
    // psi(x(t)) predicted from psi(x(0)).
    [[nodiscard]] std::complex<double> propagate(std::complex<double> value, double elapsed) const {
        return eigenvalue_.transition(elapsed) * value;
    }

private:
    Eigenvalue eigenvalue_;
};

}  // namespace dynsyslab::core
