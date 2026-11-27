#pragma once
// An Eigenfunction is an Observable with an Eigenvalue. Everything an Observable does (evaluation,
// Jacobians, partial and flow derivatives, its named representations: polynomial, interpolation,
// time-delay, ...) it does; in addition its value evolves as exp(lambda t) along solutions. It
// carries every representation it has; the first is the default.
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
    // The same eigenfunction with another representation added, or another as its default.
    [[nodiscard]] Eigenfunction withRepresentation(std::shared_ptr<const FunctionRepresentation> representation,
        const ObservableCatalog* catalog = nullptr, bool asDefault = false, std::string name = {}) const {
        return Eigenfunction(Observable::withRepresentation(std::move(representation), catalog, asDefault, std::move(name)), eigenvalue_);
    }
    [[nodiscard]] Eigenfunction as(const std::string& name) const { return Eigenfunction(Observable::as(name), eigenvalue_); }
    // psi(x(t)) predicted from psi(x(0)).
    [[nodiscard]] std::complex<double> propagate(std::complex<double> value, double elapsed) const {
        return eigenvalue_.transition(elapsed) * value;
    }

private:
    Eigenvalue eigenvalue_;
};

}  // namespace dynsyslab::core
