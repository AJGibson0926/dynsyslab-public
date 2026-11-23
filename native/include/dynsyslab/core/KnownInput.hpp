#pragma once
#include <dynsyslab/core/Integrator.hpp>
#include <type_traits>
#include <utility>

namespace dynsyslab::core {
// Applied-input history and the input associated with a sensor sample are
// distinct at a command switch. Analytic inputs use the same value for both.
// ZOH convention: Franklin, Powell & Workman, Digital Control of Dynamic
// Systems, 3rd ed. (1998). Each held interval uses the existing native solver.
class KnownScalarInput {
public:
    KnownScalarInput() = default;
    template<class F, std::enable_if_t<!std::is_same_v<std::decay_t<F>, KnownScalarInput>
        && std::is_invocable_r_v<double,F,double>,int> = 0>
    KnownScalarInput(F&& f) : applied_(std::forward<F>(f)), measured_(applied_) {}
    static KnownScalarInput recorded(std::vector<double> time, std::vector<double> values,
        double beforeFirst, bool measurementBeforeCommand = true);
    double operator()(double time) const;
    double measurement(double time) const;
    explicit operator bool() const noexcept { return bool(applied_); }
    const std::vector<double>& breakpoints() const noexcept { return breaks_; }
private:
    std::function<double(double)> applied_, measured_;
    std::vector<double> breaks_;
};

// One trajectory. Uses the ordinary integrator for continuous inputs; for a
// recorded ZOH history it restarts exactly at switches, with the old interval's
// input on its right endpoint. Output samples preserve the requested grid.
// No extrapolation outside the recorded time range and no hidden final hold.
NativeSolution integrateKnownInput(ODESystem&, const std::vector<double>& initial,
    const IntegrationOptions&, const KnownScalarInput&);
}
