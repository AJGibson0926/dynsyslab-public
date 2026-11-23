#pragma once
#include <functional>
#include <vector>
namespace dynsyslab::core {
// Deterministic sequential-importance filter. Caller owns dynamics, observation
// and noise models. No truth-state access, implicit resampling or hidden jitter.
// Doucet, de Freitas & Gordon (2001), Sequential Monte Carlo Methods in Practice.
class EnsembleBayesFilter {
  public:
    using State=std::vector<double>;
    using Transition=std::function<State(const State&)>;
    using Readout=std::function<State(const State&)>;
    explicit EnsembleBayesFilter(std::vector<State> particles,std::vector<double> prior={});
    void predict(const Transition&);
    // Independent Gaussian channels; variance must be positive. Returns log evidence
    // conditional on the existing normalized posterior (including normalizing constants).
    double correct(const State& measurement,const State& variance,const Readout&);
    State mean(const Readout&) const;
    double effectiveSampleSize() const;
    std::vector<double> weights() const;
    const std::vector<State>& particles() const { return particles_; }
  private:
    std::vector<State> particles_;
    std::vector<double> logWeights_;
};
}
