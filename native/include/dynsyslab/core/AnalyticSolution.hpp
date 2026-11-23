#pragma once
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {
using AnalyticState = std::vector<double>;
struct AnalyticSolutionDescriptor {
    std::string identifier, label, formula, conditions, evaluation;
    std::vector<std::string> references;
    std::string version{"1"};
};
struct AnalyticTrajectory {
    // Absolute native model time. Empty evaluator means inapplicable, with reason.
    std::function<AnalyticState(double)> evaluate;
    std::string reason;
    // Physical/native model-time period when known, zero when not supplied.
    // Never the auxiliary elliptic parameter period; equilibria leave this zero.
    double period{};
    // Optional invariant-torus map; angles in radians with origin at the IC.
    // Frequencies advance in native/physical model time. No common period is
    // asserted for a multifrequency torus, even if its sampled trace nearly closes.
    std::vector<double> angularFrequencies;
    std::function<AnalyticState(const std::vector<double>&)> torus;
    explicit operator bool() const { return static_cast<bool>(evaluate); }
};
struct AnalyticSolutionFamily {
    AnalyticSolutionDescriptor descriptor;
    // Bound to a system's resolved equations/parameters. IC and native initial time
    // remain explicit, allowing every trajectory in a batch to use the same property.
    std::function<AnalyticTrajectory(const AnalyticState &, double)> instantiate;
};
struct AnalyticComparison {
    std::string status{"not_applicable"}, reason;
    AnalyticSolutionDescriptor descriptor;
    std::vector<double> values; // sample-major, same coordinates and mesh as numerical data
    std::vector<double> valid;  // one entry per sample; invalid values are NaN, never extrapolated
    std::vector<double> maximumAbsoluteError;
    double period{};
};
AnalyticComparison compareAnalyticSolution(const std::vector<AnalyticSolutionFamily> &families,
                                           const AnalyticState &initial, double initialTime,
                                           const std::vector<double> &time,
                                           const std::vector<double> &numerical,
                                           const std::string &selection = "auto");
} // namespace dynsyslab::core
