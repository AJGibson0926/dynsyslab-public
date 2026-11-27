#pragma once

#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::core {

struct HarmonicStateTarget {
    std::string name;
    double offset{};
    std::vector<double> amplitudes;
    std::vector<double> angularFrequencies;
    std::vector<double> phases;

    [[nodiscard]] std::vector<double> evaluate(double time) const;
    // Analytic bounds containing the whole target: radius offset -+ sum|a_k|, velocity
    // -+ sum|a_k w_k| (the triangle inequality on the sum of sinusoids).
    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> conservativeEnvelope() const;
};

[[nodiscard]] const std::vector<HarmonicStateTarget>& rpeControlTargets();
[[nodiscard]] const HarmonicStateTarget& rpeControlTarget(const std::string& name);

}  // namespace dynsyslab::core
