#pragma once

#include <string>
#include <vector>

namespace dynsyslab::core {

struct HarmonicStateTarget {
    std::string name;
    double offset{};
    std::vector<double> amplitudes;
    std::vector<double> angularFrequencies;
    std::vector<double> phases;

    [[nodiscard]] std::vector<double> evaluate(double time) const;
};

[[nodiscard]] const std::vector<HarmonicStateTarget>& rpeControlTargets();
[[nodiscard]] const HarmonicStateTarget& rpeControlTarget(const std::string& name);

}  // namespace dynsyslab::core
