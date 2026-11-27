#pragma once
// Derivative-consistent scalar readout (docs/architecture/IDENTIFICATION.md, "Delay diagnostics").
// Each training row i carries a label g_i and a local slope, fitted so its first-order
// approximation matches the label rate along the training velocity; the slope of row i is
// a_i = dg_i dx_i / |dx_i|^2 in the asinh(feature / scale) coordinates x (NeighborReadout), where
// dx_i = v_i / sqrt(scale^2 + b_i^2) is the feature velocity pushed through the asinh map.
// A query blends the k nearest rows with Wendland weights w = (1 - r)^4 (1 + 4 r), r = d / d_{k+1}
// (Wendland, Adv. Comput. Math. 4:389-396, 1995), of the bounded local approximants
// g_i + A_i tanh(a_i . (x - x_i) / A_i), A_i = max(floor, |g_i|), and returns the value with its
// exact rate along the query velocity (the derivative of the blend, weights included).
#include <dynsyslab/core/NeighborReadout.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

class TangentReadout {
public:
    // features, rates: rows x dimension; labels, labelRates: rows; scales: dimension (empty: RMS).
    TangentReadout(std::vector<double> features, const std::vector<double>& rates, std::vector<double> labels,
                   const std::vector<double>& labelRates, std::vector<double> scales = {}, double amplitudeFloor = 1.0);
    struct Phase { double value{}, rate{}; };
    [[nodiscard]] Phase phase(const double* features, const double* rates, std::size_t k = 16) const;
    [[nodiscard]] std::size_t dimension() const noexcept { return scales_.size(); }
    [[nodiscard]] const std::vector<double>& scales() const noexcept { return scales_; }

private:
    std::vector<double> x_, slopes_, labels_, amplitude_, scales_;
    NeighborReadout search_;
};

}  // namespace dynsyslab::core
