#pragma once
// State readouts from measured features, and errors scored oscillation by oscillation
// (docs/architecture/RECEIVER_ARRAYS.md). Ported from the receiver-reconstruction side-car.
//
//   NeighborReadout  k-nearest-neighbour regression (Cover and Hart 1967, IEEE Trans. Inf. Theory
//                    13:21-27) in asinh(feature / scale) coordinates (the inverse hyperbolic sine
//                    transform: linear near zero, logarithmic in the tails; Burbidge, Magee and
//                    Robb 1988, J. Am. Stat. Assoc. 83:123-127), with inverse-square distance
//                    weights (Shepard 1968, Proc. 23rd ACM National Conference). An exact k-d tree
//                    search (Bentley 1975, Commun. ACM 18:509-517) with ties broken by row index;
//                    every estimate keeps its neighbour rows and weights.
//   oscillationMaxima  times where a velocity component changes from positive to nonpositive
//                    (maxima of its position), linearly interpolated; a resting start counts.
//   intervalErrors   per half-open interval [b_i, b_{i+1}): sample count, partial flag, error and
//                    reference sums of squares and the relative RMS percentage per component.
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

class NeighborReadout {
public:
    // features (rows x dimension), states (rows x stateDimension), positive scales (dimension).
    NeighborReadout(std::vector<double> features, std::vector<double> states, std::vector<double> scales, std::size_t stateDimension);
    struct Estimate {
        std::vector<double> values;          // stateDimension
        std::vector<std::size_t> rows;       // the k training rows, nearest first
        std::vector<double> weights;         // normalized, summing to 1
        double nearestDistance{};            // in the scaled coordinates
    };
    [[nodiscard]] Estimate predict(const double* features, std::size_t k) const;
    [[nodiscard]] std::size_t dimension() const noexcept { return scales_.size(); }
    [[nodiscard]] std::size_t stateDimension() const noexcept { return stateDimension_; }
    [[nodiscard]] std::size_t rows() const noexcept { return order_.size(); }

private:
    struct Hit { double squared; std::size_t row; };
    void build(std::size_t first, std::size_t last, std::size_t axis);
    void search(std::size_t first, std::size_t last, std::size_t axis, const std::vector<double>& query, std::size_t k,
                std::vector<Hit>& hits) const;
    std::vector<double> points_, states_, scales_;
    std::size_t stateDimension_;
    std::vector<std::size_t> order_;
};

[[nodiscard]] std::vector<double> oscillationMaxima(const std::vector<double>& time, const std::vector<double>& states,
                                                    std::size_t width, std::size_t velocity);

struct IntervalErrors {
    std::vector<double> count, partial;            // intervals
    std::vector<double> errorSquares, referenceSquares, relativeRmsPercent;  // intervals x width
    std::size_t width{};
};

// Scores `prediction` against `reference` (both samples x width on `time`) over consecutive
// intervals of `bounds`; the reference deviation is measured from `center` (empty: zero).
[[nodiscard]] IntervalErrors intervalErrors(const std::vector<double>& time, const std::vector<double>& reference,
                                            const std::vector<double>& prediction, std::size_t width,
                                            const std::vector<double>& bounds, std::vector<double> center = {});

}  // namespace dynsyslab::core
