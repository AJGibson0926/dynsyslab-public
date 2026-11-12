#pragma once
#include <cstddef>
#include <vector>

namespace receivers {
struct Neighbor {
    double squaredDistance;
    std::size_t row;
};
struct Estimate {
    std::vector<double> values;
    std::vector<Neighbor> neighbors;
    std::vector<double> weights;
};

// A training-only index for arbitrary observation and state dimensions.
// Row-major matrices, asinh(features/scales), and normalized inverse-square weights.
class Lookup {
public:
    Lookup(std::vector<double> features, std::vector<double> states,
           std::vector<double> scales, std::size_t stateDimension);
    Estimate predict(const double* features, std::size_t k) const;
private:
    std::vector<double> points_, states_, scales_;
    std::size_t stateDimension_;
    std::vector<std::size_t> order_;
    void build(std::size_t first, std::size_t last, std::size_t axis);
    void search(std::size_t first, std::size_t last, std::size_t axis,
                const std::vector<double>& query, std::size_t k, std::vector<Neighbor>& hits) const;
};
}
