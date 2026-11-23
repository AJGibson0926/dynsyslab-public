#pragma once
#include <cstddef>
#include <vector>
namespace dynsyslab::numerics {
// Sequential row-block Householder QR. R^T R = A^T W A up to roundoff;
// no normal equations, numerical rank truncation, or weight normalization.
// Append already evaluated features (and targets as extra columns if needed).
class WeightedQrAccumulator {
public:
    explicit WeightedQrAccumulator(std::size_t columns);
    void append(const std::vector<double>& values,std::size_t rows,
                const std::vector<double>& weights={});
    const std::vector<double>& factor() const {return factor_;} // row-major R
    std::size_t columns() const {return columns_;}
    std::size_t rows() const {return rows_;}
    std::size_t samples() const {return samples_;}
    double totalWeight() const {return totalWeight_;}
private:
    std::size_t columns_{},rows_{},samples_{};
    double totalWeight_{};
    std::vector<double> factor_;
};
}
