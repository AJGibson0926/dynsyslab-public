#pragma once
#include <dynsyslab/core/Hankel.hpp>
#include <cstddef>
#include <memory>
#include <vector>

namespace dynsyslab::core {

// The Hankel module has two realizations of one matrix. buildHankelMatrix
// stores X (and the shifted Y) explicitly; HankelOperator applies the same X
// without storing it, for embeddings too large to materialize. Both use the
// layout documented in Hankel.hpp: row d*C+c holds channel c at delay block d,
// oldest block first, divided by the channel scale.
//
// A HankelOperator pools independent segments (trajectories). Columns are
// segment-local windows concatenated in segment order; no window ever crosses
// a segment boundary. Segments are already on the analysis grid (any sample
// stride is applied by the caller) and use consecutive delays and snapshots.
struct HankelSegment {
    std::vector<double> values; // samples x inputs, row-major
    std::size_t samples{};
    std::size_t inputs{};
};

class HankelOperator {
public:
    HankelOperator(const std::vector<HankelSegment>& segments,
                   std::vector<HankelChannel> channels, std::size_t delays);
    ~HankelOperator();
    HankelOperator(HankelOperator&&) noexcept;
    HankelOperator& operator=(HankelOperator&&) noexcept;

    [[nodiscard]] std::size_t rows() const noexcept { return delays_ * channels_.size(); }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] std::size_t delays() const noexcept { return delays_; }
    [[nodiscard]] const std::vector<HankelChannel>& channels() const noexcept { return channels_; }
    [[nodiscard]] const std::vector<std::size_t>& segmentColumns() const noexcept { return segmentColumns_; }
    // Squared Frobenius norm of X.
    [[nodiscard]] double energy() const noexcept { return energy_; }

    // X * v for v (columns x k, row-major); result is rows x k, row-major.
    [[nodiscard]] std::vector<double> apply(const std::vector<double>& v, std::size_t k) const;
    // X^T * w for w (rows x k, row-major); result is columns x k, row-major.
    [[nodiscard]] std::vector<double> applyAdjoint(const std::vector<double>& w, std::size_t k) const;

private:
    struct Plans;
    std::size_t delays_{};
    std::vector<HankelChannel> channels_;
    std::size_t columns_{};
    std::vector<std::size_t> segmentColumns_;
    double energy_{};
    std::unique_ptr<Plans> plans_;
};

struct HankelBasisOptions {
    std::size_t rank{4};
    std::size_t oversampling{4};
    std::size_t subspaceIterations{8};
    std::size_t seed{20261002};
};

// Leading left singular vectors of X by randomized subspace iteration on X X^T.
// Signs are fixed so each vector's largest-magnitude entry is positive. The one Hankel basis type:
// a pooled history basis (HistoryObservation.hpp) is the same with rows = delays, columns = the
// training windows, the number of training trajectories and the observation's fitting scale.
struct HankelLeftBasis {
    std::vector<double> leftSingularVectors; // rows x rank, row-major
    std::vector<double> singularValues;
    std::vector<double> cumulativeEnergy;     // fraction of energy() through each mode
    double totalEnergy{};
    std::size_t rows{};
    std::size_t rank{};
    std::size_t columns{};
    std::size_t trajectories{1};              // independent runs the windows came from
    double observationScale{1.0};             // fitting normalization of the observation (not reapplied)
};

[[nodiscard]] HankelLeftBasis fitHankelLeftBasis(const HankelOperator& hankel,
                                                 const HankelBasisOptions& options = {});

} // namespace dynsyslab::core
