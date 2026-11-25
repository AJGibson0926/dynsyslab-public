#pragma once
// Hankel basis: the object H-SVD produces (Broomhead & King, Physica D 20, 1986; Vautard & Ghil,
// Physica D 35, 1989). A window h is projected onto the leading left singular vectors U of the
// Hankel matrix:
//   energy scaling    q  = U^T (h - mean)          (H-SVD scores)
//   whitened scaling  xi = Sigma^-1 U^T (h - mean)  (history-basis latent coordinates)
// H-DMD's eigenfunctions are psi = W^-1 q for the same, uncentered window matrix.
#include <dynsyslab/core/HankelRepresentations.hpp>

namespace dynsyslab::core {

enum class HankelCoordinateScaling { energy, whitened };

[[nodiscard]] HankelCoordinateScaling hankelCoordinateScaling(const std::string& name);
[[nodiscard]] std::string hankelCoordinateScalingName(HankelCoordinateScaling scaling);

inline constexpr const char* kHankelBasisSchema = "dynsyslab.hankel-basis/v1";

class HankelBasis {
public:
    // modes: rows x rank (U, window-space directions); mean: rows, or empty when uncentered.
    HankelBasis(HankelWindowLayout layout, std::vector<double> modes, std::vector<double> singularValues,
                std::vector<double> mean, DataOrigin origin);
    [[nodiscard]] const HankelWindowLayout& layout() const noexcept { return layout_; }
    [[nodiscard]] std::size_t rank() const noexcept { return singular_.size(); }
    [[nodiscard]] const std::vector<double>& modes() const noexcept { return modes_; }
    [[nodiscard]] const std::vector<double>& singularValues() const noexcept { return singular_; }
    [[nodiscard]] const std::vector<double>& mean() const noexcept { return mean_; }
    [[nodiscard]] const DataOrigin& origin() const noexcept { return origin_; }
    // Coordinates of one raw window (as from HankelWindowLayout::window_at; channel scales applied here).
    [[nodiscard]] std::vector<double> coordinates(const double* window, HankelCoordinateScaling scaling) const;
    // Coordinates of every complete window of a record labelled at `labels` (record rows).
    [[nodiscard]] std::vector<double> coordinates(const TimeSeries& record, const std::vector<std::size_t>& labels,
                                                  HankelCoordinateScaling scaling) const;
    // Metadata block; arrays are named by the caller's archive (see HankelSvdPostprocess).
    [[nodiscard]] json::Value describe() const;

private:
    HankelWindowLayout layout_;
    std::vector<double> modes_, singular_, mean_;
    DataOrigin origin_;
};

// Rebuilds a basis from an H-SVD archive's metadata block and arrays.
[[nodiscard]] HankelBasis hankelBasisFromArchive(const json::Value& block, std::vector<double> modes,
                                                 std::vector<double> singularValues, std::vector<double> mean);

}  // namespace dynsyslab::core
