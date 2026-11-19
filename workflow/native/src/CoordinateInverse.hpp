#pragma once
// The inverse of a coordinate change, retained with its archive and callable.
//
// Every coordinate change records metadata "coordinate_inverse" stating what
// inverse exists:
//   exact     identity: the coordinates are the previous coordinates.
//   readout   history coordinates: a fitted map back to the previous
//             coordinates (inverse-distance weighting over the nearest training
//             pairs). Approximate, and valid only within the training range.
//   none      no inverse is retained (with the reason), e.g. an analytic map the
//             diagnostics show is not one-to-one, or a column subset.
// loadCoordinateInverse reads it back; apply maps rows of these coordinates to
// rows of the previous coordinates. dynsyslab-cli coordinate-inverse and the
// sindy_dynamics setting rollout_inverse call it.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

class CoordinateInverse {
public:
    CoordinateInverse();
    ~CoordinateInverse();
    CoordinateInverse(CoordinateInverse&&) noexcept;
    CoordinateInverse& operator=(CoordinateInverse&&) noexcept;
    [[nodiscard]] const std::string& kind() const noexcept { return kind_; }
    [[nodiscard]] std::size_t inputs() const noexcept { return inputs_; }
    [[nodiscard]] std::size_t outputs() const noexcept { return outputs_; }
    [[nodiscard]] const std::vector<std::string>& targetNames() const noexcept { return names_; }
    [[nodiscard]] const std::vector<std::string>& targetLabels() const noexcept { return labels_; }
    [[nodiscard]] const json::Value& description() const noexcept { return description_; }
    // rows x inputs (row-major) -> rows x outputs.
    [[nodiscard]] std::vector<double> apply(const std::vector<double>& values, std::size_t rows) const;
    [[nodiscard]] dynsyslab::core::TimeSeries apply(const dynsyslab::core::TimeSeries& series) const;

    friend CoordinateInverse loadCoordinateInverse(const std::filesystem::path& archive);
private:
    struct Readout;
    std::string kind_;
    std::size_t inputs_{}, outputs_{};
    std::vector<std::string> names_, labels_;
    json::Value description_;
    std::unique_ptr<Readout> readout_;
};

// Throws when the archive records no inverse.
CoordinateInverse loadCoordinateInverse(const std::filesystem::path& archive);

// Record the inverse of a coordinate change being written: an exact identity,
// or (history kinds with an earlier coordinate archive as source) a readout
// trained on pairs of the new coordinates and the source's columns at the same
// times. `kind` "none" records only the reason.
void recordCoordinateInverse(const SourceData& source, const std::string& kind, const std::string& reason,
    const dynsyslab::core::TimeSeries& coordinates, json::Value& meta, dynsyslab::core::NpzRecord& archive);

// dynsyslab-cli coordinate-inverse COORDINATES.npz INPUT.npz OUTPUT.npz
void writeCoordinateInverse(const std::filesystem::path& coordinates, const std::filesystem::path& input,
                            const std::filesystem::path& output);

} // namespace dynsyslab::workflow::postprocess_detail
