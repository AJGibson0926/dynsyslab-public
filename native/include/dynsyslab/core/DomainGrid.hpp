#pragma once
// Functions sampled on a grid of their domain: the one service Viewer heat maps, domain plots,
// push-forwards and Koopman-action panels evaluate through.
//
//   DomainGrid  a rectangular grid on some axes of a FunctionDomain. Coordinates that are not
//               gridded stay at `anchor`, so a grid is also a slice of a higher-dimensional domain.
//               The first axis varies fastest (row-major images with axis 1 as rows).
//   GridField   complex values at every grid point (points x outputs) and a status per point:
//               valid, undefined (evaluation failed or was not finite), conflict (a push-forward
//               whose fibre carries different values) or noPreimage (no source point maps there).
#include <dynsyslab/core/FunctionDomain.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {

class FunctionRepresentation;

struct GridAxis {
    std::size_t coordinate{};  // index into the domain coordinates
    double lower{}, upper{};
    std::size_t points{101};
    [[nodiscard]] double value(std::size_t index) const noexcept {
        return points < 2 ? lower : lower + (upper - lower) * double(index) / double(points - 1);
    }
};

struct DomainGrid {
    FunctionDomain domain;
    std::vector<double> anchor;  // one value per domain coordinate
    std::vector<GridAxis> axes;

    [[nodiscard]] std::size_t dimension() const noexcept { return anchor.size(); }
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::vector<double> axisValues(std::size_t axis) const;
    void point(std::size_t index, double* out) const;
    [[nodiscard]] std::vector<double> points() const;  // size() x dimension(), row-major
    // Index of the grid cell containing `q` (by nearest node on each axis), or size() if outside.
    [[nodiscard]] std::size_t locate(const double* q) const;
    void validate() const;

    // Every coordinate gridded over [lower, upper] with `points` nodes per axis.
    [[nodiscard]] static DomainGrid box(FunctionDomain domain, std::vector<double> lower, std::vector<double> upper,
                                        std::size_t points = 101);
    // Bounds of sample rows (samples x dimension) widened by `padding` of each extent.
    [[nodiscard]] static DomainGrid around(FunctionDomain domain, const std::vector<double>& samples,
                                           std::size_t points = 101, double padding = 0.04);
};

enum class GridStatus : std::uint8_t { valid = 0, undefined = 1, conflict = 2, noPreimage = 3 };

struct GridField {
    DomainGrid grid;
    std::size_t outputs{};
    std::vector<std::string> outputNames;
    bool complexValued{};
    std::vector<std::complex<double>> values;  // points x outputs
    std::vector<GridStatus> status;            // per point
    std::vector<double> spread;                // per point and output, push-forward only
    std::vector<std::size_t> counts;           // preimage samples per point, push-forward only

    [[nodiscard]] std::complex<double> value(std::size_t point, std::size_t output = 0) const {
        return values[point * outputs + output];
    }
    [[nodiscard]] std::size_t count(GridStatus which) const noexcept;
    // One output as a scalar field (status kept), e.g. to plot one mode.
    [[nodiscard]] GridField output(std::size_t index) const;
};

using PointFunction = std::function<std::vector<std::complex<double>>(const double*)>;

// `f` at every grid point; points where it throws or is not finite are `undefined` (NaN values).
[[nodiscard]] GridField evaluateOnGrid(const PointFunction& f, std::size_t outputs, DomainGrid grid);
// An observable of the grid's coordinates (the grid lives in its state space).
[[nodiscard]] GridField evaluateOnGrid(const Observable& observable, DomainGrid grid);
// A representation on its own domain q (the grid lives in q).
[[nodiscard]] GridField evaluateOnGrid(const FunctionRepresentation& representation, DomainGrid grid);

}  // namespace dynsyslab::core
