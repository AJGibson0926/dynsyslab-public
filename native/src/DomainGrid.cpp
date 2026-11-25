#include <dynsyslab/core/DomainGrid.hpp>
#include <dynsyslab/core/FunctionRepresentation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {

std::size_t DomainGrid::size() const noexcept {
    std::size_t n = 1;
    for (const auto& axis : axes) n *= axis.points;
    return axes.empty() ? 0 : n;
}

std::vector<double> DomainGrid::axisValues(std::size_t axis) const {
    const auto& a = axes.at(axis);
    std::vector<double> out(a.points);
    for (std::size_t i = 0; i < a.points; ++i) out[i] = a.value(i);
    return out;
}

void DomainGrid::point(std::size_t index, double* out) const {
    std::copy(anchor.begin(), anchor.end(), out);
    for (const auto& axis : axes) {
        out[axis.coordinate] = axis.value(index % axis.points);
        index /= axis.points;
    }
}

std::vector<double> DomainGrid::points() const {
    std::vector<double> out(size() * dimension());
    for (std::size_t k = 0; k < size(); ++k) point(k, out.data() + k * dimension());
    return out;
}

std::size_t DomainGrid::locate(const double* q) const {
    std::size_t index = 0, stride = 1;
    for (const auto& axis : axes) {
        const double span = axis.upper - axis.lower, x = q[axis.coordinate];
        if (!std::isfinite(x)) return size();
        const double t = axis.points < 2 ? 0.0 : (x - axis.lower) / span * double(axis.points - 1);
        const double i = std::round(t);
        if (i < 0 || i > double(axis.points - 1)) return size();
        index += std::size_t(i) * stride;
        stride *= axis.points;
    }
    // Off-axis coordinates must lie on the slice (within half a cell of the anchor is not checked).
    return index;
}

void DomainGrid::validate() const {
    if (anchor.size() != domain.dimension || axes.empty())
        throw std::invalid_argument("A domain grid needs an anchor per coordinate and at least one axis");
    for (std::size_t j = 0; j < axes.size(); ++j) {
        const auto& a = axes[j];
        if (a.coordinate >= dimension() || !a.points || !std::isfinite(a.lower) || !std::isfinite(a.upper)
            || (a.points > 1 && !(a.upper > a.lower)))
            throw std::invalid_argument("Grid axes need finite increasing bounds on a domain coordinate");
        for (std::size_t k = 0; k < j; ++k)
            if (axes[k].coordinate == a.coordinate) throw std::invalid_argument("A coordinate is gridded twice");
    }
}

DomainGrid DomainGrid::box(FunctionDomain domain, std::vector<double> lower, std::vector<double> upper, std::size_t points) {
    if (lower.size() != domain.dimension || upper.size() != domain.dimension)
        throw std::invalid_argument("Grid bounds need one value per domain coordinate");
    DomainGrid grid;
    grid.anchor = lower;
    for (std::size_t j = 0; j < lower.size(); ++j) grid.axes.push_back({j, lower[j], upper[j], points});
    grid.domain = std::move(domain);
    grid.validate();
    return grid;
}

DomainGrid DomainGrid::around(FunctionDomain domain, const std::vector<double>& samples, std::size_t points, double padding) {
    const auto d = domain.dimension;
    if (!d || samples.empty() || samples.size() % d) throw std::invalid_argument("Grid bounds need sample rows of the domain");
    std::vector<double> lower(d, std::numeric_limits<double>::infinity()), upper(d, -std::numeric_limits<double>::infinity());
    for (std::size_t k = 0; k < samples.size(); ++k) {
        if (!std::isfinite(samples[k])) continue;
        lower[k % d] = std::min(lower[k % d], samples[k]);
        upper[k % d] = std::max(upper[k % d], samples[k]);
    }
    for (std::size_t j = 0; j < d; ++j) {
        if (!std::isfinite(lower[j])) throw std::invalid_argument("Grid bounds need finite samples");
        const double pad = padding * std::max(upper[j] - lower[j], 1e-12);
        lower[j] -= pad;
        upper[j] += pad;
    }
    return box(std::move(domain), std::move(lower), std::move(upper), points);
}

std::size_t GridField::count(GridStatus which) const noexcept {
    return std::size_t(std::count(status.begin(), status.end(), which));
}

GridField GridField::output(std::size_t index) const {
    if (index >= outputs) throw std::out_of_range("grid field output");
    GridField out;
    out.grid = grid;
    out.outputs = 1;
    if (index < outputNames.size()) out.outputNames = {outputNames[index]};
    out.complexValued = complexValued;
    out.status = status;
    out.counts = counts;
    for (std::size_t k = 0; k < status.size(); ++k) {
        out.values.push_back(values[k * outputs + index]);
        if (!spread.empty()) out.spread.push_back(spread[k * outputs + index]);
    }
    return out;
}

GridField evaluateOnGrid(const PointFunction& f, std::size_t outputs, DomainGrid grid) {
    grid.validate();
    if (!outputs) throw std::invalid_argument("A grid field needs at least one output");
    GridField out;
    const auto n = grid.size(), d = grid.dimension();
    out.outputs = outputs;
    out.values.assign(n * outputs, {});
    out.status.assign(n, GridStatus::valid);
    std::vector<double> q(d);
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < n; ++k) {
        grid.point(k, q.data());
        std::vector<std::complex<double>> v;
        bool ok = true;
        try { v = f(q.data()); } catch (const std::exception&) { ok = false; }
        ok = ok && v.size() == outputs
             && std::all_of(v.begin(), v.end(), [](auto z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); });
        if (!ok) {
            out.status[k] = GridStatus::undefined;
            std::fill_n(out.values.begin() + std::ptrdiff_t(k * outputs), outputs, std::complex<double>(nan, nan));
            continue;
        }
        std::copy(v.begin(), v.end(), out.values.begin() + std::ptrdiff_t(k * outputs));
    }
    out.grid = std::move(grid);
    return out;
}

GridField evaluateOnGrid(const Observable& observable, DomainGrid grid) {
    if (observable.stateDimension() != grid.dimension())
        throw std::invalid_argument("The grid lives in a space of another dimension than the observable's state");
    auto out = evaluateOnGrid([&](const double* x) { return observable.evaluate(x); }, observable.outputDimension(), std::move(grid));
    out.complexValued = observable.descriptor().complexValued;
    out.outputNames = observable.outputDimension() == 1 ? std::vector<std::string>{observable.descriptor().name}
                                                        : observable.descriptor().outputNames;
    return out;
}

GridField evaluateOnGrid(const FunctionRepresentation& representation, DomainGrid grid) {
    if (representation.coordinateDimension() != grid.dimension())
        throw std::invalid_argument("The grid and the representation domain differ in dimension");
    auto out = evaluateOnGrid([&](const double* q) { return representation.evaluate(q); }, representation.outputCount(),
                              std::move(grid));
    out.complexValued = true;
    return out;
}

}  // namespace dynsyslab::core
