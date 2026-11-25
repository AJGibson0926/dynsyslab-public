// Pull-back, push-forward, flow maps and the Koopman action on grids (CoordinateMaps.hpp).
#include <dynsyslab/core/CoordinateMaps.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// One classical RK4 integration of the unforced field over `time`, in equal steps <= maximumStep.
bool flow(const ODESystem& system, double* x, double time, double maximumStep) {
    const std::size_t n = system.descriptor().stateDimension;
    const auto steps = std::max<std::size_t>(1, std::size_t(std::ceil(std::abs(time) / maximumStep)));
    const double h = time / double(steps);
    std::vector<double> k1(n), k2(n), k3(n), k4(n), y(n);
    for (std::size_t s = 0; s < steps; ++s) {
        system.rhs(0, x, k1.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k1[i];
        system.rhs(0, y.data(), k2.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k2[i];
        system.rhs(0, y.data(), k3.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h * k3[i];
        system.rhs(0, y.data(), k4.data());
        for (std::size_t i = 0; i < n; ++i) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        if (system.invalid(0, x)) return false;
        for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
    }
    return true;
}

double determinant(std::vector<double> a, std::size_t n) {  // Gaussian elimination with pivoting
    double det = 1;
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r) if (std::abs(a[r * n + c]) > std::abs(a[p * n + c])) p = r;
        if (a[p * n + c] == 0) return 0;
        if (p != c) { for (std::size_t k = 0; k < n; ++k) std::swap(a[p * n + k], a[c * n + k]); det = -det; }
        det *= a[c * n + c];
        for (std::size_t r = c + 1; r < n; ++r) {
            const double f = a[r * n + c] / a[c * n + c];
            for (std::size_t k = c; k < n; ++k) a[r * n + k] -= f * a[c * n + k];
        }
    }
    return det;
}

}  // namespace

Observable pullBack(const Observable& function, const Observable& map) {
    if (map.outputDimension() != function.stateDimension())
        throw std::invalid_argument("pull-back: the map's outputs are not the function's coordinates");
    const std::size_t outputs = function.outputDimension(), q = function.stateDimension(), x = map.stateDimension();
    auto descriptor = function.descriptor();
    descriptor.name += "_pullback";
    descriptor.label = descriptor.label + " pulled back by " + map.descriptor().label;
    Observable::Evaluator evaluate = [function, map](const double* state) {
        const auto image = map.evaluateReal(state);
        return function.evaluate(image.data());
    };
    Observable::Jacobian jacobian;
    if (function.hasAnalyticJacobian() && map.hasAnalyticJacobian())
        jacobian = [function, map, outputs, q, x](const double* state) {
            const auto image = map.evaluateReal(state);
            const auto outer = function.jacobian(image.data(), ObservableJacobianBackend::analytic);  // outputs x q
            const auto inner = map.jacobian(state, ObservableJacobianBackend::analytic);              // q x x
            Observable::Values out(outputs * x);
            for (std::size_t r = 0; r < outputs; ++r)
                for (std::size_t a = 0; a < q; ++a)
                    for (std::size_t i = 0; i < x; ++i) out[r * x + i] += outer[r * q + a] * inner[a * x + i];
            return out;
        };
    return Observable(std::move(descriptor), x, outputs, std::move(evaluate), std::move(jacobian));
}

Observable flowMap(std::shared_ptr<const ODESystem> system, double time, double maximumStep) {
    if (!system || !std::isfinite(time) || !(maximumStep > 0)) throw std::invalid_argument("flow map needs a system, a finite time and a positive step");
    const auto n = system->descriptor().stateDimension;
    std::vector<std::string> names;
    for (const auto& s : system->descriptor().stateNames) names.push_back(s);
    if (names.size() != n) { names.clear(); for (std::size_t i = 0; i < n; ++i) names.push_back("x_" + std::to_string(i)); }
    ObservableDescriptor descriptor{"flow", "Flow map over t = " + std::to_string(time), "", system->descriptor().preferredCoordinates, names, false};
    return Observable(std::move(descriptor), n, n, [system, time, maximumStep, n](const double* x) {
        std::vector<double> y(x, x + n);
        if (!flow(*system, y.data(), time, maximumStep)) y.assign(n, nan);
        return Observable::Values(y.begin(), y.end());
    });
}

std::vector<GridField> koopmanAction(const Observable& function, const ODESystem& system, const DomainGrid& grid,
                                     double timeStep, std::size_t steps, double maximumStep) {
    const auto n = system.descriptor().stateDimension;
    if (grid.dimension() != n || function.stateDimension() != n)
        throw std::invalid_argument("Koopman action: the grid and the function must live in the system's state space");
    // Flow every grid point once per step and evaluate g along the way.
    const auto points = grid.points();
    std::vector<double> state = points;
    std::vector<bool> alive(grid.size(), true);
    std::vector<GridField> out;
    for (std::size_t k = 0; k <= steps; ++k) {
        if (k)
            for (std::size_t p = 0; p < grid.size(); ++p)
                if (alive[p]) alive[p] = flow(system, state.data() + p * n, timeStep, maximumStep);
        GridField field;
        field.grid = grid;
        field.outputs = function.outputDimension();
        field.values.resize(grid.size() * field.outputs);
        field.status.resize(grid.size());
        for (std::size_t p = 0; p < grid.size(); ++p) {  // g at the flowed state of each grid point
            std::vector<std::complex<double>> v;
            bool ok = alive[p];
            if (ok) {
                try { v = function.evaluate(state.data() + p * n); } catch (const std::exception&) { ok = false; }
                ok = ok && std::all_of(v.begin(), v.end(), [](auto z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); });
            }
            field.status[p] = ok ? GridStatus::valid : GridStatus::undefined;
            for (std::size_t r = 0; r < field.outputs; ++r)
                field.values[p * field.outputs + r] = ok ? v[r] : std::complex<double>(nan, nan);
        }
        field.complexValued = function.descriptor().complexValued;
        field.outputNames = function.outputDimension() == 1 ? std::vector<std::string>{function.descriptor().name}
                                                            : function.descriptor().outputNames;
        out.push_back(std::move(field));
    }
    return out;
}

GridField pushForward(const Observable& function, const Observable& map, const DomainGrid& sourceGrid,
                      const DomainGrid& target, PushForwardOptions options) {
    sourceGrid.validate();
    auto source = sourceGrid;
    if (options.oversampling > 1)
        for (auto& axis : source.axes) axis.points = (axis.points - 1) * options.oversampling + 1;
    target.validate();
    if (function.stateDimension() != source.dimension() || map.stateDimension() != source.dimension()
        || map.outputDimension() != target.dimension())
        throw std::invalid_argument("push-forward: the function and map must act on the source grid, the map into the target grid");
    const std::size_t outputs = function.outputDimension(), cells = target.size(), d = source.dimension();
    std::vector<std::complex<double>> sum(cells * outputs);
    std::vector<std::size_t> count(cells);
    std::vector<std::size_t> cellOf(source.size(), cells);
    std::vector<std::complex<double>> samples(source.size() * outputs);
    std::vector<double> x(d);
    double scale = 0;
    for (std::size_t s = 0; s < source.size(); ++s) {
        source.point(s, x.data());
        try {
            const auto image = map.evaluateReal(x.data());
            const auto value = function.evaluate(x.data());
            const auto cell = target.locate(image.data());
            if (cell >= cells) continue;
            bool finite = true;
            for (const auto v : value) finite = finite && std::isfinite(v.real()) && std::isfinite(v.imag());
            if (!finite) continue;
            cellOf[s] = cell;
            ++count[cell];
            for (std::size_t r = 0; r < outputs; ++r) {
                samples[s * outputs + r] = value[r];
                sum[cell * outputs + r] += value[r];
                scale = std::max(scale, std::abs(value[r]));
            }
        } catch (const std::exception&) {
        }
    }
    GridField out;
    out.grid = target;
    out.outputs = outputs;
    out.complexValued = function.descriptor().complexValued;
    out.outputNames = outputs == 1 ? std::vector<std::string>{function.descriptor().name} : function.descriptor().outputNames;
    out.values.assign(cells * outputs, {nan, nan});
    out.spread.assign(cells * outputs, nan);
    out.status.assign(cells, GridStatus::noPreimage);
    out.counts = count;
    for (std::size_t c = 0; c < cells; ++c)
        if (count[c])
            for (std::size_t r = 0; r < outputs; ++r) {
                out.values[c * outputs + r] = sum[c * outputs + r] / double(count[c]);
                out.spread[c * outputs + r] = 0;
            }
    for (std::size_t s = 0; s < source.size(); ++s) {
        const auto c = cellOf[s];
        if (c >= cells) continue;
        for (std::size_t r = 0; r < outputs; ++r)
            out.spread[c * outputs + r] = std::max(out.spread[c * outputs + r], std::abs(samples[s * outputs + r] - out.values[c * outputs + r]));
    }
    const double limit = options.tolerance * std::max(scale, 1e-300);
    for (std::size_t c = 0; c < cells; ++c) {
        if (!count[c]) continue;
        bool consistent = true;
        for (std::size_t r = 0; r < outputs; ++r) consistent = consistent && out.spread[c * outputs + r] <= limit;
        out.status[c] = consistent ? GridStatus::valid : GridStatus::conflict;
    }
    return out;
}

GridField roundTripError(const Observable& function, const Observable& map, const DomainGrid& source,
                         const DomainGrid& target, PushForwardOptions options) {
    const auto forward = pushForward(function, map, source, target, options);
    const auto outputs = function.outputDimension();
    auto out = evaluateOnGrid([&](const double* x) {
        const auto image = map.evaluateReal(x);
        const auto cell = target.locate(image.data());
        const auto value = function.evaluate(x);
        std::vector<std::complex<double>> error(outputs, {nan, nan});
        if (cell >= target.size() || forward.status[cell] == GridStatus::noPreimage) return error;
        for (std::size_t r = 0; r < outputs; ++r) error[r] = std::abs(forward.value(cell, r) - value[r]);
        return error;
    }, outputs, source);
    // Carry the target status back: where the cell is a conflict or has no preimage, so is x.
    std::vector<double> x(source.dimension());
    for (std::size_t s = 0; s < source.size(); ++s) {
        source.point(s, x.data());
        try {
            const auto cell = target.locate(map.evaluateReal(x.data()).data());
            out.status[s] = cell >= target.size() ? GridStatus::noPreimage
                          : forward.status[cell] == GridStatus::valid ? out.status[s] : forward.status[cell];
        } catch (const std::exception&) {
            out.status[s] = GridStatus::undefined;
        }
    }
    out.outputNames = forward.outputNames;
    return out;
}

GridField mapVolumeFactor(const Observable& map, const DomainGrid& grid) {
    const std::size_t m = map.outputDimension(), n = map.stateDimension();
    if (n != grid.dimension()) throw std::invalid_argument("volume factor: the grid lives in another space than the map");
    auto out = evaluateOnGrid([&](const double* x) {
        const auto J = map.jacobian(x);  // m x n
        std::vector<double> gram(m * m);
        for (std::size_t a = 0; a < m; ++a)
            for (std::size_t b = 0; b < m; ++b)
                for (std::size_t i = 0; i < n; ++i) gram[a * m + b] += (J[a * n + i] * std::conj(J[b * n + i])).real();
        return std::vector<std::complex<double>>{std::sqrt(std::max(0.0, determinant(gram, m)))};
    }, 1, grid);
    out.outputNames = {"volume_factor"};
    return out;
}

}  // namespace dynsyslab::core
