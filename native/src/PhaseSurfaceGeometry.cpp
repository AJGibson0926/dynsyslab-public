#include "PhaseSurfaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::phase_detail::surface {
double dot(const PhaseState &a, const PhaseState &b, const PhaseSpaceOptions &o) {
    double result = 0;
    for (std::size_t j = 0; j < a.size(); ++j) {
        const double w = o.upper[j] - o.lower[j];
        result += a[j] * b[j] / (w * w);
    }
    return result;
}
double length(const PhaseState &a, const PhaseSpaceOptions &o) {
    return std::sqrt(std::max(0., dot(a, a, o)));
}
PhaseState subtract(PhaseState a, const PhaseState &b) {
    for (std::size_t j = 0; j < a.size(); ++j)
        a[j] -= b[j];
    return a;
}
PhaseState unit(PhaseState a, const PhaseSpaceOptions &o) {
    const double size = length(a, o);
    if (size < 1e-14)
        throw std::runtime_error("Degenerate manifold tangent or mesh edge");
    for (auto &v : a)
        v /= size;
    return a;
}
PhaseState interpolate(const Ring &ring, double position) {
    if (ring.empty() || !std::isfinite(position))
        throw std::runtime_error("Invalid surface ring parameter");
    position = std::fmod(position, double(ring.size()));
    if (position < 0)
        position += ring.size();
    // A tiny negative parameter can round to exactly ring.size() after wrapping.
    if (position >= double(ring.size()))
        position = 0;
    const auto i = std::size_t(std::floor(position)), j = (i + 1) % ring.size();
    const double t = position - i;
    auto x = ring[i];
    for (std::size_t k = 0; k < x.size(); ++k)
        x[k] += (ring[j][k] - x[k]) * t;
    return x;
}
PhaseState tangent(const Ring &ring, double p, const PhaseSpaceOptions &o) {
    const auto r = interpolate(ring, p), left = interpolate(ring, p - 1),
               right = interpolate(ring, p + 1);
    auto a = unit(subtract(r, left), o), b = unit(subtract(right, r), o);
    for (std::size_t j = 0; j < a.size(); ++j)
        a[j] += b[j];
    return unit(a, o);
}
double angle(const PhaseState &a, const PhaseState &b, const PhaseSpaceOptions &o) {
    return std::acos(std::clamp(dot(unit(a, o), unit(b, o), o), -1., 1.));
}
Ring tangentBasis(const PhaseFixedPoint &point, bool stable, const PhaseSpaceOptions &o) {
    Ring basis;
    const auto n = point.state.size();
    for (std::size_t k = 0; k < n; ++k) {
        const auto value = point.eigenvalues[k];
        if (stable ? value.real() >= -o.eigenvalueTolerance : value.real() <= o.eigenvalueTolerance)
            continue;
        if (value.imag() < -o.eigenvalueTolerance)
            continue;
        PhaseState real(n), imaginary(n);
        for (std::size_t j = 0; j < n; ++j) {
            real[j] = point.eigenvectors[j * n + k].real();
            imaginary[j] = point.eigenvectors[j * n + k].imag();
        }
        if (std::abs(value.imag()) > o.eigenvalueTolerance) {
            // One common scale preserves the real/imaginary eigenvector ellipse.
            const double size = std::sqrt(dot(real, real, o) + dot(imaginary, imaginary, o));
            if (size < 1e-14)
                throw std::runtime_error("Degenerate complex surface eigenspace");
            for (std::size_t j = 0; j < n; ++j) {
                real[j] /= size;
                imaginary[j] /= size;
            }
            basis.push_back(real);
            basis.push_back(imaginary);
        } else
            basis.push_back(unit(real, o));
    }
    if (basis.size() != 2)
        throw std::runtime_error(
            "A rank-two real invariant eigenspace is required for surface continuation");
    const double aa = dot(basis[0], basis[0], o), bb = dot(basis[1], basis[1], o),
                 ab = dot(basis[0], basis[1], o);
    if (aa * bb - ab * ab < 1e-8 * aa * bb)
        throw std::runtime_error("Surface eigenvector ellipse is ill-conditioned");
    return basis;
}
} // namespace dynsyslab::core::phase_detail::surface
