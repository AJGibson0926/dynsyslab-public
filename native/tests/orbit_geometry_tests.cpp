// Orbit geometry for injectivity audits (core/OrbitGeometry.hpp) on curves with known answers.
#include <dynsyslab/core/OrbitGeometry.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
constexpr double pi = 3.141592653589793;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void circle(double cx, double cy, double r, std::size_t n, std::vector<double>& x, std::vector<double>& y) {
    x.clear(); y.clear();
    for (std::size_t k = 0; k <= n; ++k) {
        const double t = 2 * pi * double(k) / double(n);
        x.push_back(cx + r * std::cos(t)); y.push_back(cy + r * std::sin(t));
    }
}

}  // namespace

int main() {
    try {
        // Two unit circles one radius apart cross at (1/2, +-sqrt(3)/2).
        std::vector<double> ax, ay, bx, by;
        circle(0, 0, 1, 2000, ax, ay);
        circle(1, 0, 1, 2000, bx, by);
        const auto crossings = c::findCurveIntersections(ax, ay, bx, by, 1e-6);
        require(crossings.size() == 2, "two crossings of the offset circles");
        for (const auto& x : crossings)
            require(std::abs(x.x - 0.5) < 1e-5 && std::abs(std::abs(x.y) - std::sqrt(3.0) / 2) < 1e-5, "crossing points");
        // Concentric circles never cross; their nearest separation is the radius gap.
        circle(0, 0, 2, 2000, bx, by);
        require(c::findCurveIntersections(ax, ay, bx, by).empty(), "concentric circles do not cross");
        std::vector<double> a, b;
        for (std::size_t k = 0; k < ax.size(); ++k) { a.push_back(ax[k]); a.push_back(ay[k]); b.push_back(bx[k]); b.push_back(by[k]); }
        require(std::abs(c::nearestSeparation(a, b, 2) - 1) < 1e-9, "nearest separation of concentric circles");
        require(std::abs(c::polygonArea(bx, by) - 4 * pi) < 1e-4, "enclosed area of a circle");

        // Peaks of a sinusoid with a small ripple: one per period once prominence filters the ripple.
        std::vector<double> y, path, spiral;
        const double dt = 0.01;
        for (std::size_t k = 0; k < 6000; ++k) {
            const double t = double(k) * dt;
            y.push_back(std::cos(t) + 0.01 * std::cos(37 * t));
            path.push_back(std::cos(t)); path.push_back(std::sin(t));
            const double r = std::exp(-0.05 * t);
            spiral.push_back(r * std::cos(t)); spiral.push_back(r * std::sin(t));
        }
        const auto peaks = c::findPeaks(y, std::size_t(0.5 / dt), 0.15 * 2);
        require(peaks.size() == 9, "one peak per period");  // t = 2 pi k, k = 1..9 within 60
        for (std::size_t k = 0; k < peaks.size(); ++k) require(std::abs(double(peaks[k]) * dt - 2 * pi * double(k + 1)) < 0.05, "peak times");
        require(c::peakClosure(path, 2, peaks) < 0.01, "a periodic orbit closes at its peaks");
        require(c::peakClosure(spiral, 2, peaks) > 0.05, "a decaying spiral does not close");
        const auto cycle = c::representativeCycle(path, 2, peaks, 400);
        std::vector<double> cx, cy;
        for (std::size_t s = 0; s < 400; ++s) { cx.push_back(cycle[2 * s]); cy.push_back(cycle[2 * s + 1]); }
        require(std::abs(c::polygonArea(cx, cy) - pi) < 1e-3, "representative cycle encloses the unit disc");
        // The Lissajous curve (sin t, sin 2t) revisits the origin at t = pi and 2 pi.
        std::vector<double> ft, fx, fy;
        for (int k = 0; k <= 20000; ++k) { const double t = 0.5 + 0.0005 * k; ft.push_back(t); fx.push_back(std::sin(t)); fy.push_back(std::sin(2 * t)); }
        const auto visit = c::refineVisits(ft, fx, fy, 3.0, 3.3, 6.1, 6.4, 1, 1, 1e-10);
        require(visit.matched && std::abs(visit.first - pi) < 1e-6 && std::abs(visit.second - 2 * pi) < 1e-6, "refined visits of a crossing");
        const auto matching = c::matchVisits({{1, 5}, {2, 6}}, {{2.05, 6.02}, {1.01, 4.99}, {9, 9}}, 0.1, 0.1);
        require(matching.matches.size() == 2 && matching.matches[0] == std::make_pair(std::size_t(0), std::size_t(1))
                && matching.matches[1] == std::make_pair(std::size_t(1), std::size_t(0)) && matching.unmatchedSecond == std::vector<std::size_t>{2},
                "one-to-one visit matching");
        std::cout << "orbit geometry tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
