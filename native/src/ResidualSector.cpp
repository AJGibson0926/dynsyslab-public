// Exact two-variable LP for the residual sector:
//   minimize w_x c_x + w_u c_u  s.t.  c_x a_k + c_u b_k >= rho_k,  c >= 0,
// with a = ||z||, b = ||u||, rho = ||r||. For fixed t = c_u the best c_x is the
// upper envelope of the lines (rho_k - t b_k) / a_k and 0, so the objective is
// convex piecewise linear and its minimum lies at a breakpoint. Error structure:
// SafEDMD (Strasser, Worthmann et al., arXiv:2402.03145). Scenario bound: Campi &
// Garatti, SIAM J. Optim. 19 (2008) 1211-1230.
#include <dynsyslab/core/ResidualSector.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Line { double slope, intercept; };

double norm(const std::vector<double>& v, std::size_t start, std::size_t count) {
    double sum = 0.0;
    for (std::size_t i = start; i < start + count; ++i) sum += v[i] * v[i];
    return std::sqrt(sum);
}

void checkArrays(const std::vector<double>& z, const std::vector<double>& u, const std::vector<double>& r,
    std::size_t samples, std::size_t n, std::size_t m) {
    if (samples == 0 || n == 0 || m == 0) throw std::invalid_argument("Residual sector needs samples and dimensions");
    if (z.size() != samples * n || u.size() != samples * m || r.size() != samples * n)
        throw std::invalid_argument("Residual sector arrays have inconsistent sizes");
    for (const auto* values : {&z, &u, &r})
        for (double v : *values)
            if (!std::isfinite(v)) throw std::invalid_argument("Residual sector data must be finite");
}

// Upper envelope of lines sorted by increasing slope.
std::vector<Line> upperEnvelope(std::vector<Line> lines) {
    std::sort(lines.begin(), lines.end(), [](const Line& a, const Line& b) {
        return a.slope != b.slope ? a.slope < b.slope : a.intercept < b.intercept;
    });
    std::vector<Line> hull;
    for (const auto& line : lines) {
        if (!hull.empty() && hull.back().slope == line.slope) hull.pop_back();
        while (hull.size() >= 2) {
            const auto& p = hull[hull.size() - 2];
            const auto& q = hull.back();
            // q is redundant when line overtakes p no later than q does.
            if ((line.intercept - p.intercept) * (q.slope - p.slope) >=
                (q.intercept - p.intercept) * (line.slope - p.slope)) hull.pop_back();
            else break;
        }
        hull.push_back(line);
    }
    return hull;
}

}  // namespace

double scenarioViolationBound(std::size_t samples, std::size_t d, double beta) {
    if (samples == 0 || d == 0 || d > samples || !(beta > 0.0 && beta < 1.0))
        throw std::invalid_argument("Invalid scenario bound arguments");
    const double N = double(samples);
    const auto logTail = [&](double eps) {  // log of the binomial tail
        double best = -std::numeric_limits<double>::infinity();
        std::vector<double> terms;
        for (std::size_t i = 0; i < d; ++i) {
            const double k = double(i);
            const double value = std::lgamma(N + 1) - std::lgamma(k + 1) - std::lgamma(N - k + 1) +
                k * std::log(eps) + (N - k) * std::log1p(-eps);
            terms.push_back(value);
            best = std::max(best, value);
        }
        double sum = 0.0;
        for (double value : terms) sum += std::exp(value - best);
        return best + std::log(sum);
    };
    const double target = std::log(beta);
    double low = 1e-15, high = 1.0 - 1e-15;
    if (logTail(high) > target) return 1.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (low + high);
        (logTail(mid) > target ? low : high) = mid;
    }
    return high;  // Rounded toward the conservative side.
}

ResidualSectorBound estimateResidualSector(const std::vector<double>& z, const std::vector<double>& u,
    const std::vector<double>& residual, std::size_t samples, std::size_t n, std::size_t m,
    const ResidualSectorOptions& options) {
    checkArrays(z, u, residual, samples, n, m);
    if (!(options.inflation >= 1.0) || !std::isfinite(options.inflation) || !(options.stateWeight >= 0.0) ||
        !(options.inputWeight >= 0.0) || !(options.zeroTolerance >= 0.0))
        throw std::invalid_argument("Invalid residual sector options");
    // Data-dependent weights would change the program with the sample and void the scenario theorem.
    if (options.scenario && !(options.stateWeight > 0.0 && options.inputWeight > 0.0))
        throw std::invalid_argument("A scenario bound needs state and input weights fixed before sampling");

    std::vector<double> a(samples), b(samples), rho(samples);
    double scale = 0.0, meanA = 0.0, meanB = 0.0, maxRho = 0.0;
    for (std::size_t k = 0; k < samples; ++k) {
        a[k] = norm(z, k * n, n); b[k] = norm(u, k * m, m); rho[k] = norm(residual, k * n, n);
        scale = std::max({scale, a[k], b[k]});
        meanA += a[k] / double(samples); meanB += b[k] / double(samples);
        maxRho = std::max(maxRho, rho[k]);
    }
    const double zero = options.zeroTolerance * std::max(scale, 1.0);

    ResidualSectorBound out;
    out.samples = samples;
    out.maximumResidual = maxRho;
    out.inflation = options.inflation;
    out.guarantee = "data envelope over the supplied samples; not a probabilistic bound";

    std::vector<Line> lines{{0.0, 0.0}};
    double lowerT = 0.0;
    for (std::size_t k = 0; k < samples; ++k) {
        if (rho[k] <= zero) continue;
        if (a[k] > zero) lines.push_back({-b[k] / a[k], rho[k] / a[k]});
        else if (b[k] > zero) lowerT = std::max(lowerT, rho[k] / b[k]);
        else {
            out.message = "nonzero residual at z = 0, u = 0: the data violates the zero-target convention";
            return out;
        }
    }
    const auto hull = upperEnvelope(std::move(lines));
    const double wx = options.stateWeight > 0.0 ? options.stateWeight : (meanA > 0.0 ? meanA : 1.0);
    const double wu = options.inputWeight > 0.0 ? options.inputWeight : (meanB > 0.0 ? meanB : 1.0);
    const auto envelope = [&](double t) {
        double value = 0.0;
        for (const auto& line : hull) value = std::max(value, line.slope * t + line.intercept);
        return value;
    };
    double bestT = lowerT, bestX = envelope(lowerT), bestValue = wx * bestX + wu * bestT;
    for (std::size_t i = 0; i + 1 < hull.size(); ++i) {
        const double t = (hull[i].intercept - hull[i + 1].intercept) / (hull[i + 1].slope - hull[i].slope);
        if (!(t > lowerT) || !std::isfinite(t)) continue;
        const double x = std::max(hull[i].slope * t + hull[i].intercept, 0.0);
        const double value = wx * x + wu * t;
        if (value < bestValue) { bestValue = value; bestT = t; bestX = x; }
    }
    // Recompute c_x at the chosen t so the constraints hold, then round outward.
    bestX = envelope(bestT);
    const double outward = 1.0 + 8.0 * std::numeric_limits<double>::epsilon();
    out.cState = bestX * outward * options.inflation;
    out.cInput = bestT * outward * options.inflation;
    out.feasible = true;
    if (options.scenario) {
        if (options.scenario->declaration.empty())
            throw std::invalid_argument("A scenario bound needs a sampling declaration");
        const double confidence = options.scenario->confidence;
        if (!(confidence > 0.0 && confidence < 1.0)) throw std::invalid_argument("Scenario confidence must be in (0, 1)");
        if (samples < 2) throw std::invalid_argument("A scenario bound needs at least two samples");
        out.violationProbability = scenarioViolationBound(samples, 2, 1.0 - confidence);
        out.confidence = confidence;
        out.guarantee = "scenario bound, conditional on the caller's declaration: " + options.scenario->declaration;
    }
    return out;
}

ResidualSectorBound estimateResidualSector(const BilinearKrom& model, const BilinearDataset& data,
    const ResidualSectorOptions& options) {
    const auto residual = bilinearResiduals(model, data);
    return estimateResidualSector(data.z, data.u, residual, data.samples, model.stateDimension,
        model.inputDimension, options);
}

ResidualSectorCoverage residualSectorCoverage(const ResidualSectorBound& bound, const std::vector<double>& z,
    const std::vector<double>& u, const std::vector<double>& residual, std::size_t samples,
    std::size_t n, std::size_t m) {
    checkArrays(z, u, residual, samples, n, m);
    if (!bound.feasible) throw std::invalid_argument("Coverage needs a feasible sector bound");
    ResidualSectorCoverage out;
    out.samples = samples;
    for (std::size_t k = 0; k < samples; ++k) {
        const double excess = norm(residual, k * n, n) - bound.cState * norm(z, k * n, n) - bound.cInput * norm(u, k * m, m);
        if (excess > 0.0) { ++out.violations; out.worstExcess = std::max(out.worstExcess, excess); }
    }
    out.fraction = 1.0 - double(out.violations) / double(samples);
    return out;
}

}  // namespace dynsyslab::core
