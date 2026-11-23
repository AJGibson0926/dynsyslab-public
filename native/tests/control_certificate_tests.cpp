// Identity tests for the certified bilinear control tools. No research data.
#include <dynsyslab/core/BilinearKrom.hpp>
#include <dynsyslab/core/CertifiedDesign.hpp>
#include <dynsyslab/core/DenseSdp.hpp>
#include <dynsyslab/core/ResidualSector.hpp>
#include <dynsyslab/math/ControlCertificateQuery.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

using namespace dynsyslab::core;
namespace {

int checks = 0;
void check(bool ok, const std::string& what) { ++checks; if (!ok) throw std::runtime_error(what); }
void near(double a, double b, double tol, const std::string& what) {
    check(std::abs(a - b) <= tol, what + " (" + std::to_string(a) + " vs " + std::to_string(b) + ")");
}
template <class F> void rejects(F f, const std::string& what) {
    bool caught = false;
    try { (void)f(); } catch (const std::exception&) { caught = true; }
    check(caught, what);
}

// P = [p11 p12; p12 p22]: I >= P >= d I and P - A'PA >= d I.
SdpProblem lyapunovProblem(const std::array<double, 4>& a, double d) {
    SdpProblem p;
    p.variables = 3;
    const std::array<std::array<double, 4>, 3> basis{{{1, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 1}}};
    const auto lyap = [&](const std::array<double, 4>& P) {
        std::vector<double> out(4);
        for (std::size_t i = 0; i < 2; ++i) for (std::size_t j = 0; j < 2; ++j) {
            double s = P[i * 2 + j];
            for (std::size_t k = 0; k < 2; ++k) for (std::size_t l = 0; l < 2; ++l) s -= a[k * 2 + i] * P[k * 2 + l] * a[l * 2 + j];
            out[i * 2 + j] = s;
        }
        return out;
    };
    LmiBlock lower{2, {-d, 0, 0, -d}, {}}, upper{2, {1, 0, 0, 1}, {}}, decrease{2, {-d, 0, 0, -d}, {}};
    for (const auto& e : basis) {
        lower.coefficients.emplace_back(e.begin(), e.end());
        upper.coefficients.push_back({-e[0], -e[1], -e[2], -e[3]});
        decrease.coefficients.push_back(lyap(e));
    }
    p.blocks = {lower, upper, decrease};
    return p;
}

void lyapunovFeasibility() {
    const auto stable = lyapunovProblem({0.9, 0.2, 0.0, 0.8}, 1e-3);
    const auto r = solveDenseSdp(stable);
    check(r.status == SdpStatus::converged, "stable Lyapunov LMI should converge");
    // The solver enforced F >= 1e-3 I, so the unshifted LMI must clear a positive margin.
    check(verifyLmiCertificate(lyapunovProblem({0.9, 0.2, 0.0, 0.8}, 0.0), r.x, 5e-4).certified, "stable Lyapunov solution should certify");
    const auto unstable = lyapunovProblem({1.1, 0.2, 0.0, 0.8}, 1e-3);
    const auto u = solveDenseSdp(unstable);
    check(u.status == SdpStatus::primalInfeasible, std::string("unstable Lyapunov LMI: ") + sdpStatusName(u.status));
    check(u.x.empty(), "infeasible result carries no solution");
}

LmiBlock block(std::size_t d, std::vector<double> constant, std::vector<std::vector<double>> terms) {
    LmiBlock b; b.dimension = d; b.constant = std::move(constant); b.coefficients = std::move(terms); return b;
}

void knownValues() {
    // min x s.t. [x 1; 1 x] >= 0 has value 1.
    SdpProblem a;
    a.variables = 1; a.objective = {1.0};
    a.blocks.push_back(block(2, {0, 1, 1, 0}, {{1, 0, 0, 1}}));
    const auto ra = solveDenseSdp(a);
    check(ra.status == SdpStatus::converged, "known SDP converges");
    near(ra.objective, 1.0, 1e-5, "known SDP value");
    // max t s.t. M - t I >= 0 is lambda_min([2 1; 1 2]) = 1.
    SdpProblem b;
    b.variables = 1; b.objective = {-1.0};
    b.blocks.push_back(block(2, {2, 1, 1, 2}, {{-1, 0, 0, -1}}));
    const auto rb = solveDenseSdp(b);
    check(rb.status == SdpStatus::converged, "eigenvalue SDP converges");
    near(rb.x[0], 1.0, 1e-5, "minimum eigenvalue by SDP");
    // Unbounded: min -x s.t. x >= 0.
    SdpProblem c;
    c.variables = 1; c.objective = {-1.0};
    c.blocks.push_back(block(1, {0}, {{1}}));
    check(solveDenseSdp(c).status == SdpStatus::dualInfeasible, "unbounded SDP is reported");
    // min trace P s.t. P - A'PA >= I is attained by the discrete Lyapunov solution.
    const std::array<double, 4> stableA{0.7, 0.3, -0.2, 0.5};
    SdpProblem ly;
    ly.variables = 3; ly.objective = {1.0, 0.0, 1.0};
    ly.blocks = {lyapunovProblem(stableA, 1.0).blocks[2]};
    const auto rl = solveDenseSdp(ly);
    check(rl.status == SdpStatus::converged, "Lyapunov trace SDP converges");
    std::array<double, 4> term{1, 0, 0, 1}, sum{0, 0, 0, 0};
    for (int k = 0; k < 400; ++k) {  // P = sum_k (A')^k A^k
        for (std::size_t i = 0; i < 4; ++i) sum[i] += term[i];
        std::array<double, 4> next{};
        for (std::size_t i = 0; i < 2; ++i) for (std::size_t j = 0; j < 2; ++j)
            for (std::size_t p = 0; p < 2; ++p) for (std::size_t q = 0; q < 2; ++q)
                next[i * 2 + j] += stableA[p * 2 + i] * term[p * 2 + q] * stableA[q * 2 + j];
        term = next;
    }
    near(rl.objective, sum[0] + sum[3], 2e-5, "Lyapunov trace value");
    near(rl.x[1], sum[1], 2e-5, "Lyapunov off-diagonal");
    // The verifier rejects a point outside the cone and asymmetry in the data.
    check(!verifyLmiCertificate(a, {0.5}, 0.0).certified, "verifier rejects infeasible point");
    check(verifyLmiCertificate(a, {1.5}, 0.1).certified, "verifier accepts interior point");
    check(!verifyLmiCertificate(a, {1.5}, 0.6).certified, "verifier enforces the margin");
    auto bad = a; bad.blocks[0].constant = {0, 1, 0, 0};
    check(!verifyLmiCertificate(bad, {2.0}, 0.0).certified, "verifier rejects asymmetric LMI");
    rejects([&] { return solveDenseSdp(bad); }, "solver rejects asymmetric LMI");
}

BilinearKrom plant() {
    BilinearKrom m;
    m.stateDimension = 2; m.inputDimension = 1; m.sampleTime = 0.1;
    m.A = {1.02, 0.1, 0.0, 0.98};
    m.B0 = {0.0, 0.1};
    m.B = {{0.0, 0.0, 0.02, 0.0}};
    m.provenance = "unit test";
    return m;
}

BilinearDataset sampled(const BilinearKrom& m, std::size_t samples, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    BilinearDataset d;
    d.samples = samples;
    for (std::size_t k = 0; k < samples; ++k) {
        const std::vector<double> z{uniform(rng), uniform(rng)}, u{uniform(rng)};
        const auto next = stepBilinearKrom(m, z, u);
        d.z.insert(d.z.end(), z.begin(), z.end());
        d.u.insert(d.u.end(), u.begin(), u.end());
        d.zNext.insert(d.zNext.end(), next.begin(), next.end());
    }
    return d;
}

void bilinearRecovery() {
    const auto truth = plant();
    const auto fit = fitBilinearKrom(sampled(truth, 60, 3), 2, 1);
    check(fit.identifiable && fit.rank == 5, "bilinear design has full rank");
    for (std::size_t i = 0; i < 4; ++i) near(fit.model.A[i], truth.A[i], 1e-11, "A recovered");
    for (std::size_t i = 0; i < 2; ++i) near(fit.model.B0[i], truth.B0[i], 1e-11, "B0 recovered");
    for (std::size_t i = 0; i < 4; ++i) near(fit.model.B[0][i], truth.B[0][i], 1e-11, "B1 recovered");
    check(fit.stats.maximum < 1e-11, "noise-free residual vanishes");

    const auto lin = linearizeBilinearKrom(truth, {0.3, -0.2}, {0.5});
    near(lin.A[2], 0.01, 1e-14, "frozen A");
    near(lin.B[1], 0.1 + 0.02 * 0.3, 1e-14, "frozen B");
    const auto next = stepBilinearKrom(truth, {0.3, -0.2}, {0.5});
    near(lin.A[2] * 0.3 + lin.A[3] * -0.2 + lin.B[1] * 0.5 + lin.offset[1], next[1], 1e-14, "linearization is exact at the point");

    const auto text = json::dump(bilinearKromToJson(truth));
    const auto back = bilinearKromFromJson(json::parse(text));
    check(back.A == truth.A && back.B0 == truth.B0 && back.B == truth.B, "JSON round trip is exact");
    auto deficient = sampled(truth, 60, 4);
    std::fill(deficient.u.begin(), deficient.u.end(), 0.0);
    rejects([&] { return fitBilinearKrom(deficient, 2, 1); }, "rank-deficient fit is rejected");
}

void sectorEstimation() {
    // ||r|| = 0.03 ||z|| + 0.01 ||u|| exactly: the LP must recover a tight envelope.
    std::mt19937 rng(9);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::vector<double> z, u, r;
    const std::size_t samples = 400;
    for (std::size_t k = 0; k < samples; ++k) {
        const double a = uniform(rng), b = uniform(rng), v = uniform(rng), angle = 3.0 * uniform(rng);
        const double size = 0.03 * std::hypot(a, b) + 0.01 * std::abs(v);
        z.insert(z.end(), {a, b}); u.push_back(v);
        r.insert(r.end(), {size * std::cos(angle), size * std::sin(angle)});
    }
    const auto bound = estimateResidualSector(z, u, r, samples, 2, 1);
    check(bound.feasible, "sector is feasible");
    near(bound.cState, 0.03, 1e-9, "c_x recovered");
    near(bound.cInput, 0.01, 1e-9, "c_u recovered");
    check(residualSectorCoverage(bound, z, u, r, samples, 2, 1).violations == 0, "in-sample coverage is complete");
    check(!bound.violationProbability, "no probability without a sampling declaration");

    ResidualSectorOptions declared;
    declared.scenario = ScenarioAssumption{"i.i.d. uniform test samples", 0.99};
    rejects([&] { return estimateResidualSector(z, u, r, samples, 2, 1, declared); }, "scenario needs pre-set weights");
    declared.stateWeight = 1.0; declared.inputWeight = 1.0;
    const auto scenario = estimateResidualSector(z, u, r, samples, 2, 1, declared);
    check(scenario.violationProbability && *scenario.violationProbability > 0.0 && *scenario.violationProbability < 0.05, "scenario bound");
    const double eps = *scenario.violationProbability, N = double(samples);
    near(std::pow(1 - eps, N) + N * eps * std::pow(1 - eps, N - 1), 0.01, 1e-6, "scenario bound solves the binomial tail");

    // The envelope LP matches brute-force vertex enumeration on irregular data.
    for (unsigned seed = 0; seed < 25; ++seed) {
        std::mt19937 local(100 + seed);
        std::vector<double> zs, us, rs;
        const std::size_t count = 24;
        for (std::size_t k = 0; k < count; ++k) {
            const double z1 = uniform(local), z2 = uniform(local), v = uniform(local), r1 = uniform(local), r2 = uniform(local);
            zs.insert(zs.end(), {z1, z2}); us.push_back(v);
            rs.insert(rs.end(), {0.1 * r1, 0.1 * r2});
        }
        ResidualSectorOptions weights; weights.stateWeight = 0.7; weights.inputWeight = 1.3;
        const auto lp = estimateResidualSector(zs, us, rs, count, 2, 1, weights);
        std::vector<std::array<double, 3>> rows;
        for (std::size_t k = 0; k < count; ++k)
            rows.push_back({std::hypot(zs[2 * k], zs[2 * k + 1]), std::abs(us[k]), std::hypot(rs[2 * k], rs[2 * k + 1])});
        const auto admissible = [&](double x, double y) {
            if (x < -1e-12 || y < -1e-12) return false;
            for (const auto& row : rows) if (x * row[0] + y * row[1] < row[2] - 1e-11) return false;
            return true;
        };
        double best = HUGE_VAL;
        const auto consider = [&](double x, double y) { if (admissible(x, y)) best = std::min(best, 0.7 * x + 1.3 * y); };
        for (const auto& row : rows) { consider(row[2] / row[0], 0.0); consider(0.0, row[2] / row[1]); }
        for (std::size_t i = 0; i < count; ++i) for (std::size_t j = 0; j < i; ++j) {
            const double det = rows[i][0] * rows[j][1] - rows[j][0] * rows[i][1];
            if (std::abs(det) < 1e-12) continue;
            consider((rows[i][2] * rows[j][1] - rows[j][2] * rows[i][1]) / det, (rows[i][0] * rows[j][2] - rows[j][0] * rows[i][2]) / det);
        }
        near(0.7 * lp.cState + 1.3 * lp.cInput, best, 1e-9, "sector LP is optimal");
    }

    const auto offset = estimateResidualSector(std::vector<double>{0.0, 0.0}, std::vector<double>{0.0},
        std::vector<double>{0.1, 0.0}, 1, 2, 1);
    check(!offset.feasible, "residual at the target is rejected");
}

ResidualSectorBound sectorOf(double cx, double cu) {
    ResidualSectorBound s; s.feasible = true; s.cState = cx; s.cInput = cu; return s;
}

void certifiedStabilization() {
    const auto model = plant();
    const auto sector = sectorOf(0.004, 0.002);
    CertifiedDesignSpec spec;
    spec.inputBounds = {1.0};
    spec.stateRadius = 1.0;
    spec.contraction = 0.995;
    const auto design = designCertifiedBilinearController(model, sector, spec);
    check(design.certified, "design should certify: " + design.message + " / " + sdpStatusName(design.solverStatus));
    check(verifyCertifiedDesign(model, sector, spec, design.gain, design.X, design.mu).certified, "verifier accepts the design");

    // Certificate rejection when a solution is perturbed.
    auto badGain = design.gain; badGain[0] += 5.0;
    check(!verifyCertifiedDesign(model, sector, spec, badGain, design.X, design.mu).certified, "perturbed gain is rejected");
    auto badX = design.X; badX[0] *= 40.0;
    check(!verifyCertifiedDesign(model, sector, spec, design.gain, badX, design.mu).certified, "perturbed X is rejected");
    auto asymmetric = design.X; asymmetric[1] += 1e-3;
    check(!verifyCertifiedDesign(model, sector, spec, design.gain, asymmetric, design.mu).certified, "asymmetric X is rejected");
    check(!verifyCertifiedDesign(model, sector, spec, design.gain, design.X, -design.mu).certified, "negative mu is rejected");

    // Closed-loop invariance from the boundary of E, with adversarial sector residuals.
    const double det = design.X[0] * design.X[3] - design.X[1] * design.X[2];
    const std::array<double, 4> P{design.X[3] / det, -design.X[1] / det, -design.X[2] / det, design.X[0] / det};
    const auto V = [&](double a, double b) { return P[0] * a * a + (P[1] + P[2]) * a * b + P[3] * b * b; };
    // The Cholesky factor of X maps the unit circle onto the boundary of E.
    const double l11 = std::sqrt(design.X[0]), l21 = design.X[2] / l11, l22 = std::sqrt(design.X[3] - l21 * l21);
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> angle(0.0, 6.283185307179586);
    std::size_t steps = 0;
    bool invariant = true, decays = true, inBox = true;
    for (int trial = 0; trial < 400; ++trial) {
        const double t = angle(rng);
        std::vector<double> z{l11 * std::cos(t), l21 * std::cos(t) + l22 * std::sin(t)};
        near(V(z[0], z[1]), 1.0, 1e-9, "start on the boundary");
        for (int k = 0; k < 400; ++k, ++steps) {
            const std::vector<double> u{design.gain[0] * z[0] + design.gain[1] * z[1]};
            inBox = inBox && std::abs(u[0]) <= spec.inputBounds[0] * (1 + 1e-9);
            auto next = stepBilinearKrom(model, z, u);
            const double size = sector.cState * std::hypot(z[0], z[1]) + sector.cInput * std::abs(u[0]);
            const double direction = trial % 2 ? angle(rng) : std::atan2(P[2] * next[0] + P[3] * next[1], P[0] * next[0] + P[1] * next[1]);
            next[0] += size * std::cos(direction);
            next[1] += size * std::sin(direction);
            const double before = V(z[0], z[1]), after = V(next[0], next[1]);
            decays = decays && after <= spec.contraction * before * (1 + 1e-9) + 1e-14;
            invariant = invariant && after <= 1.0 + 1e-9;
            z = next;
        }
        decays = decays && V(z[0], z[1]) <= std::pow(spec.contraction, 400) * (1 + 1e-6);
    }
    check(inBox, "input stays in the box");
    check(decays, "certified decay holds along closed-loop trajectories");
    check(invariant, "certified set is invariant");
    check(steps == 160000, "simulation ran");

    // An impossible request is an explicit infeasible result with no controller.
    const auto impossible = designCertifiedBilinearController(model, sectorOf(0.6, 0.2), spec);
    check(!impossible.certified && impossible.gain.empty() && impossible.X.empty(), "no controller without a certificate");
    check(impossible.status == CertifiedDesignStatus::infeasible, std::string("large sector: ") + certifiedDesignStatusName(impossible.status));
}

// Exercises state scaling, two inputs (four vertices) and the one-sided sector branches.
void scaledMultiInput() {
    BilinearKrom model;
    model.stateDimension = 2; model.inputDimension = 2;
    model.A = {1.01, 0.08, -0.02, 0.97};
    model.B0 = {0.5, 0.0, 0.1, 0.4};
    model.B = {{0.01, 0.0, 0.0, -0.01}, {0.0, 0.015, 0.01, 0.0}};
    const std::array<std::pair<double, double>, 4> sectors{{{0.003, 0.01}, {0.004, 0.0}, {0.0, 0.01}, {0.0, 0.0}}};
    for (const auto& [cx, cu] : sectors) {
        const double sx = cx, su = cu;
        const auto sector = sectorOf(sx, su);
        CertifiedDesignSpec spec;
        spec.inputBounds = {0.6, 0.8};
        spec.stateRadius = 5.0;
        spec.contraction = 0.99;
        const auto design = designCertifiedBilinearController(model, sector, spec);
        check(design.certified, "scaled design: " + design.message);
        const double det = design.X[0] * design.X[3] - design.X[1] * design.X[1];
        const auto V = [&](const std::vector<double>& z) {
            return (design.X[3] * z[0] * z[0] - 2 * design.X[1] * z[0] * z[1] + design.X[0] * z[1] * z[1]) / det;
        };
        check(design.X[0] <= 25.0 + 1e-6 && design.X[3] <= 25.0 + 1e-6, "E lies in the trusted ball");
        const double l11 = std::sqrt(design.X[0]), l21 = design.X[2] / l11, l22 = std::sqrt(design.X[3] - l21 * l21);
        std::mt19937 rng(23);
        std::uniform_real_distribution<double> angle(0.0, 6.283185307179586);
        bool holds = true;
        for (int trial = 0; trial < 200; ++trial) {
            const double t = angle(rng);
            std::vector<double> z{l11 * std::cos(t), l21 * std::cos(t) + l22 * std::sin(t)};
            for (int k = 0; k < 200; ++k) {
                const std::vector<double> u{design.gain[0] * z[0] + design.gain[1] * z[1], design.gain[2] * z[0] + design.gain[3] * z[1]};
                holds = holds && std::abs(u[0]) <= 0.6 * (1 + 1e-9) && std::abs(u[1]) <= 0.8 * (1 + 1e-9);
                auto next = stepBilinearKrom(model, z, u);
                const double size = sx * std::hypot(z[0], z[1]) + su * std::hypot(u[0], u[1]), d = angle(rng);
                next[0] += size * std::cos(d); next[1] += size * std::sin(d);
                holds = holds && V(next) <= spec.contraction * V(z) * (1 + 1e-9) + 1e-13;
                z = next;
            }
        }
        check(holds, "scaled closed loop keeps box, decay and invariance");
    }
}

void queries() {
    using dynsyslab::math::evaluateControlCertificateQuery;
    check(!evaluateControlCertificateQuery("not_ours", json::Value::Object{}), "foreign operations fall through");
    const auto sdp = evaluateControlCertificateQuery("sdp_feasibility", json::parse(
        R"({"variables":1,"objective":[1],"margin":0,"blocks":[{"dimension":2,"constant":[[0.1,1],[1,0.1]],"coefficients":[[[1,0],[0,1]]]}]})"));
    check(sdp && sdp->at("status").string() == "converged", "sdp query converges");
    near(sdp->at("objective").number(), 0.9, 1e-5, "sdp query value");

    const auto truth = plant();
    const auto d = sampled(truth, 40, 5);
    json::Value::Object request{{"state_dimension", 2.0}, {"input_dimension", 1.0},
        {"z", denseRowsToJson(d.z, d.samples, 2)}, {"z_next", denseRowsToJson(d.zNext, d.samples, 2)},
        {"u", denseRowsToJson(d.u, d.samples, 1)}};
    const auto fit = evaluateControlCertificateQuery("bilinear_fit", request);
    check(fit && fit->at("identifiable").boolean(), "fit query");
    auto sectorRequest = request;
    sectorRequest["model"] = fit->at("model");
    const auto sector = evaluateControlCertificateQuery("residual_sector", sectorRequest);
    check(sector && sector->at("feasible").boolean() && sector->at("c_x").number() < 1e-9, "sector query");
    const auto design = evaluateControlCertificateQuery("certified_design", json::Value::Object{
        {"model", fit->at("model")}, {"sector", json::Value::Object{{"c_x", 0.004}, {"c_u", 0.002}}},
        {"input_bounds", json::Value::Array{1.0}}, {"state_radius", 1.0}, {"contraction", 0.995}});
    check(design && design->at("certified").boolean() && design->at("certificate").at("certified").boolean(), "design query certifies");
    // A joint-LP sector is marked unsupported and cannot produce a gain.
    check(!sector->at("supported").boolean(true), "joint LP sector is marked unsupported");
    const auto refused = evaluateControlCertificateQuery("certified_design", json::Value::Object{
        {"model", fit->at("model")}, {"sector", *sector},
        {"input_bounds", json::Value::Array{1.0}}, {"state_radius", 1.0}, {"contraction", 0.995}});
    check(refused && refused->at("status").string() == "sector_unsupported" && !refused->at("certified").boolean(true)
        && refused->at("K").isNull(), "design refuses an unsupported sector");
    auto twoStage = sectorRequest;
    twoStage["method"] = "zero_input_first";
    const auto unsupported = evaluateControlCertificateQuery("residual_sector", twoStage);
    check(unsupported && !unsupported->at("supported").boolean(true) && !unsupported->at("feasible").boolean(true),
        "two-stage sector needs zero-input holds");
    rejects([] { return evaluateControlCertificateQuery("certified_design", json::Value::Object{}); }, "malformed design query is rejected");
}

}  // namespace

int main() {
    const std::pair<const char*, void (*)()> suites[] = {{"lyapunov", lyapunovFeasibility}, {"known_values", knownValues},
        {"bilinear", bilinearRecovery}, {"sector", sectorEstimation}, {"certified", certifiedStabilization},
        {"scaled_multi_input", scaledMultiInput}, {"queries", queries}};
    int failures = 0;
    for (const auto& [name, run] : suites) {
        try { run(); std::cout << "PASS " << name << "\n"; }
        catch (const std::exception& e) { ++failures; std::cout << "FAIL " << name << ": " << e.what() << "\n"; }
    }
    std::cout << checks << " checks, " << failures << " failed suites\n";
    return failures ? 1 : 0;
}
