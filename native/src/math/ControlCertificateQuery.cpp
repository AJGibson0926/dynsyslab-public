// Stateless JSON adapter over the certified bilinear control tools. Algorithms and
// citations live in core/DenseSdp, BilinearKrom, ResidualSector and CertifiedDesign.
#include <dynsyslab/math/ControlCertificateQuery.hpp>

#include <dynsyslab/core/BilinearKrom.hpp>
#include <dynsyslab/core/CertifiedDesign.hpp>
#include <dynsyslab/core/DenseSdp.hpp>
#include <dynsyslab/core/ResidualSector.hpp>
#include <dynsyslab/core/SupportedResidualSector.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::math {
namespace {

namespace c = dynsyslab::core;
using J = core::json::Value;
using O = J::Object;
using A = J::Array;

double number(const J& v, const char* name) {
    if (!std::holds_alternative<double>(v.storage()) || !std::isfinite(v.number()))
        throw std::invalid_argument(std::string(name) + " must be a finite number");
    return v.number();
}
double numberOr(const J& q, const char* key, double fallback) { return q.at(key).isNull() ? fallback : number(q.at(key), key); }
std::size_t count(const J& v, const char* name) {
    const double x = number(v, name);
    if (x < 0.0 || x != std::floor(x) || x > 1e9) throw std::invalid_argument(std::string(name) + " must be a non-negative integer");
    return static_cast<std::size_t>(x);
}
std::vector<double> reals(const J& v, const char* name) {
    std::vector<double> out;
    for (const auto& x : v.array()) out.push_back(number(x, name));
    return out;
}
J vectorJson(const std::vector<double>& v) { A a; a.reserve(v.size()); for (double x : v) a.emplace_back(x); return a; }

std::vector<double> matrix(const J& v, std::size_t expectedRows, std::size_t expectedColumns, const char* name) {
    std::size_t rows = 0, columns = 0;
    auto out = c::denseRowsFromJson(v, rows, columns);
    if (rows != expectedRows || columns != expectedColumns) throw std::invalid_argument(std::string(name) + " has the wrong shape");
    return out;
}

c::SdpOptions solverOptions(const J& q) {
    c::SdpOptions o;
    o.tolerance = numberOr(q, "tolerance", o.tolerance);
    o.infeasibilityTolerance = numberOr(q, "infeasibility_tolerance", o.infeasibilityTolerance);
    if (!q.at("max_iterations").isNull()) o.maxIterations = count(q.at("max_iterations"), "max_iterations");
    return o;
}

J certificateJson(const c::LmiCertificate& x) {
    return O{{"certified", x.certified}, {"margin", x.margin}, {"worst_min_eigenvalue", x.worstMinEigenvalue},
        {"rounding_allowance", x.roundingAllowance}, {"block_min_eigenvalues", vectorJson(x.blockMinEigenvalues)},
        {"reason", x.reason}, {"scope", x.scope}};
}

c::SdpProblem sdpProblem(const J& q) {
    c::SdpProblem p;
    p.variables = count(q.at("variables"), "variables");
    if (!q.at("objective").isNull()) p.objective = reals(q.at("objective"), "objective");
    for (const auto& item : q.at("blocks").array()) {
        c::LmiBlock block;
        block.dimension = count(item.at("dimension"), "dimension");
        block.constant = matrix(item.at("constant"), block.dimension, block.dimension, "constant");
        const auto& terms = item.at("coefficients").array();
        if (terms.size() != p.variables) throw std::invalid_argument("Each block needs one coefficient per variable");
        for (const auto& term : terms)
            block.coefficients.push_back(term.isNull() ? std::vector<double>{}
                : matrix(term, block.dimension, block.dimension, "coefficient"));
        p.blocks.push_back(std::move(block));
    }
    return p;
}

J sdpFeasibility(const J& q) {
    const auto problem = sdpProblem(q);
    const double margin = numberOr(q, "margin", 0.0);
    const auto result = c::solveDenseSdp(problem, solverOptions(q));
    O out{{"status", c::sdpStatusName(result.status)}, {"iterations", double(result.iterations)},
        {"primal_residual", result.primalResidual}, {"dual_residual", result.dualResidual}, {"gap", result.gap}};
    if (result.status == c::SdpStatus::converged) {
        out["x"] = vectorJson(result.x);
        out["objective"] = result.objective;
        out["certificate"] = certificateJson(c::verifyLmiCertificate(problem, result.x, margin));
    } else {
        out["certificate"] = O{{"certified", false}, {"reason", "solver did not return a converged solution"}};
    }
    return out;
}

c::BilinearDataset dataset(const J& q, std::size_t n, std::size_t m) {
    c::BilinearDataset d;
    std::size_t rows = 0, columns = 0;
    d.z = c::denseRowsFromJson(q.at("z"), rows, columns);
    if (columns != n) throw std::invalid_argument("z has the wrong width");
    d.samples = rows;
    d.zNext = matrix(q.at("z_next"), rows, n, "z_next");
    d.u = matrix(q.at("u"), rows, m, "u");
    return d;
}

J residualStats(const c::BilinearResidualStats& s) {
    return O{{"root_mean_square", s.rootMeanSquare}, {"maximum", s.maximum}, {"relative", s.relative}};
}

J bilinearFit(const J& q) {
    const auto n = count(q.at("state_dimension"), "state_dimension");
    const auto m = count(q.at("input_dimension"), "input_dimension");
    c::BilinearFitOptions options;
    options.ridge = numberOr(q, "ridge", 0.0);
    options.sampleTime = numberOr(q, "sample_time", 0.0);
    options.provenance = q.at("provenance").string();
    options.allowRankDeficient = q.at("allow_rank_deficient").boolean(false);
    // Optional frozen_A (n x n rows): regress only B0 and B_i on z_next - A z.
    const bool frozen = !q.at("frozen_A").isNull();
    const auto fit = frozen
        ? c::fitBilinearKromFrozenA(dataset(q, n, m), matrix(q.at("frozen_A"), n, n, "frozen_A"), n, m, options)
        : c::fitBilinearKrom(dataset(q, n, m), n, m, options);
    return O{{"model", c::bilinearKromToJson(fit.model)}, {"regressors", double(fit.regressors)}, {"rank", double(fit.rank)},
        {"identifiable", fit.identifiable}, {"frozen_A", frozen},
        {"condition_number", std::isfinite(fit.conditionNumber) ? J(fit.conditionNumber) : J()},
        {"residual", residualStats(fit.stats)}};
}

c::ResidualSectorOptions sectorOptions(const J& q) {
    c::ResidualSectorOptions o;
    o.inflation = numberOr(q, "inflation", 1.0);
    o.stateWeight = numberOr(q, "state_weight", 0.0);
    o.inputWeight = numberOr(q, "input_weight", 0.0);
    if (!q.at("sampling_assumption").isNull()) {
        const auto& s = q.at("sampling_assumption");
        o.scenario = c::ScenarioAssumption{s.at("declaration").string(), numberOr(s, "confidence", 0.999)};
    }
    return o;
}

J sectorJson(const c::ResidualSectorBound& b) {
    O out{{"feasible", b.feasible}, {"c_x", b.cState}, {"c_u", b.cInput}, {"inflation", b.inflation},
        {"samples", double(b.samples)}, {"maximum_residual", b.maximumResidual},
        {"guarantee", b.guarantee}, {"message", b.message}};
    if (b.violationProbability) out["violation_probability"] = *b.violationProbability;
    if (b.confidence) out["confidence"] = *b.confidence;
    // The joint LP may put the whole residual on one channel, so it cannot feed a certificate.
    out["supported"] = false;
    out["identification"] = "joint LP envelope; use method zero_input_first for a certified design";
    return out;
}

c::SupportedSectorOptions supportedOptions(const J& q) {
    c::SupportedSectorOptions o;
    o.inflation = numberOr(q, "inflation", o.inflation);
    o.zeroInputTolerance = numberOr(q, "zero_input_tolerance", o.zeroInputTolerance);
    o.zeroResidualTolerance = numberOr(q, "zero_residual_tolerance", o.zeroResidualTolerance);
    o.stateRadius = numberOr(q, "state_radius", 0.0);
    o.inputBound = numberOr(q, "input_bound", 0.0);
    o.minimumCoverage = numberOr(q, "minimum_coverage", o.minimumCoverage);
    if (!q.at("minimum_zero_input_pairs").isNull()) o.minimumZeroInputPairs = count(q.at("minimum_zero_input_pairs"), "minimum_zero_input_pairs");
    if (!q.at("minimum_forced_pairs").isNull()) o.minimumForcedPairs = count(q.at("minimum_forced_pairs"), "minimum_forced_pairs");
    if (!q.at("one_sided_sector").isNull())
        o.oneSided = c::OneSidedSectorDeclaration{q.at("one_sided_sector").at("channel").string(), q.at("one_sided_sector").at("reason").string()};
    return o;
}

J supportedSectorJson(const c::SupportedResidualSector& s) {
    return O{{"feasible", s.bound.feasible}, {"supported", s.supported}, {"c_x", s.bound.cState}, {"c_u", s.bound.cInput},
        {"inflation", s.bound.inflation}, {"samples", double(s.bound.samples)}, {"maximum_residual", s.bound.maximumResidual},
        {"method", "zero_input_first"}, {"reason", s.reason}, {"one_sided", s.oneSided},
        {"declaration_reason", s.declarationReason}, {"zero_input_pairs", double(s.zeroInputPairs)},
        {"forced_pairs", double(s.forcedPairs)}, {"largest_zero_input_state_norm", s.largestZeroInputStateNorm},
        {"largest_state_norm", s.largestStateNorm}, {"largest_input_norm", s.largestInputNorm},
        {"guarantee", s.bound.guarantee}};
}

J residualSector(const J& q) {
    if (q.at("method").string() == "zero_input_first") {
        const auto model = c::bilinearKromFromJson(q.at("model"));
        return supportedSectorJson(c::estimateSupportedResidualSector(model,
            dataset(q, model.stateDimension, model.inputDimension), supportedOptions(q)));
    }
    const auto options = sectorOptions(q);
    if (!q.at("model").isNull()) {
        const auto model = c::bilinearKromFromJson(q.at("model"));
        return sectorJson(c::estimateResidualSector(model, dataset(q, model.stateDimension, model.inputDimension), options));
    }
    const auto n = count(q.at("state_dimension"), "state_dimension");
    const auto m = count(q.at("input_dimension"), "input_dimension");
    std::size_t rows = 0, columns = 0;
    const auto z = c::denseRowsFromJson(q.at("z"), rows, columns);
    if (columns != n) throw std::invalid_argument("z has the wrong width");
    const auto u = matrix(q.at("u"), rows, m, "u");
    const auto r = matrix(q.at("residual"), rows, n, "residual");
    return sectorJson(c::estimateResidualSector(z, u, r, rows, n, m, options));
}

J certifiedDesign(const J& q) {
    // A sector estimated from data and marked unsupported must not produce a gain.
    const auto& supported = q.at("sector").at("supported");
    if (!supported.isNull() && !supported.boolean(false))
        return O{{"status", "sector_unsupported"}, {"certified", false},
            {"message", "This data-estimated sector does not identify both residual channels; use method zero_input_first."}};
    const auto model = c::bilinearKromFromJson(q.at("model"));
    c::ResidualSectorBound sector;
    sector.feasible = true;
    sector.cState = number(q.at("sector").at("c_x"), "c_x");
    sector.cInput = number(q.at("sector").at("c_u"), "c_u");
    c::CertifiedDesignSpec spec;
    spec.inputBounds = reals(q.at("input_bounds"), "input_bounds");
    spec.stateRadius = number(q.at("state_radius"), "state_radius");
    spec.contraction = numberOr(q, "contraction", spec.contraction);
    spec.solverMargin = numberOr(q, "solver_margin", spec.solverMargin);
    spec.certificateMargin = numberOr(q, "certificate_margin", spec.certificateMargin);
    spec.solver = solverOptions(q);
    const auto r = c::designCertifiedBilinearController(model, sector, spec);
    O out{{"status", c::certifiedDesignStatusName(r.status)}, {"certified", r.certified}, {"message", r.message},
        {"solver_status", c::sdpStatusName(r.solverStatus)}, {"iterations", double(r.iterations)},
        {"primal_residual", r.primalResidual}, {"dual_residual", r.dualResidual}, {"gap", r.gap},
        {"control_law", "u = K z"}};
    if (r.certified) {
        out["K"] = c::denseRowsToJson(r.gain, model.inputDimension, model.stateDimension);
        out["X"] = c::denseRowsToJson(r.X, model.stateDimension, model.stateDimension);
        out["mu"] = r.mu;
        out["inscribed_radius"] = r.inscribedRadius;
        out["certificate"] = certificateJson(r.certificate);
    }
    return out;
}

}  // namespace

std::optional<core::json::Value> evaluateControlCertificateQuery(const std::string& op, const core::json::Value& q) {
    if (op == "sdp_feasibility") return sdpFeasibility(q);
    if (op == "bilinear_fit") return bilinearFit(q);
    if (op == "residual_sector") return residualSector(q);
    if (op == "certified_design") return certifiedDesign(q);
    return std::nullopt;
}

}  // namespace dynsyslab::math
