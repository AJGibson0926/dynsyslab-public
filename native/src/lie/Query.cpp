#include "Exact.hpp"
#include <cmath>

namespace dynsyslab::core::lie {
namespace {
using J = json::Value;
std::size_t count(const J &v, std::size_t fallback, std::size_t maximum) {
  if (v.isNull())
    return fallback;
  if (!std::holds_alternative<double>(v.storage()))
    throw std::invalid_argument("Count must be a JSON integer");
  auto n = v.number();
  if (!std::isfinite(n) || n < 0 || n > double(maximum) || std::floor(n) != n)
    throw std::invalid_argument("Count out of range");
  return std::size_t(n);
}
const J::Array &list(const J &v) {
  if (!std::holds_alternative<J::Array>(v.storage()))
    throw std::invalid_argument("Expected JSON array");
  return v.array();
}
std::vector<std::string> strings(const J &v) {
  std::vector<std::string> out;
  for (const auto &a : list(v)) {
    if (!std::holds_alternative<std::string>(a.storage()))
      throw std::invalid_argument("Exact coefficients/points must be strings");
    detail::rational(a.string());
    out.push_back(a.string());
  }
  return out;
}
Polynomial polynomial(const J &v, std::size_t n) {
  Polynomial p{n, {}};
  for (const auto &term : list(v)) {
    Term t;
    for (const auto &x : list(term.at("powers")))
      t.powers.push_back(unsigned(count(x, 0, 128)));
    if (!std::holds_alternative<std::string>(term.at("coefficient").storage()))
      throw std::invalid_argument(
          "Exact coefficient must be an integer or p/q string");
    t.coefficient = term.at("coefficient").string();
    p.terms.push_back(t);
  }
  return detail::write(detail::read(p), n);
}
VectorField field(const J &v, std::size_t n) {
  VectorField f;
  for (const auto &p : list(v))
    f.push_back(polynomial(p, n));
  if (f.size() != n)
    throw std::invalid_argument("Vector field must have dimension components");
  return f;
}
J array(const std::vector<std::string> &a) {
  J::Array out;
  for (const auto &v : a)
    out.emplace_back(v);
  return out;
}
J dimensions(const std::vector<std::size_t> &a) {
  J::Array out;
  for (auto v : a)
    out.emplace_back(double(v));
  return out;
}
J encode(const Polynomial &p) {
  J::Array out;
  for (const auto &t : p.terms) {
    J::Array powers;
    for (auto k : t.powers)
      powers.emplace_back(double(k));
    out.emplace_back(
        J::Object{{"coefficient", t.coefficient}, {"powers", powers}});
  }
  return out;
}
J encode(const VectorField &f) {
  J::Array out;
  for (const auto &p : f)
    out.push_back(encode(p));
  return out;
}
J matrices(const std::vector<std::vector<std::string>> &a) {
  J::Array out;
  for (const auto &v : a)
    out.push_back(array(v));
  return out;
}
} // namespace
J analyze(const J &request) {
  if (request.at("schema").string() != "dynsyslab.lie-algebra-query/v1")
    throw std::invalid_argument("Unknown Lie algebra query schema");
  auto n = count(request.at("dimension"), 0, 16);
  if (!n)
    throw std::invalid_argument("Require dimension 1..16");
  auto action = request.at("action").string("algebra");
  const auto &budget = request.at("limits");
  Limits limits;
  limits.maximumBasis = count(budget.at("maximum_basis"), 32, 64);
  limits.maximumOperations =
      count(budget.at("maximum_operations"), 4096, 100000);
  limits.maximumDegree = unsigned(count(budget.at("maximum_degree"), 24, 64));
  detail::limitsValid(limits);
  J::Object out{{"schema", "dynsyslab.lie-algebra-result/v1"},
                {"action", action},
                {"dimension", double(n)},
                {"arithmetic", "exact rational"},
                {"bracket_convention", "[X,Y]=DY*X-DX*Y; [L_X,L_Y]=L_[X,Y]"},
                {"coordinate_names", request.at("coordinate_names")},
                {"provenance", request.at("provenance")}};
  try {
    if (action == "box_bound") {
      std::vector<Polynomial> residuals;
      for (const auto &p : list(request.at("residuals")))
        residuals.push_back(polynomial(p, n));
      auto b = residualBoxBound(residuals, strings(request.at("radii")));
      out["uniform_l1_bound"] = b.uniformL1;
      out["lipschitz_l1_bound"] = b.lipschitzL1;
      out["vanishes_at_origin"] = b.vanishesAtOrigin;
      out["proportional_state_bound_certified"] = b.vanishesAtOrigin;
      out["scope"] = "Known polynomial residual on the closed origin-centered "
                     "box. Induced 1-norm bound; not a learned-model "
                     "confidence bound or stability certificate.";
    } else if (action == "bracket")
      out["value"] =
          encode(bracket(field(request.at("x"), n), field(request.at("y"), n)));
    else if (action == "poisson")
      out["value"] = encode(poissonBracket(polynomial(request.at("a"), n),
                                           polynomial(request.at("b"), n),
                                           strings(request.at("tensor"))));
    else if (action == "flow_bracket")
      out["value"] = encode(flowBracket(field(request.at("field"), n),
                                        polynomial(request.at("a"), n),
                                        polynomial(request.at("b"), n)));
    else if (action == "derivative")
      out["value"] =
          encode(differentiate(polynomial(request.at("observable"), n),
                               count(request.at("axis"), 0, n - 1)));
    else if (action == "lie_derivative")
      out["value"] =
          encode(lieDerivative(field(request.at("field"), n),
                               polynomial(request.at("observable"), n)));
    else if (action == "product")
      out["value"] = encode(multiply(polynomial(request.at("a"), n),
                                     polynomial(request.at("b"), n)));
    else if (action == "algebra") {
      std::vector<VectorField> generators;
      for (const auto &f : list(request.at("generators")))
        generators.push_back(field(f, n));
      const auto result = generatedAlgebra(generators, limits);
      J::Array basis;
      for (const auto &f : result.basis)
        basis.push_back(encode(f));
      out["status"] = result.status;
      out["reason"] = result.reason;
      out["basis"] = basis;
      out["basis_words"] = array(result.words);
      out["growth"] = dimensions(result.growth);
      out["basis_dimension"] = double(result.basis.size());
      out["operations"] = double(result.operations);
      if (result.status == "closed") {
        out["structure_constants"] = array(result.structureConstants);
        out["structure_layout"] =
            "i,j,k: [basis_i,basis_j]=sum_k c_ijk basis_k";
        out["derived_dimensions"] = dimensions(result.derivedDimensions);
        out["lower_central_dimensions"] =
            dimensions(result.lowerCentralDimensions);
        out["center_basis_coefficients"] = matrices(result.center);
        out["killing_form"] = array(result.killingForm);
        out["solvable"] = result.solvable;
        out["nilpotent"] = result.nilpotent;
        if (!request.at("ideal_generators").isNull()) {
          std::vector<VectorField> seeds;
          for (const auto &f : list(request.at("ideal_generators")))
            seeds.push_back(field(f, n));
          auto ideal = generatedIdeal(result, seeds, limits);
          J::Array fields;
          for (const auto &f : ideal.basis)
            fields.push_back(encode(f));
          J::Object info{{"status", ideal.status},
                         {"reason", ideal.reason},
                         {"basis", fields},
                         {"dimension", double(ideal.basis.size())}};
          if (ideal.status == "closed") {
            info["nilpotent"] = ideal.nilpotent;
            info["solvable"] = ideal.solvable;
            info["lower_central_dimensions"] =
                dimensions(ideal.lowerCentralDimensions);
          }
          out["generated_ideal"] = info;
        }
      }
      if (!request.at("points").isNull()) {
        J::Array ranks;
        for (const auto &point : list(request.at("points")))
          ranks.emplace_back(double(rankAt(result.basis, strings(point))));
        out["pointwise_distribution_ranks"] = ranks;
        out["rank_scope"] =
            "Exact at supplied rational points; lower bounds if basis is "
            "incomplete. Not a global controllability certificate.";
      }
    } else if (action == "module" || action == "defect") {
      std::vector<Polynomial> seeds;
      for (const auto &p : list(request.at("observables")))
        seeds.push_back(polynomial(p, n));
      if (action == "defect") {
        auto r = representationDefect(seeds, field(request.at("field"), n),
                                      strings(request.at("matrix")));
        out["residuals"] = encode(r);
        bool zero = true;
        for (const auto &p : r)
          zero = zero && p.terms.empty();
        out["identically_zero"] = zero;
      } else {
        std::vector<VectorField> ops;
        for (const auto &f : list(request.at("operators")))
          ops.push_back(field(f, n));
        if (request.at("include_partials").boolean())
          for (std::size_t j = 0; j < n; ++j) {
            VectorField f(n, Polynomial{n, {}});
            f[j].terms.push_back({std::vector<unsigned>(n), "1"});
            ops.push_back(f);
          }
        auto r = invariantModule(seeds, ops, request.at("extend").boolean(true),
                                 limits);
        out["status"] = r.status;
        out["reason"] = r.reason;
        out["basis"] = encode(r.basis);
        out["basis_dimension"] = double(r.basis.size());
        out["matrices"] = matrices(r.matrices);
        out["matrix_layout"] =
            "operator, output observable, predictor observable: L_X Psi=A_X "
            "Psi; coefficient-column representation is A_X^T";
        J::Array residuals;
        for (const auto &a : r.residuals)
          residuals.push_back(encode(a));
        out["residuals"] = residuals;
        out["global_linear_state_readout_certified"] =
            r.linearStateReconstruction;
        out["linear_state_readout"] = array(r.linearStateReadout);
        out["eigenfunction_basis_claimed"] = false;
      }
    } else
      throw std::invalid_argument("Unknown Lie algebra action");
  } catch (const detail::Limit &e) {
    out["status"] = "limit_reached";
    out["reason"] = e.what();
  }
  return out;
}
} // namespace dynsyslab::core::lie
