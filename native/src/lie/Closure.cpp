#include "Exact.hpp"

namespace dynsyslab::core::lie::detail {
// Our second-order cascade certificate: H=[F,q^2 d_p],
// [H,q^n d_p]=-(n+2)q^(n+1)d_p whenever F=p d_q+(A(q)p+C(q))d_p.
// Combined with the faithful representation obstruction of Katayama/Susuki,
// arXiv:2611.20438, Theorem 1. A bounded unsuccessful search alone is NOT a
// proof.
bool powerCascade(const std::vector<F> &generators, const Span &basis) {
  F y(2);
  y[1][{2, 0}] = 1;
  if (!basis.reduce(flatten(y)).second.empty())
    return false;
  for (const auto &f : generators)
    if (f.size() == 2 && f[0] == P{{{0, 1}, Q(1)}}) {
      bool affineVelocity = true;
      for (const auto &[k, v] : f[1])
        if (k[1] > 1)
          affineVelocity = false;
      if (affineVelocity)
        return true;
    }
  return false;
}
} // namespace dynsyslab::core::lie::detail

namespace dynsyslab::core::lie {
using namespace detail;
AlgebraResult generatedAlgebra(const std::vector<VectorField> &generators,
                               const Limits &limits) {
  limitsValid(limits);
  if (generators.empty() || generators.size() > 64)
    throw std::invalid_argument("Require 1..64 generators");
  const auto n = generators[0].size();
  std::vector<F> g;
  Span span;
  AlgebraResult out;
  for (const auto &f : generators) {
    if (f.size() != n)
      throw std::invalid_argument("Generator dimension mismatch");
    g.push_back(read(f));
  }
  auto insert = [&](const F &f, const std::string &word) {
    auto flat = flatten(f);
    if (span.reduce(flat).second.empty())
      return;
    if (span.rows.size() >= limits.maximumBasis)
      throw Limit("Independent basis budget reached");
    for (const auto &p : f)
      if (degree(p) > limits.maximumDegree)
        throw Limit("Polynomial degree budget reached");
    span.insert(flat);
    out.basis.push_back(write(f));
    out.words.push_back(word);
    out.growth.push_back(span.rows.size());
  };
  out.status = "closed";
  try {
    for (std::size_t i = 0; i < g.size(); ++i)
      insert(g[i], "g" + std::to_string(i));
    for (std::size_t cursor = 0; cursor < span.rows.size(); ++cursor) {
      if (n == 3 && sinusoidalMechanicalCascade(g, span)) {
        out.status = "proven_infinite";
        out.reason = "Sinusoidal mechanical harmonic cascade on s^2+c^2=1: "
                     "Y_n=sin(n*q)*d_v; "
                     "Y_(n+1)=(2*[Y_1,[F,Y_n]]+(n-1)*Y_(n-1))/(n+1). "
                     "Y_0=0; nonzero sine/cosine restoring term and d_v membership "
                     "verified exactly. Independent on every open physical patch. "
                     "Finite smooth spaces invariant under F and d_v contain only "
                     "constants on each connected physical domain.";
        return out;
      }
      if (n == 3 && inverseRadiusCascade(g, span)) {
        out.status = "proven_infinite";
        out.reason = "Inverse-radius cascade in coordinates (r,v,w): G=w*d_v, "
                     "H=[F,G], [H,w^n*d_v]=(n-2)*w^(n+2)*d_v for odd n>=1. "
                     "Infinitely many independent fields, including on the "
                     "invariant physical surface r*w=1, r>0.";
        return out;
      }
      if (n == 2 && powerCascade(g, span)) {
        out.status = "proven_infinite";
        out.reason = "Second-order power cascade: "
                     "[H,x1^n*d_x2]=-(n+2)*x1^(n+1)*d_x2, n>=2; infinitely "
                     "many independent fields on every open state domain.";
        return out;
      }
      const auto v = unflatten(span.originals[cursor], n);
      const auto word = out.words[cursor];
      for (std::size_t j = 0; j < g.size(); ++j) {
        if (out.operations++ >= limits.maximumOperations)
          throw Limit("Bracket operation budget reached");
        insert(commutator(g[j], v),
               "[g" + std::to_string(j) + "," + word + "]");
      }
    }
    // All generator adjoints preserve the completed span. Jacobi implies Lie
    // closure.
    const auto r = span.rows.size();
    out.structureConstants.assign(r * r * r, "0");
    for (std::size_t i = 0; i < r; ++i)
      for (std::size_t j = 0; j < r; ++j) {
        auto [c, residual] = span.reduce(flatten(commutator(
            unflatten(span.originals[i], n), unflatten(span.originals[j], n))));
        if (!residual.empty())
          throw std::logic_error("Internal closure inconsistency");
        for (std::size_t k = 0; k < r; ++k)
          out.structureConstants[(i * r + j) * r + k] = str(c[k]);
      }
    properties(out);
    out.reason =
        "Exact rational polynomial identities certify finite Lie closure; this "
        "alone does not prove an eigenfunction embedding or controllability.";
  } catch (const Limit &e) {
    out.status = "limit_reached";
    out.reason = e.what();
  }
  return out;
}

AlgebraResult generatedIdeal(const AlgebraResult &ambient,
                             const std::vector<VectorField> &seeds,
                             const Limits &limits) {
  limitsValid(limits);
  if (ambient.status != "closed" || ambient.basis.empty() || seeds.empty() ||
      seeds.size() > 64)
    throw std::invalid_argument("Ideal requires a nonempty certified finite "
                                "ambient algebra and 1..64 seeds");
  const auto n = ambient.basis[0].size();
  Span space, ideal;
  std::vector<F> fields;
  for (const auto &f : ambient.basis) {
    fields.push_back(read(f));
    space.insert(flatten(fields.back()));
  }
  for (const auto &f : seeds) {
    if (f.size() != n)
      throw std::invalid_argument("Ideal seed dimension mismatch");
    auto v = flatten(read(f));
    if (!space.reduce(v).second.empty())
      throw std::invalid_argument("Ideal seed is outside ambient algebra");
    ideal.insert(v);
  }
  AlgebraResult out;
  std::size_t operations = 0;
  try {
    if (ideal.rows.size() > limits.maximumBasis)
      throw Limit("Ideal basis budget reached");
    for (std::size_t i = 0; i < ideal.rows.size(); ++i) {
      auto f = unflatten(ideal.originals[i], n);
      for (const auto &g : fields) {
        if (operations++ >= limits.maximumOperations)
          throw Limit("Ideal operation budget reached");
        auto v = flatten(commutator(g, f));
        if (!ideal.reduce(v).second.empty()) {
          if (ideal.rows.size() >= limits.maximumBasis)
            throw Limit("Ideal basis budget reached");
          ideal.insert(v);
        }
      }
    }
    std::vector<VectorField> result;
    for (const auto &v : ideal.originals)
      result.push_back(write(unflatten(v, n)));
    if (result.empty())
      result.push_back(write(F(n)));
    out = generatedAlgebra(result, limits);
    if (out.status == "closed")
      out.reason = "Smallest ideal containing the supplied seeds; closure "
                   "under all ambient adjoints checked exactly.";
  } catch (const Limit &e) {
    out.status = "limit_reached";
    out.reason = e.what();
    for (const auto &v : ideal.originals)
      out.basis.push_back(write(unflatten(v, n)));
  }
  return out;
}

ModuleResult invariantModule(const std::vector<Polynomial> &seeds,
                             const std::vector<VectorField> &operators,
                             bool extend, const Limits &limits) {
  limitsValid(limits);
  if (seeds.empty() || seeds.size() > 64 || operators.empty() ||
      operators.size() > 64)
    throw std::invalid_argument("Require 1..64 observables and operators");
  auto n = seeds[0].dimension;
  Span span;
  std::vector<F> ops;
  ModuleResult out;
  for (const auto &f : operators) {
    if (f.size() != n)
      throw std::invalid_argument("Module operator dimension mismatch");
    ops.push_back(read(f));
  }
  for (const auto &p : seeds) {
    if (p.dimension != n)
      throw std::invalid_argument("Observable dimension mismatch");
    if (span.rows.size() >= limits.maximumBasis)
      throw std::invalid_argument("Seed count exceeds basis budget");
    if (!span.insert(read(p)))
      throw std::invalid_argument(
          "Observable seeds must be linearly independent");
  }
  out.status = "closed";
  std::size_t operations = 0;
  try {
    if (extend)
      for (std::size_t cursor = 0; cursor < span.rows.size(); ++cursor) {
        auto p = span.originals[cursor];
        for (const auto &f : ops) {
          if (operations++ >= limits.maximumOperations)
            throw Limit("Module operation budget reached");
          auto next = action(f, p);
          if (span.reduce(next).second.empty())
            continue;
          if (span.rows.size() >= limits.maximumBasis)
            throw Limit("Module dimension budget reached");
          if (degree(next) > limits.maximumDegree)
            throw Limit("Module degree budget reached");
          span.insert(next);
        }
      }
  } catch (const Limit &e) {
    out.status = "limit_reached";
    out.reason = e.what();
  }
  for (const auto &p : span.originals)
    out.basis.push_back(write(p, n));
  for (const auto &f : ops) {
    std::vector<std::string> matrix;
    std::vector<Polynomial> residuals;
    for (const auto &p : span.originals) {
      auto [c, r] = span.reduce(action(f, p));
      for (const auto &a : c)
        matrix.push_back(str(a));
      residuals.push_back(write(r, n));
      if (!r.empty() && out.status == "closed")
        out.status = "not_invariant";
    }
    out.matrices.push_back(std::move(matrix));
    out.residuals.push_back(std::move(residuals));
  }
  out.linearStateReconstruction = true;
  for (std::size_t j = 0; j < n; ++j) {
    Key k(n);
    k[j] = 1;
    auto [c, r] = span.reduce(P{{k, Q(1)}});
    if (!r.empty())
      out.linearStateReconstruction = false;
    for (const auto &a : c)
      out.linearStateReadout.push_back(str(a));
  }
  if (!out.linearStateReconstruction)
    out.linearStateReadout.clear();
  if (out.reason.empty())
    out.reason =
        out.status == "closed"
            ? "Exact common-invariant function space. Eigenfunction basis "
              "requires diagonalizing the drift action; not assumed."
            : "Nonzero exact residuals refute invariance of the supplied span, "
              "not existence of a larger finite span.";
  return out;
}

std::vector<Polynomial>
representationDefect(const std::vector<Polynomial> &basis, const VectorField &f,
                     const std::vector<std::string> &matrix) {
  const auto r = basis.size();
  if (!r || matrix.size() != r * r)
    throw std::invalid_argument("Defect matrix shape mismatch");
  auto field = read(f);
  std::vector<P> p;
  for (const auto &q : basis) {
    if (q.dimension != f.size())
      throw std::invalid_argument("Defect dimension mismatch");
    p.push_back(read(q));
  }
  std::vector<Polynomial> out;
  for (std::size_t i = 0; i < r; ++i) {
    auto v = action(field, p[i]);
    for (std::size_t j = 0; j < r; ++j)
      axpy(v, p[j], -rational(matrix[i * r + j]));
    out.push_back(write(v, f.size()));
  }
  return out;
}
} // namespace dynsyslab::core::lie
