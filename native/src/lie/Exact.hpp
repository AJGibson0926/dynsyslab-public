#pragma once
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <dynsyslab/core/LieAlgebra.hpp>
#include <map>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core::lie::detail {
using Q = boost::multiprecision::cpp_rational;
using Key = std::vector<unsigned>;
using P = std::map<Key, Q>;
using F = std::vector<P>;
struct Limit : std::runtime_error {
  using std::runtime_error::runtime_error;
};
Q rational(const std::string &);
std::string str(const Q &);
void add(P &, const Key &, const Q &);
void axpy(P &, const P &, const Q &);
P read(const Polynomial &);
Polynomial write(const P &, std::size_t dimension);
F read(const VectorField &);
VectorField write(const F &);
P derivative(const P &, std::size_t);
P product(const P &, const P &);
P action(const F &, const P &);
F commutator(const F &, const F &);
P flatten(const F &);
F unflatten(const P &, std::size_t);
unsigned degree(const P &);
Q evaluate(const P &, const std::vector<Q> &);
void limitsValid(const Limits &);

// Sparse exact echelon basis, retaining coordinates in the ORIGINAL input
// basis. Gaussian elimination over Q; no sampled rank or tolerance enters a
// certificate.
struct Span {
  std::vector<P> originals, rows;
  std::vector<std::vector<Q>> transforms;
  std::pair<std::vector<Q>, P> reduce(P v) const {
    std::vector<Q> c(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
      const auto at = v.find(rows[i].begin()->first);
      if (at == v.end())
        continue;
      Q a = at->second;
      axpy(v, rows[i], -a);
      for (std::size_t j = 0; j < transforms[i].size(); ++j)
        c[j] += a * transforms[i][j];
    }
    return {c, v};
  }
  bool insert(const P &v) {
    auto [c, r] = reduce(v);
    if (r.empty())
      return false;
    Q pivot = r.begin()->second;
    for (auto &entry : r)
      entry.second /= pivot;
    for (auto &a : c)
      a = -a / pivot;
    c.push_back(1 / pivot);
    originals.push_back(v);
    rows.push_back(std::move(r));
    transforms.push_back(std::move(c));
    return true;
  }
};
std::vector<std::vector<Q>> nullspace(std::vector<std::vector<Q>> rows,
                                      std::size_t n);
void properties(AlgebraResult &);
bool powerCascade(const std::vector<F> &generators, const Span &basis);
bool inverseRadiusCascade(const std::vector<F> &generators, const Span &basis);
bool sinusoidalMechanicalCascade(const std::vector<F>& generators,const Span& basis);
} // namespace dynsyslab::core::lie::detail
