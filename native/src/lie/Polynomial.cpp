#include "Exact.hpp"
#include <cctype>
#include <cmath>

namespace dynsyslab::core::lie::detail {
Q rational(const std::string &s) {
  if (s.empty() || s.size() > 512)
    throw std::invalid_argument("Expected bounded integer or p/q string");
  auto integer = [](const std::string &t) {
    std::size_t i = (!t.empty() && (t[0] == '-' || t[0] == '+')) ? 1 : 0;
    if (i == t.size())
      throw std::invalid_argument("Empty integer");
    boost::multiprecision::cpp_int n = 0;
    for (; i < t.size(); ++i) {
      if (!std::isdigit(static_cast<unsigned char>(t[i])))
        throw std::invalid_argument(
            "Exact coefficient requires integer or p/q");
      n = n * 10 + (t[i] - '0');
    }
    return t[0] == '-' ? -n : n;
  };
  auto slash = s.find('/');
  auto n = integer(s.substr(0, slash));
  if (slash == std::string::npos)
    return Q(n);
  auto d = integer(s.substr(slash + 1));
  if (d == 0)
    throw std::invalid_argument("Zero rational denominator");
  return Q(n) / d;
}
std::string str(const Q &q) { return q.str(); }
void add(P &p, const Key &k, const Q &q) {
  if (q == 0)
    return;
  auto &value = p[k];
  value += q;
  if (value != 0) {
    if (boost::multiprecision::msb(abs(numerator(value))) > 16384 ||
        boost::multiprecision::msb(denominator(value)) > 16384)
      throw Limit("Rational coefficient bit budget exceeded");
  } else
    p.erase(k);
  if (p.size() > 8192)
    throw Limit("Polynomial term budget exceeded");
}
void axpy(P &p, const P &a, const Q &q) {
  for (const auto &[k, v] : a)
    add(p, k, q * v);
}
P read(const Polynomial &p) {
  if (p.dimension == 0 || p.dimension > 16 || p.terms.size() > 8192)
    throw std::invalid_argument("Polynomial dimension/term count out of range");
  P out;
  for (const auto &t : p.terms) {
    if (t.powers.size() != p.dimension)
      throw std::invalid_argument("Monomial dimension mismatch");
    if (std::any_of(t.powers.begin(), t.powers.end(),
                    [](unsigned v) { return v > 128; }))
      throw std::invalid_argument("Exponent exceeds 128");
    add(out, t.powers, rational(t.coefficient));
  }
  return out;
}
Polynomial write(const P &p, std::size_t n) {
  Polynomial out{n, {}};
  for (const auto &[k, v] : p)
    out.terms.push_back({k, str(v)});
  return out;
}
F read(const VectorField &f) {
  if (f.empty() || f.size() > 16)
    throw std::invalid_argument("Vector field dimension out of range");
  F out;
  for (const auto &p : f) {
    if (p.dimension != f.size())
      throw std::invalid_argument("Vector field shape mismatch");
    out.push_back(read(p));
  }
  return out;
}
VectorField write(const F &f) {
  VectorField out;
  for (const auto &p : f)
    out.push_back(write(p, f.size()));
  return out;
}
P derivative(const P &p, std::size_t j) {
  P out;
  for (const auto &[k, v] : p) {
    if (j >= k.size())
      throw std::invalid_argument("Derivative axis out of range");
    if (k[j]) {
      auto a = k;
      --a[j];
      add(out, a, v * k[j]);
    }
  }
  return out;
}
P product(const P &a, const P &b) {
  P out;
  if (a.size() * b.size() > 1000000)
    throw Limit("Polynomial product budget exceeded");
  for (const auto &[k, v] : a)
    for (const auto &[l, w] : b) {
      Key n = k;
      for (std::size_t j = 0; j < n.size(); ++j) {
        n[j] += l[j];
        if (n[j] > 128)
          throw Limit("Exponent budget exceeded");
      }
      add(out, n, v * w);
    }
  return out;
}
P action(const F &f, const P &p) {
  P out;
  for (std::size_t j = 0; j < f.size(); ++j)
    axpy(out, product(f[j], derivative(p, j)), 1);
  return out;
}
F commutator(const F &x, const F &y) {
  if (x.size() != y.size())
    throw std::invalid_argument("Bracket dimension mismatch");
  F out(x.size());
  for (std::size_t j = 0; j < x.size(); ++j) {
    out[j] = action(x, y[j]);
    axpy(out[j], action(y, x[j]), -1);
  }
  return out;
}
P flatten(const F &f) {
  P out;
  for (std::size_t j = 0; j < f.size(); ++j)
    for (const auto &[k, v] : f[j]) {
      Key a = k;
      a.insert(a.begin(), unsigned(j));
      add(out, a, v);
    }
  return out;
}
F unflatten(const P &p, std::size_t n) {
  F f(n);
  for (const auto &[k, v] : p)
    f.at(k[0])[Key(k.begin() + 1, k.end())] = v;
  return f;
}
unsigned degree(const P &p) {
  unsigned d = 0;
  for (const auto &[k, v] : p)
    d = std::max(d, std::accumulate(k.begin(), k.end(), 0u));
  return d;
}
Q evaluate(const P &p, const std::vector<Q> &x) {
  Q out = 0;
  for (const auto &[k, v] : p) {
    Q a = v;
    for (std::size_t j = 0; j < x.size(); ++j)
      for (unsigned e = 0; e < k[j]; ++e)
        a *= x[j];
    out += a;
  }
  return out;
}
void limitsValid(const Limits &l) {
  if (!l.maximumBasis || l.maximumBasis > 64 || !l.maximumOperations ||
      l.maximumOperations > 100000 || !l.maximumDegree || l.maximumDegree > 64)
    throw std::invalid_argument(
        "Invalid algebra budgets (basis<=64, operations<=100000, degree<=64)");
}
} // namespace dynsyslab::core::lie::detail

namespace dynsyslab::core::lie {
using namespace detail;
std::string binaryRational(double value) {
  if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite coefficient");
  if (value == 0) return "0";
  int exponent;
  const auto mantissa=std::frexp(std::abs(value),&exponent);
  boost::multiprecision::cpp_int integer=static_cast<unsigned long long>(std::ldexp(mantissa,53));
  Q q(integer); exponent-=53;
  if(exponent>=0)q*=boost::multiprecision::cpp_int(1)<<exponent;
  else q/=boost::multiprecision::cpp_int(1)<<(-exponent);
  return str(value<0?-q:q);
}
std::vector<Polynomial> independentPolynomials(const std::vector<Polynomial>& polynomials) {
  if(polynomials.empty())return {};
  const auto n=polynomials[0].dimension;
  Span span; std::vector<Polynomial> result;
  for(const auto& p:polynomials){
    if(p.dimension!=n)throw std::invalid_argument("Independent polynomial dimension mismatch");
    auto exact=read(p);if(span.insert(exact))result.push_back(write(exact,n));
  }
  return result;
}
Polynomial differentiate(const Polynomial &p, std::size_t j) {
  if (j >= p.dimension)
    throw std::invalid_argument("Derivative axis out of range");
  return write(derivative(read(p), j), p.dimension);
}
Polynomial multiply(const Polynomial &a, const Polynomial &b) {
  if (a.dimension != b.dimension)
    throw std::invalid_argument("Product dimension mismatch");
  return write(product(read(a), read(b)), a.dimension);
}
Polynomial lieDerivative(const VectorField &f, const Polynomial &p) {
  if (f.size() != p.dimension)
    throw std::invalid_argument("Lie derivative dimension mismatch");
  return write(action(read(f), read(p)), p.dimension);
}
VectorField bracket(const VectorField &x, const VectorField &y) {
  return write(commutator(read(x), read(y)));
}
Polynomial poissonBracket(const Polynomial &a, const Polynomial &b,
                          const std::vector<std::string> &tensor) {
  auto n = a.dimension;
  if (n != b.dimension || tensor.size() != n * n)
    throw std::invalid_argument("Poisson tensor shape mismatch");
  auto p = read(a), q = read(b);
  P out;
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      Q c = rational(tensor[i * n + j]);
      if (c != -rational(tensor[j * n + i]))
        throw std::invalid_argument("Poisson tensor must be antisymmetric");
      axpy(out, product(derivative(p, i), derivative(q, j)), c);
    }
  return write(out, n);
}
Polynomial flowBracket(const VectorField &f, const Polynomial &a,
                       const Polynomial &b) {
  if (f.size() != a.dimension || a.dimension != b.dimension)
    throw std::invalid_argument("Flow bracket shape mismatch");
  auto x = read(f);
  auto p = read(a), q = read(b);
  P out = product(p, action(x, q));
  axpy(out, product(q, action(x, p)), -1);
  return write(out, a.dimension);
}
std::size_t rankAt(const std::vector<VectorField> &fields,
                   const std::vector<std::string> &point) {
  std::vector<Q> x;
  for (const auto &a : point)
    x.push_back(rational(a));
  Span s;
  for (const auto &f : fields) {
    if (f.size() != x.size())
      throw std::invalid_argument("Rank point dimension mismatch");
    auto v = read(f);
    P row;
    for (std::size_t j = 0; j < v.size(); ++j)
      add(row, {unsigned(j)}, evaluate(v[j], x));
    s.insert(row);
  }
  return s.rows.size();
}
} // namespace dynsyslab::core::lie
