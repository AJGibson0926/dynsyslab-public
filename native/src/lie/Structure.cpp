#include "Exact.hpp"

namespace dynsyslab::core::lie::detail {
// Classical adjoint/derived/lower-central/center/Killing computations over Q.
// See Serre, Lie Algebras and Lie Groups (1965), and the representation-based
// embedding criterion in Katayama/Susuki, arXiv:2611.20438.
std::vector<std::vector<Q>> nullspace(std::vector<std::vector<Q>> a,
                                      std::size_t n) {
  std::vector<std::size_t> pivots;
  std::size_t row = 0;
  for (std::size_t col = 0; col < n && row < a.size(); ++col) {
    std::size_t pivot = row;
    while (pivot < a.size() && a[pivot][col] == 0)
      ++pivot;
    if (pivot == a.size())
      continue;
    std::swap(a[row], a[pivot]);
    Q scale = a[row][col];
    for (auto &x : a[row])
      x /= scale;
    for (std::size_t i = 0; i < a.size(); ++i)
      if (i != row) {
        Q v = a[i][col];
        for (std::size_t j = 0; j < n; ++j)
          a[i][j] -= v * a[row][j];
      }
    pivots.push_back(col);
    ++row;
  }
  std::vector<std::vector<Q>> out;
  for (std::size_t j = 0; j < n; ++j)
    if (std::find(pivots.begin(), pivots.end(), j) == pivots.end()) {
      std::vector<Q> v(n);
      v[j] = 1;
      for (std::size_t i = 0; i < pivots.size(); ++i)
        v[pivots[i]] = -a[i][j];
      out.push_back(v);
    }
  return out;
}
void properties(AlgebraResult &out) {
  const auto r = out.basis.size();
  std::vector<Q> c;
  for (const auto &a : out.structureConstants)
    c.push_back(rational(a));
  auto at = [&](std::size_t i, std::size_t j, std::size_t k) -> const Q & {
    return c[(i * r + j) * r + k];
  };
  auto comm = [&](const P &x, const P &y) {
    P z;
    for (const auto &[i, a] : x)
      for (const auto &[j, b] : y)
        for (std::size_t k = 0; k < r; ++k)
          add(z, {unsigned(k)}, a * b * at(i[0], j[0], k));
    return z;
  };
  Span all;
  for (std::size_t i = 0; i < r; ++i)
    all.insert(P{{{unsigned(i)}, Q(1)}});
  auto series = [&](bool derived) {
    Span current = all;
    std::vector<std::size_t> dims{r};
    while (!current.rows.empty()) {
      Span next;
      const auto &left = derived ? current : all;
      for (const auto &a : left.originals)
        for (const auto &b : current.originals)
          next.insert(comm(a, b));
      dims.push_back(next.rows.size());
      if (next.rows.size() == current.rows.size())
        break;
      current = std::move(next);
    }
    return dims;
  };
  out.derivedDimensions = series(true);
  out.lowerCentralDimensions = series(false);
  out.solvable = out.derivedDimensions.back() == 0;
  out.nilpotent = out.lowerCentralDimensions.back() == 0;
  std::vector<std::vector<Q>> centerRows;
  for (std::size_t j = 0; j < r; ++j)
    for (std::size_t k = 0; k < r; ++k) {
      std::vector<Q> row(r);
      for (std::size_t i = 0; i < r; ++i)
        row[i] = at(i, j, k);
      centerRows.push_back(std::move(row));
    }
  for (const auto &v : nullspace(centerRows, r)) {
    std::vector<std::string> row;
    for (const auto &a : v)
      row.push_back(str(a));
    out.center.push_back(std::move(row));
  }
  out.killingForm.resize(r * r, "0");
  for (std::size_t i = 0; i < r; ++i)
    for (std::size_t j = 0; j < r; ++j) {
      Q value = 0;
      for (std::size_t k = 0; k < r; ++k)
        for (std::size_t l = 0; l < r; ++l)
          value += at(i, l, k) * at(j, k, l);
      out.killingForm[i * r + j] = str(value);
    }
}
} // namespace dynsyslab::core::lie::detail
