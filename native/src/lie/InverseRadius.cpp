#include "Exact.hpp"

namespace dynsyslab::core::lie::detail {
// Exact certificate for inverse-radius second-order mechanics. Introduce w=1/r:
// F=(v,a(r,w)+b(r,w)*v-3*w*v^2/2,-v*w^2), G=(0,w,0).
// The surface r*w=1 is invariant. H=[F,G] has H_w=w^3 and d_v H_v=2*w^2.
// Thus [H,w^n*d_v]=(n-2)*w^(n+2)*d_v. Odd n never makes the coefficient zero.
// The witnesses remain independent ON r*w=1, w>0, not merely in ambient R^3.
// Project derivation (RPE control review, 2026-11-21); finite-representation
// obstruction: Katayama/Susuki, arXiv:2611.20438, Theorem 1.
bool inverseRadiusCascade(const std::vector<F> &generators, const Span &basis) {
  F input(3);
  input[1][{0, 0, 1}] = 1;
  if (!basis.reduce(flatten(input)).second.empty())
    return false;
  for (const auto &f : generators) {
    if (f.size() != 3 || f[0] != P{{{0, 1, 0}, Q(1)}} ||
        f[2] != P{{{0, 1, 2}, Q(-1)}})
      continue;
    P quadratic;
    bool valid = true;
    for (const auto &[powers, value] : f[1]) {
      if (powers[1] > 2)
        valid = false;
      if (powers[1] == 2)
        quadratic[powers] = value;
    }
    if (valid && quadratic == P{{{0, 2, 1}, Q(-3) / 2}})
      return true;
  }
  return false;
}
} // namespace dynsyslab::core::lie::detail
