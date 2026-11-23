#include "Exact.hpp"

namespace dynsyslab::core::lie {
// Triangle inequality + mean-value integral on a convex box. This is a
// deterministic bound for KNOWN polynomial residuals, not a sampled EDMD bound.
// Motivated by proportional residual sectors in Strasser et al., SafEDMD,
// arXiv:2402.03145, Corollary 3.2; no SDP/LMI certificate is claimed here.
BoxBound residualBoxBound(const std::vector<Polynomial> &residual,
                          const std::vector<std::string> &radii) {
  using namespace detail;
  if (radii.empty() || radii.size() > 16 || residual.empty() ||
      residual.size() > 64)
    throw std::invalid_argument("Invalid residual box dimensions");
  std::vector<Q> rho;
  for (const auto &s : radii) {
    auto q = rational(s);
    if (q <= 0)
      throw std::invalid_argument("Box radii must be positive");
    rho.push_back(q);
  }
  auto bound = [&](const P &p) {
    Q value = 0;
    for (const auto &[k, c] : p) {
      Q v = abs(c);
      for (std::size_t j = 0; j < rho.size(); ++j)
        for (unsigned e = 0; e < k[j]; ++e)
          v *= rho[j];
      value += v;
    }
    return value;
  };
  Q uniform = 0;
  std::vector<Q> columnSums(rho.size());
  bool vanishes = true;
  for (const auto &p : residual) {
    if (p.dimension != rho.size())
      throw std::invalid_argument("Residual box dimension mismatch");
    auto poly = read(p);
    uniform += bound(poly);
    if (poly.count(Key(rho.size())))
      vanishes = false;
    for (std::size_t j = 0; j < rho.size(); ++j)
      columnSums[j] += bound(derivative(poly, j));
  }
  return {str(uniform),
          str(*std::max_element(columnSums.begin(), columnSums.end())),
          vanishes};
}
} // namespace dynsyslab::core::lie
