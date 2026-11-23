#include <dynsyslab/core/KoopmanModel.hpp>
#include <stdexcept>

namespace dynsyslab::core {
// Chain rule on the existing frozen observable and polynomial contracts.
// See Goswami/Paley (CDC 2017), lifted input map; no new derivative backend.
std::vector<std::complex<double>>
KoopmanModel::stateJacobian(const double *x) const {
  if (requiresExternalCoordinates())
    throw std::invalid_argument(
        "Instantaneous state Jacobian is unavailable for history coordinates");
  if (scattered_)
    return scattered_->gradient(x);
  const auto q = coordinateObservable().evaluateReal(x);
  const auto h = coordinateObservable().jacobian(x);
  const auto d = stateDimension(), k = coordinateDimension(),
             m = liftedDimension();
  std::vector<std::complex<double>> result(m * d);
  for (std::size_t i = 0; i < m; ++i) {
    const auto row = polynomialEigenfunctions().gradient(i, q.data());
    for (std::size_t a = 0; a < d; ++a)
      for (std::size_t b = 0; b < k; ++b)
        result[i * d + a] += row[b] * h[b * d + a];
  }
  return result;
}
std::vector<double> KoopmanModel::stateJacobianSplit(const double *x) const {
  const auto complex = stateJacobian(x);
  const auto d = stateDimension(), m = liftedDimension(), n = splitDimension();
  std::vector<double> out(n * d);
  for (std::size_t a = 0; a < d; ++a) {
    std::vector<std::complex<double>> column(m);
    for (std::size_t i = 0; i < m; ++i)
      column[i] = complex[i * d + a];
    const auto real = splitValues(column);
    for (std::size_t i = 0; i < n; ++i)
      out[i * d + a] = real[i];
  }
  return out;
}
} // namespace dynsyslab::core
