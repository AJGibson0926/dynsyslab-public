#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

[[nodiscard]] std::vector<double> multiply(
    const std::vector<double>& left,std::size_t rows,std::size_t inner,
    const std::vector<double>& right,std::size_t columns);

[[nodiscard]] std::vector<std::complex<double>> multiply(
    const std::vector<std::complex<double>>& left,std::size_t rows,
    std::size_t inner,const std::vector<std::complex<double>>& right,
    std::size_t columns);

}  // namespace dynsyslab::numerics
