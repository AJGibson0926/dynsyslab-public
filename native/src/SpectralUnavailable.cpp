#include <dynsyslab/core/Spectral.hpp>

#include <stdexcept>

namespace dynsyslab::core {
namespace { [[noreturn]] void unavailable() {
    throw std::runtime_error("native linear algebra was disabled in this build");
} }
SVDResult singularValueDecomposition(const std::vector<double>&, std::size_t,
                                     std::size_t, std::size_t) { unavailable(); }
DMDResult dynamicModeDecomposition(const std::vector<double>&, std::size_t,
                                   std::size_t, double, const DMDOptions&) { unavailable(); }
HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>&, std::size_t,
    std::size_t, double, const HankelDMDOptions&) { unavailable(); }
}  // namespace dynsyslab::core
