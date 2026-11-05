#include <dynsyslab/core/SINDy.hpp>
#include <stdexcept>
namespace dynsyslab::core {
namespace { [[noreturn]] void unavailable() { throw std::runtime_error("native linear algebra was disabled in this build"); } }
KoopmanSINDyResult fitKoopmanEigenfunctions(const std::vector<double>&, std::size_t,
    std::size_t, const std::vector<double>&, const std::vector<std::complex<double>>&,
    const std::vector<double>&, const KoopmanSINDyOptions&) { unavailable(); }
KoopmanSINDyResult fitKoopmanEigenfunctionsInCoordinates(const std::vector<double>&,
    std::size_t,std::size_t,const std::vector<double>&,const std::vector<std::complex<double>>&,
    const Observable&,const std::vector<double>&,const KoopmanSINDyOptions&) { unavailable(); }
}  // namespace dynsyslab::core
