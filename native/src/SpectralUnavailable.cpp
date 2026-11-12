#include <dynsyslab/core/Spectral.hpp>

#include <stdexcept>

namespace dynsyslab::core {
namespace { [[noreturn]] void unavailable() {
    throw std::runtime_error("native linear algebra was disabled in this build");
} }
SVDResult singularValueDecomposition(const std::vector<double>&, std::size_t,
                                     std::size_t, std::size_t) { unavailable(); }
HankelDMDResult fitHankelDMD(HankelMatrix,const std::vector<std::string>&,double,const HankelDMDOptions&) { unavailable(); }
DMDResult dynamicModeDecomposition(const std::vector<double>&, std::size_t,
                                   std::size_t, double, const DMDOptions&) { unavailable(); }
HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>&, std::size_t,
    std::size_t, double, const HankelDMDOptions&) { unavailable(); }
HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>&,std::size_t,
    const std::vector<std::string>&,double,const HankelDMDOptions&) { unavailable(); }
HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>&,std::size_t,
    std::size_t,double,const Observable&,const HankelDMDOptions&) { unavailable(); }
HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>&,std::size_t,
    std::size_t,double,const std::vector<Observable>&,const HankelDMDOptions&) { unavailable(); }
}  // namespace dynsyslab::core
