#pragma once
// Post-processing methods for invariant geometry (docs/architecture/LCS_AND_SSM.md):
//   ftle  FTLE fields, their ridges and variational LCS of the source's system (LcsPostprocess.cpp)
//   ssm   spectral submanifolds and reduced dynamics, from the vector field or from the sources'
//         trajectories (SsmPostprocess.cpp)
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runFtle(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                              const std::filesystem::path& root, const std::string& label);
// orbit_geometry: cross-orbit injectivity and closure of coordinate records (OrbitGeometryPostprocess.cpp)
std::filesystem::path runOrbitGeometry(const std::vector<std::filesystem::path>& sources, const dynsyslab::core::json::Value& settings,
                                       const std::filesystem::path& root, const std::string& label);
// run_invariant: a conserved polynomial from run identity (RunInvariantPostprocess.cpp)
std::filesystem::path runRunInvariant(const std::vector<std::filesystem::path>& sources, const dynsyslab::core::json::Value& settings,
                                      const std::filesystem::path& root, const std::string& label);
std::filesystem::path runSsm(const std::vector<std::filesystem::path>& sources, const dynsyslab::core::json::Value& settings,
                             const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
