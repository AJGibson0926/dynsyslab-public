#pragma once
// Proportional residual envelope ||r|| <= c_x ||z|| + c_u ||u|| from residual data,
// the error structure used by SafEDMD (Strasser, Worthmann et al., arXiv:2402.03145).
//
// The estimate is the exact optimum of a two-variable linear program over the
// supplied samples. It is a DATA ENVELOPE: it says nothing about unseen states
// unless the caller declares sampling assumptions. With a declaration of i.i.d.
// samples, a scenario bound on the violation probability is reported (Campi &
// Garatti, SIAM J. Optim. 19 (2008) 1211-1230). DynSysLab cannot check that
// declaration; trajectory data is generally not i.i.d. The theorem needs an
// objective fixed before sampling, so a scenario request must set both weights.
#include <dynsyslab/core/BilinearKrom.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ScenarioAssumption {
    std::string declaration;  // required: how the samples were drawn (i.i.d. from what)
    double confidence{0.999}; // 1 - beta, in (0, 1)
};

struct ResidualSectorOptions {
    double inflation{1.0};       // >= 1, multiplies both coefficients after the fit
    double stateWeight{0.0};     // LP objective weight; 0 selects the mean ||z|| (not with a scenario)
    double inputWeight{0.0};     // LP objective weight; 0 selects the mean ||u|| (not with a scenario)
    double zeroTolerance{1e-13}; // norms below this (relative to the data scale) count as zero
    std::optional<ScenarioAssumption> scenario;
};

struct ResidualSectorBound {
    bool feasible{false};   // false: a nonzero residual occurs at z = 0, u = 0
    double cState{}, cInput{};
    double inflation{1.0};
    std::size_t samples{};
    double maximumResidual{};
    std::optional<double> violationProbability;  // only with a declared scenario
    std::optional<double> confidence;
    std::string guarantee;
    std::string message;
};

struct ResidualSectorCoverage {
    std::size_t samples{}, violations{};
    double fraction{};      // fraction of samples satisfying the bound
    double worstExcess{};   // max(||r|| - bound, 0)
};

[[nodiscard]] ResidualSectorBound estimateResidualSector(const std::vector<double>& z,
    const std::vector<double>& u, const std::vector<double>& residual, std::size_t samples,
    std::size_t stateDimension, std::size_t inputDimension, const ResidualSectorOptions& options = {});
[[nodiscard]] ResidualSectorBound estimateResidualSector(const BilinearKrom& model,
    const BilinearDataset& data, const ResidualSectorOptions& options = {});
// Use on held-out data to test an estimated bound.
[[nodiscard]] ResidualSectorCoverage residualSectorCoverage(const ResidualSectorBound& bound,
    const std::vector<double>& z, const std::vector<double>& u, const std::vector<double>& residual,
    std::size_t samples, std::size_t stateDimension, std::size_t inputDimension);
// Smallest epsilon with sum_{i<d} C(N,i) eps^i (1-eps)^(N-i) <= beta (Campi & Garatti 2008).
[[nodiscard]] double scenarioViolationBound(std::size_t samples, std::size_t decisionVariables, double beta);

}  // namespace dynsyslab::core
