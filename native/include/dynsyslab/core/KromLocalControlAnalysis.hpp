#pragma once
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/LoopMargins.hpp>

namespace dynsyslab::core {
struct KromLocalControlOptions {
  std::vector<double> point;
  std::string measurement, controlMap{"constant"};
  double statePenalty{5000}, inputPenalty{.01}, processCovariance{1},
      measurementCovariance{1e-4};
  double holdingInput{};
  std::size_t frequencyPoints{1200};
};
struct KromLocalControlAnalysis {
  StateSpaceModel surrogate, physical, compensator;
  std::vector<double> gain, physicalGain, estimatorGain;
  std::vector<double> physicalEquilibriumResidual;
  // Affine consistency at the operating point. A stable frozen deviation
  // loop does not remove these generator/readout biases.
  std::vector<double> liftedAtPoint, generatorAtPoint;
  std::vector<double> measuredAtPoint, decodedMeasurementAtPoint,
      decodedStateAtPoint;
  LoopRobustness surrogateLqr, physicalLqr, surrogateLqg, physicalLqg;
};
// Frozen local LTI loops. Physical feedback uses K*J_Psi at the operating
// point; valid as the exact feedback linearization when the lifted tracking
// error is zero. Does not certify a periodic, saturated, delayed, sampled, or
// scheduled controller.
KromLocalControlAnalysis
analyzeKromLocalControl(const KoopmanModel &, const KromLocalControlOptions &);
} // namespace dynsyslab::core
