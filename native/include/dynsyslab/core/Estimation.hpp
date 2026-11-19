#pragma once

#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/Observability.hpp>

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

// reconstructionSplitMatrix is declared in Observability.hpp (linear analysis).

struct ContinuousLQESolution {
    std::vector<double> covariance;          // state x state
    std::vector<double> gain;                // state x measurement
    std::vector<std::complex<double>> observerPoles;
    double relativeResidual{};
    bool detectable{};
    std::size_t stateDimension{};
    std::size_t measurementDimension{};
};

[[nodiscard]] ContinuousLQESolution solveContinuousLQE(
    const std::vector<double>& generator, const std::vector<double>& measurement,
    const std::vector<double>& processCovariance,
    const std::vector<double>& measurementCovariance,
    std::size_t stateDimension, std::size_t measurementDimension,
    bool requireDetectable = true);

struct KoopmanLQEDesign {
    std::string measurementObservable;
    std::string stateObservable;
    std::vector<double> generator;
    std::vector<double> measurementMatrix;
    std::vector<double> stateReconstructionMatrix;
    std::vector<double> measurementOffset;   // y = y0 + C_y z (reconstruction intercept; zeros if none)
    std::vector<double> outputOffset;        // x = x0 + C_x z
    ContinuousLQESolution solution;
    std::size_t physicalStateDimension{};
    std::size_t outputDimension{}; // optional decoded output, need not be a physical state
    std::string modelIdentity;
};

// Measurement-only design. No physical state or decoder is required to estimate
// lifted coordinates; an optional named reconstruction supplies output values.
[[nodiscard]] KoopmanLQEDesign designKoopmanLQEFromMeasurements(
    const KoopmanModel& model, const std::string& measurementObservable,
    const std::vector<double>& processCovariance,
    const std::vector<double>& measurementCovariance,
    const std::string& outputObservable = {}, bool requireDetectable = true);

[[nodiscard]] KoopmanLQEDesign designKoopmanLQE(
    const KoopmanModel& model, const std::string& measurementObservable,
    const std::string& stateObservable,
    const std::vector<double>& processCovariance,
    const std::vector<double>& measurementCovariance,
    bool requireDetectable = true);

struct KoopmanLQEResult {
    std::vector<double> time;
    std::vector<double> liftedEstimate;      // samples x lifted split state
    std::vector<double> stateEstimate;       // samples x physical state
    std::vector<double> measurementEstimate;// samples x measurement
    std::vector<double> innovation;          // samples x measurement
    std::vector<double> stateError;          // optional samples x physical state
    std::vector<double> stateRmse;
    std::vector<double> innovationRmse;
    std::size_t samples{};
    std::size_t liftedDimension{};
    std::size_t stateDimension{};
    std::size_t measurementDimension{};
};

// Measurements are zero-order held between samples, matching the retained KLQE workflow.
[[nodiscard]] KoopmanLQEResult runKoopmanLQEFromLifted(
    const KoopmanModel& model, const KoopmanLQEDesign& design,
    const std::vector<double>& time, const std::vector<double>& measurements,
    const std::vector<double>& initialRealLiftedValues,
    const std::vector<double>& trueOutputs = {});

// Physical-state initialization adapter, retained for established clients.
[[nodiscard]] KoopmanLQEResult runKoopmanLQE(
    const KoopmanModel& model, const KoopmanLQEDesign& design,
    const std::vector<double>& time, const std::vector<double>& measurements,
    const std::vector<double>& initialPhysicalState,
    const std::vector<double>& trueStates = {});

}  // namespace dynsyslab::core
