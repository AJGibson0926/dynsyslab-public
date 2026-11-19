#pragma once
// mode_selection: choose the eigenfunction subset of a candidate K-ROM pool for
// echo-only KLQE -> KLQR (docs/architecture/MODE_SELECTION.md).
//
// Every quantity is offline. Online, the estimator measures the echo alone
// (y = e0 + C_e psi); psi is never evaluated on live data.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/KoopmanModel.hpp>

#include <complex>
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail::mode_selection {

using Complex = std::complex<double>;

// One record, thinned, with psi for every pool mode.
struct Record {
    std::string path;
    bool forced{};
    bool readout{false};   // used to fit C_e, C_x (readout_records)
    double dt{};                      // after thinning
    std::size_t samples{}, modes{}, stateDimension{};
    std::vector<double> time, echo, states, input;  // states: samples x stateDimension
    std::vector<Complex> psi;                        // samples x modes
    std::string derivativeEvidence;                  // how e' entered (coordinate records)
};

// A real mode or a conjugate pair: the unit of selection.
struct Unit {
    std::vector<std::size_t> modes;  // pool indices (1, or 2 for a pair: [upper, lower])
    Complex lambda;                  // continuous rate of modes[0]
    bool pair{};
    double scale{1.0};               // rms |psi| on unforced records
    Complex inputMap;                // constant B_psi of modes[0], regressed on forced records
    double inputMapResidual{1.0};    // forced: rms(d - B g) / rms(d), d = psi(t+h) - e^{lambda h} psi(t), g = int e^{lambda(h-s)} u
    double invariance{1.0};          // unforced: rms(psi(t+h) - e^{lambda h} psi(t)) / rms(psi(t+h) - psi(t))
    double processIntensity{};       // unforced miss variance / h / scale^2: the KLQE process-noise intensity
    double echoShare{};              // echo variance explained by this unit alone (R^2)
};

struct Options {
    std::string measurement{"far_field_acoustic_echo"};
    double sensorNoise{0.01};      // relative to the echo's standard deviation
    double ridge{1e-10};
    double stabilityMargin{1e-6};
};

// Linear echo and state readouts of one subset, fitted on the readout records.
struct Readouts {
    std::vector<double> echo;      // 1 x k (scaled split coordinates)
    double echoOffset{};
    std::vector<double> state;     // d x k
    std::vector<double> stateOffset;
    double echoError{};            // rms residual / std(e)
    std::vector<double> stateError;// per state component, rms residual / std
    double echoStd{};
    std::vector<double> stateStd;
};

struct Score {
    Readouts readouts;
    double predictedStateError{};  // rms over components of sqrt(eps_x^2 + (C_x P C_x^T)_ii / var x_i)
    std::vector<double> predictedStateErrors;
    double inputMapResidual{};     // worst unit
    double stateReach{};           // sqrt(lambda_min(D^-1 C_x W_c C_x^T D^-1)), larger is better
    double observabilityBalance{}; // lambda_min / lambda_max of W_o(A, C_e)
    double pbhMargin{};            // smallest relative PBH observability margin
    double hsvEchoMin{}, hsvEchoMax{}, hsvStateMin{}, hsvStateMax{};
    bool valid{true};
    std::string failure;
};

std::vector<Record> loadRecords(const core::KoopmanModel& model,
    const std::vector<std::filesystem::path>& sources, const json::Value& settings, const Options& options);
std::vector<Unit> buildUnits(const core::KoopmanModel& model, const std::vector<Record>& records,
    const json::Value& settings, const Options& options);
Score scoreSubset(const std::vector<Unit>& units, const std::vector<std::size_t>& subset,
    const std::vector<Record>& records, const Options& options);
Readouts fitReadouts(const std::vector<Unit>& units, const std::vector<std::size_t>& subset,
    const std::vector<Record>& records, double ridge);

}  // namespace dynsyslab::workflow::postprocess_detail::mode_selection
