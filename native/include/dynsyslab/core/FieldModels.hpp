#pragma once
// Data-driven vector fields and maps beyond a single dictionary regression
// (docs/architecture/FIELD_MODELS.md). Ported from the miscellaneous-computes side-car.
//
// A field model is a JSON definition with a "kind"; fieldFromDefinition turns any of them into an
// Observable (state -> rate, or inputs -> outputs for a map), so fits, rollouts and archives treat
// every kind alike. Kinds:
//   polynomial_model  a dictionary and its coefficients (SindyDynamicsModel definition)
//   sum               the sum of several fields (a baseline and a fitted residual)
//   compact_radial    Wendland-weighted (1 - r)^4 (1 + 4 r) partition over the k nearest training
//                     centers (Wendland 1995, Adv. Comput. Math. 4:389-396), optional affine tail
//   local_tangent     sparse quadratic fields in two-dimensional local PCA tangent charts, blended
//                     by the same weights, with damping transverse to each chart
//   paired_phase      radial and angular rates of complex modal pairs from amplitude and harmonic
//                     relative-phase features (phase locking of pair j to pair 0 at winding h_j)
//   second_order      imposed kinematic rows x_{2i}' = x_{2i+1} and a fitted acceleration map
//   pair_network      second-order taps whose accelerations are shared by lag separation
//   interior          outer taps by a map, interior taps as critically damped responses to a drive
// Regression everywhere is sequentially thresholded ridge least squares on normalized moments
// (SINDy: Brunton, Proctor and Kutz 2016, PNAS 113:3932-3937; ridge: Hoerl and Kennard 1970).
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Observable.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

// One record: times, states (samples x dimension) and labels (samples x outputs: rates, or a
// map's outputs).
struct FieldRecord {
    std::vector<double> time, state, label;
    std::size_t dimension{}, outputs{};
    [[nodiscard]] std::size_t samples() const noexcept { return time.size(); }
};

// Coefficients (features x outputs) from the moments gram = [Theta, Y]^T [Theta, Y] of size
// (features + outputs)^2: features and targets normalized, ridge added, terms below `threshold`
// (in normalized units) removed and the rest refitted until the support is stable.
[[nodiscard]] std::vector<double> thresholdedLeastSquares(const std::vector<double>& gram, std::size_t features,
                                                          std::size_t outputs, double ridge, double threshold);

[[nodiscard]] Observable fieldFromDefinition(const json::Value& definition);

struct Rollout {
    std::vector<double> time, state;  // computed samples only
    bool complete{};
};
// x' = f(x) from x0 on `times` by classical RK4 with steps <= maximumStep (0: a quarter interval);
// stops at the first non-finite state (complete false).
[[nodiscard]] Rollout rolloutField(const Observable& field, const std::vector<double>& x0, const std::vector<double>& times,
                                   double maximumStep = 0);

struct MapFitOptions {
    std::string family{"polynomial"};  // core/ModelLibraries.hpp families
    std::size_t degree{2}, centers{24}, harmonics{5};
    std::vector<double> rationalWidths{0.25, 1.0, 4.0};
    double ridge{1e-6}, threshold{0.003};
    bool timeWeighted{false};          // trapezoid quadrature weights
};
// Y = Theta(X) C over a library frozen on the records; with a baseline field the fit is of the
// residual and the result is their sum.
[[nodiscard]] json::Value fitMap(const std::vector<FieldRecord>& records, const MapFitOptions& options,
                                 const json::Value& baseline = {});

struct CompactRadialOptions {
    std::size_t stride{40}, neighbors{10};
    double asinhScale{0};
    bool affineTail{false};
    double ridge{0.01}, threshold{1e-5};
    bool timeWeighted{false};
    double speedScale{0};  // > 0: rows weighted by 1 / hypot(1, |y| / speedScale)
};
[[nodiscard]] json::Value fitCompactRadialField(const std::vector<FieldRecord>& records, const CompactRadialOptions& options = {});

struct LocalTangentOptions {
    std::size_t stride{36}, samples{160}, neighbors{14};
    double decay{20}, ridge{2e-4}, threshold{0.006};
};
[[nodiscard]] json::Value fitLocalTangentField(const std::vector<FieldRecord>& records, const LocalTangentOptions& options = {});

struct PairedPhaseOptions {
    std::vector<int> harmonics;            // winding numbers per pair; empty: from the first record
    bool phase{true};
    std::vector<std::size_t> weakWindows;  // integral rows over these sample counts
    double instantWeight{1}, ridge{3e-5}, threshold{0.003};
};
[[nodiscard]] json::Value fitPairedPhaseField(const std::vector<FieldRecord>& records, const PairedPhaseOptions& options = {});

// Records carry the full rates (samples x dimension) as labels.
[[nodiscard]] json::Value fitSecondOrderField(const std::vector<FieldRecord>& records, const MapFitOptions& options);
[[nodiscard]] json::Value fitPairNetwork(const std::vector<FieldRecord>& records, const std::vector<double>& lags,
                                         const std::vector<double>& weights, const MapFitOptions& options, double tolerance = 1e-12);
struct InteriorOptions {
    std::vector<double> gamma;          // one value or one per interior tap
    MapFitOptions outer, drive;
    bool waveformDrive{false};          // fit the drive from waveforms (fitDampedDrive)
    bool feedback{false};
    MapFitOptions feedbackOptions;
};
[[nodiscard]] json::Value fitInteriorField(const std::vector<FieldRecord>& records, const InteriorOptions& options);

// Filtered-library waveform fit of x'' + 2 g x' + g^2 x = g^2 phi(u): the responses of every library
// column through the critically damped filter (exact for inputs held at interval midpoints), less
// each record's free transient, regressed on the observed positions and velocities. Returns the
// drive f = g^2 phi (the interior kind's forcing).
struct DriveRecord {
    std::vector<double> time, input, position, velocity;  // input samples x m, position/velocity samples x p
    std::size_t inputs{}, outputs{};
};
[[nodiscard]] json::Value fitDampedDrive(const std::vector<DriveRecord>& records, double gamma, const MapFitOptions& options);

// Waveform objective: sample-by-channel weights combining a global-RMS and fixed-interval local-RMS
// normalized squared error with trapezoid quadrature (samples x channels).
[[nodiscard]] std::vector<double> waveformWeights(const std::vector<double>& time, const std::vector<double>& observations,
                                                  std::size_t channels, const std::vector<double>& globalWeights,
                                                  const std::vector<double>& localWeights, double duration,
                                                  const std::vector<double>& floors, std::vector<double> globalFloors = {});
// Rotate coefficient directions (coefficients x p) by the SVD of the training response (rows x p),
// keeping resolved combinations; returns the new directions (coefficients x kept) and fills
// `responses` (rows x kept).
[[nodiscard]] std::vector<double> resolvedDirections(const std::vector<double>& response, std::size_t rows,
                                                     const std::vector<double>& directions, std::size_t coefficients,
                                                     std::size_t perturbations, std::vector<double>& responses,
                                                     double relativeCutoff = 0.001, double singularFloor = 0.015, double scale = 0.035);

struct CalibrationOptions {
    double duration{1};
    std::size_t startsPerRecord{2}, maximumEvaluations{40};
    double coefficientPenalty{1e-4}, maximumStep{0};
    double changeBound{50};  // relative coefficient change allowed in each direction
    // Waveform weighting of the residual (empty global weights: plain RMS scaling).
    std::vector<double> globalWeights, localWeights, floors;
    double waveformDuration{0};
};
struct CalibrationResult {
    json::Value definition;  // the calibrated polynomial model, or the original when rejected
    bool accepted{};
    double initialLoss{}, candidateLoss{};
    std::size_t shots{}, terms{}, evaluations{};
};
// Refit the active coefficients of a polynomial model against short free trajectories of the
// records (Levenberg-Marquardt; Nocedal and Wright 2006, chapter 10).
[[nodiscard]] CalibrationResult calibrateField(const json::Value& polynomialModel, const std::vector<FieldRecord>& records,
                                               const CalibrationOptions& options);

}  // namespace dynsyslab::core
