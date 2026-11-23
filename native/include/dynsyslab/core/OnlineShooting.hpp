#pragma once
#include <dynsyslab/core/ShootingIdentification.hpp>
#include <optional>

namespace dynsyslab::core {
struct OnlineShootingOptions {
    double sampleTime{};
    std::size_t firstFitIndex{1},fitEvery{1},maximumSamples{1000000};
    bool retainOriginalStarts{false};
};
struct OnlineShootingSample {
    double time{};
    std::size_t index{};
    bool ready{};
    std::vector<double> state; // empty until the first past-window fit
    std::optional<double> informationTime;
    std::optional<ShootingIdentificationResult> fit;
};
using ShootingStartProvider=std::function<std::vector<ShootingStart>(
    const std::vector<double>& times,const std::vector<double>& measurements)>;

// Causal expanding-window shooting observer. Each call receives ONLY the new
// sensor sample and the command held during the preceding interval. The observer
// has no plant reference, true initial state, controller, or future input data.
// A first-fit start provider, when supplied, sees only that acquired prefix.
// Standard single shooting/Gauss-Newton: Nocedal & Wright (2006), ch.10.
// Not a Kalman filter, process-noise model, or real-time computation guarantee.
class OnlineShootingObserver {
public:
    OnlineShootingObserver(ODESystem& model,std::vector<double> measurementVariances,
        std::vector<ShootingStart>,InitialStateMap,InputAwareReadout,
        OnlineShootingOptions,IntegrationOptions = {},
        numerics::NonlinearLeastSquaresOptions = {},ShootingStartProvider = {});
    OnlineShootingSample observe(double time,const std::vector<double>& measurement,
                                double previousAppliedInput);
    // Does not mutate the observer or assimilate data. Caller supplies the
    // newly commanded constant input for this prospective interval.
    NativeSolution forecast(double finalTime,double appliedInput,std::size_t intervals) const;
private:
    ODESystem& model_;
    std::vector<double> variances_,time_,observed_,commands_,state_;
    std::vector<ShootingStart> starts_,original_;
    InitialStateMap initialize_;
    InputAwareReadout read_;
    OnlineShootingOptions options_;
    IntegrationOptions solver_;
    numerics::NonlinearLeastSquaresOptions optimizer_;
    ShootingStartProvider provider_;
    double beforeFirst_{};
    std::optional<double> informationTime_;
};
}
