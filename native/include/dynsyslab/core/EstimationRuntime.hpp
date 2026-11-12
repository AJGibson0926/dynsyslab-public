#pragma once

#include <dynsyslab/core/Estimation.hpp>
#include <memory>

namespace dynsyslab::core {

struct KoopmanEstimate {
    double time{};
    std::vector<double> lifted;      // KROM real split layout
    std::vector<double> measurement; // y0 + Cy * lifted
    std::vector<double> output;      // x0 + Cx * lifted (optional)
    std::vector<double> innovation;  // measured y - estimated y
};

// Continuous-gain KLQE driven by sampled, held measurements and known inputs.
// No physical state, measurement function, controller, or solver is required.
class KoopmanLQEObserver {
public:
    KoopmanLQEObserver(const KoopmanModel& model, KoopmanLQEDesign design,
        std::string controlRepresentation, std::vector<double> initialLifted,
        std::vector<double> measurementOffset = {},
        std::vector<double> outputOffset = {});
    // At the first call only initialize the clock and measurement hold.
    // Later calls advance [previousTime,time) using the PREVIOUS measurement
    // and the input ACTUALLY APPLIED over that interval. Current y is then
    // latched for the next interval. This is not a discrete Kalman update.
    KoopmanEstimate observe(double time, const std::vector<double>& measurement,
        const std::vector<double>& previousAppliedInput);
    const KoopmanModel& model() const noexcept { return *model_; }
    const KoopmanEstimate& estimate() const noexcept { return estimate_; }

private:
    std::shared_ptr<const KoopmanModel> model_; // immutable deployment snapshot
    KoopmanLQEDesign design_;
    std::vector<double> offsetY_, offsetX_, previousMeasurement_;
    std::vector<double> baseGenerator_, offsetInput_, inputOperators_;
    std::vector<double> heldResponse_, heldTransition_, heldInput_;
    double heldDt_{-1};
    bool initialized_{false}, bilinear_{false};
    KoopmanEstimate estimate_;
};
} // namespace dynsyslab::core
