#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/ShootingIdentification.hpp>
#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <array>
namespace dynsyslab::workflow {
struct PeriodicFeedbackBinding {
    core::ODESystem& plant;
    core::ODESystem& model;
    core::InitialStateMap initializePlant,initializeEstimate;
    std::function<std::vector<double>(const double*)> decodePlant,decodeEstimate;
    core::InputAwareReadout measurePlant,measureEstimate;
    std::array<std::size_t,2> stateChannels;
};
core::json::Value evaluatePeriodicFeedback(const core::json::Value&,
    const math::PeriodicAtlas&,const PeriodicFeedbackBinding&,
    const std::vector<core::ShootingStart>&,core::IntegrationOptions,
    const numerics::NonlinearLeastSquaresOptions&);
}
