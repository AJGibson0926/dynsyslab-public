#include "BubbleSupport.hpp"
#include "DynamicalKnowledgeInternal.hpp"
#include "AnalyticSystems.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core::detail {
namespace {

class MarmottantNative final : public ODESystem {
public:
    explicit MarmottantNative(SystemConfiguration configuration)
        : configuration_(std::move(configuration)), descriptor_(descriptorFor("marmottant")) {
        BubbleMedium reference{1000.0, 1.0e-3, 0.0725, 0.0, 1.095};
        medium_ = bubbleMedium(configuration_, reference);
        if (configuration_.preset == "matlab_cleanup_long") {
            radius0_ = 0.82e-6; shellViscosity_ = 7.2e-9; criticalTension_ = 0.13;
        }
        radius0_ = resolvedParameter(configuration_, "equilibrium_radius_m", radius0_);
        ambient_ = resolvedParameter(configuration_, "ambient_pressure_pa", ambient_);
        shellViscosity_ = resolvedParameter(configuration_, "shell_viscosity_kg_s", shellViscosity_);
        bucklingRadius_ = resolvedParameter(configuration_, "buckling_radius_m", radius0_);
        criticalTension_ = resolvedParameter(configuration_, "critical_surface_tension_n_m", criticalTension_);
        shellElasticity_ = resolvedParameter(configuration_, "shell_elasticity_n_m", shellElasticity_);
        liquidSoundSpeed_ = resolvedParameter(configuration_, "liquid_sound_speed_m_s", liquidSoundSpeed_);
        gasSoundSpeed_ = resolvedParameter(configuration_, "gas_sound_speed_m_s", gasSoundSpeed_);
        if (!(radius0_ > 0.0 && ambient_ > 0.0 && shellViscosity_ > 0.0
              && bucklingRadius_ > 0.0 && criticalTension_ > 0.0 && shellElasticity_ > 0.0))
            throw std::invalid_argument("Marmottant physical parameters must be positive");
        initializeGroups();
        initial_ = selectedInitial(configuration_, {1.0, 0.0});
        if (initial_.size() != 2) throw std::invalid_argument("Marmottant initial state must have size 2");
    }

    void resetRuntimeState() override { ruptured_ = false; }
    const SystemDescriptor& descriptor() const override { return descriptor_; }
    const SystemConfiguration& configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return initial_; }
    DynamicalKnowledge dynamicalKnowledge() const override {
        return builtInDynamicalKnowledge(descriptor_, configuration_);
    }
    std::vector<AnalyticSolutionFamily> analyticSolutions() const override {
        return bubbleLimitAnalyticSolutions(configuration_);
    }
    void rhs(double time, const double* state, double* derivative) const override {
        const auto scaling = coordinates();
        const double physical = configuration_.driver.value(
            scaling.physicalTimeAvailable ? scaling.physicalTime(time) : time);
        rhsWithInput(time, state, scaling.physicalInputAvailable
            ? scaling.modelInput(physical) : physical, derivative);
    }
    void rhsWithInput(double, const double* state, double input, double* out) const override {
        const double radius = state[0], velocity = state[1];
        const double sigma = surfaceTension(radius, true);
        const double inertial = -1.5 * velocity * velocity;
        const double gas = (1.0 - 3.0 * medium_.kappa * liquidMach_ * velocity)
            * euler_ * gasPressureRatio_ * std::pow(radius, -3.0 * medium_.kappa);
        const double laplace = euler_ * (2.0 * sigma / (ambient_ * radius0_)) / radius;
        const double viscous = euler_ * viscosity_ * velocity / radius;
        const double shell = euler_ * shellViscosityGroup_ * velocity / (radius * radius);
        out[0] = velocity;
        out[1] = (inertial + gas - (euler_ + laplace + viscous + shell)) / radius
            - euler_ * input / radius;
    }
    std::vector<double> jacobian(double time, const double* state) const override {
        const bool held = ruptured_;
        std::vector<double> matrix(4), plus(state, state + 2), minus(plus), fp(2), fm(2);
        try {
            for (std::size_t column = 0; column < 2; ++column) {
                const double step = 1.0e-6 * std::max(1.0, std::abs(state[column]));
                plus[column] += step; minus[column] -= step;
                ruptured_ = held; rhs(time, plus.data(), fp.data());
                ruptured_ = held; rhs(time, minus.data(), fm.data());
                for (std::size_t row = 0; row < 2; ++row)
                    matrix[row * 2 + column] = (fp[row] - fm[row]) / (2.0 * step);
                plus[column] = state[column]; minus[column] = state[column];
            }
        } catch (...) {
            ruptured_ = held;
            throw;
        }
        ruptured_ = held;
        return matrix;
    }
    std::vector<double> inputMatrix(const double* state) const override {
        return {0.0, -euler_ / state[0]};
    }
    bool invalid(double, const double* state) const override {
        return genericInvalid(state, 2) || state[0] <= 0.0;
    }
    std::vector<Observable> modelObservables() const override {
        static const char* names[]{"surface_tension","shell_viscosity",
            "inertial_kinetic_energy","internal_gas_pressure","laplace_pressure",
            "viscous_damping","momentum_balance"};
        static const char* labels[]{"Surface tension","Shell viscosity",
            "Inertial kinetic energy","Internal gas pressure","Laplace pressure",
            "Viscous damping","Momentum balance"};
        std::vector<Observable> result;
        for (std::size_t index=0;index<7;++index)
            result.emplace_back(ObservableDescriptor{names[index],labels[index],
                index==0 ? "N/m" : "dimensionless","dimensionless",{},false},2,1,
                [this,index](const double* state) {
                    return Observable::Values{observableTerms(state)[index]};
                });
        return result;
    }
    std::map<std::string,std::string> observableAliases() const override {
        return {{"r","radius"},{"R","radius"},{"rdot","radial_velocity"},
                {"Rdot","radial_velocity"}};
    }

private:
    void initializeGroups() {
        const auto scaling = bubbleScaling(medium_, radius0_, ambient_, 0.5);
        const double sigma0 = shellElasticity_ * (std::pow(radius0_ / bucklingRadius_, 2.0) - 1.0);
        const double pg0 = ambient_ + 2.0 * sigma0 / radius0_;
        liquidMach_ = liquidSoundSpeed_ == 0.0 ? 0.0 : radius0_ / (liquidSoundSpeed_ * scaling.timeScale);
        gasMach_ = gasSoundSpeed_ == 0.0 ? 0.0 : radius0_ / (3.0 * gasSoundSpeed_ * scaling.timeScale);
        static_cast<void>(gasMach_);
        euler_ = ambient_ * scaling.timeScale * scaling.timeScale
            / (medium_.density * radius0_ * radius0_);
        viscosity_ = 4.0 * medium_.viscosity / (ambient_ * scaling.timeScale);
        gasPressureRatio_ = pg0 / ambient_;
        shellViscosityGroup_ = 4.0 * shellViscosity_ / (ambient_ * radius0_ * scaling.timeScale);
        buckling_ = bucklingRadius_ / radius0_;
        rupture_ = bucklingRadius_ * std::sqrt(1.0 + medium_.surfaceTension / shellElasticity_) / radius0_;
        breaking_ = bucklingRadius_ * std::sqrt(1.0 + criticalTension_ / shellElasticity_) / radius0_;
        configuration_.dimensionlessGroups = {
            {"liquid_mach", liquidMach_}, {"gas_mach", gasMach_}, {"euler", euler_},
            {"viscosity", viscosity_}, {"gas_pressure_ratio", gasPressureRatio_},
            {"shell_viscosity", shellViscosityGroup_}, {"buckling_radius_ratio", buckling_},
            {"rupture_radius_ratio", rupture_}, {"breaking_radius_ratio", breaking_},
            {"natural_frequency_hz", scaling.naturalFrequency},
            {"characteristic_time_s", scaling.timeScale}};
        configuration_.dimensionlessGroupUnits = {
            {"natural_frequency_hz", "Hz"}, {"characteristic_time_s", "s"}};
        configuration_.coordinateScaling = {scaling.timeScale, ambient_, true, true,
            "tau", "p_control/p0", "Pa"};
    }

    double surfaceTension(double radius, bool update) const {
        if (update) {
            if (radius < rupture_) ruptured_ = false;
            else if (radius >= breaking_) ruptured_ = true;
        }
        if (radius < buckling_) return 0.0;
        if (radius <= rupture_)
            return shellElasticity_ * (std::pow(radius / buckling_, 2.0) - 1.0);
        if (radius < breaking_)
            return ruptured_ ? medium_.surfaceTension
                : shellElasticity_ * (std::pow(radius / buckling_, 2.0) - 1.0);
        if (radius > breaking_) return medium_.surfaceTension;
        return 0.0;
    }

    std::vector<double> observableTerms(const double* state) const {
        const double radius=state[0], velocity=state[1];
        const double sigma=surfaceTension(radius,false);
        const double inertial=-1.5*velocity*velocity;
        const double gas=(1.0-3.0*medium_.kappa*liquidMach_*velocity)*euler_
            *gasPressureRatio_*std::pow(radius,-3.0*medium_.kappa);
        const double laplace=euler_*(2.0*sigma/(ambient_*radius0_))/radius;
        const double viscous=euler_*viscosity_*velocity/radius;
        const double shell=euler_*shellViscosityGroup_*velocity/(radius*radius);
        return {sigma,shell,inertial,gas,laplace,viscous,
                inertial+gas-(euler_+laplace+viscous+shell)};
    }

    SystemConfiguration configuration_;
    SystemDescriptor descriptor_;
    BubbleMedium medium_;
    std::vector<double> initial_;
    double radius0_{0.975e-6}, ambient_{101.0e3}, shellViscosity_{15.0e-9};
    double bucklingRadius_{0.975e-6}, criticalTension_{1.0}, shellElasticity_{1.0};
    double liquidSoundSpeed_{1480.0}, gasSoundSpeed_{0.0};
    double liquidMach_{}, gasMach_{}, euler_{}, viscosity_{}, gasPressureRatio_{};
    double shellViscosityGroup_{}, buckling_{}, rupture_{}, breaking_{};
    mutable bool ruptured_{false};
};

}  // namespace

std::unique_ptr<ODESystem> makeMarmottant(const SystemConfiguration& configuration) {
    return std::make_unique<MarmottantNative>(configuration);
}

}  // namespace dynsyslab::core::detail
