#include "BubbleSupport.hpp"
#include "DynamicalKnowledgeInternal.hpp"
#include "AnalyticSystems.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace dynsyslab::core::detail {
namespace {

struct ReddyGroups {
    double liquidMach{}, gasMach{}, euler{}, capillarity{}, viscosity{}, gasPressure{};
    double pressureAmplitude{}, pressureDerivativeAmplitude{}, omega{}, waveNumber{};
};

class ReddySzeriNative final : public ODESystem {
public:
    explicit ReddySzeriNative(SystemConfiguration configuration)
        : configuration_(std::move(configuration)), descriptor_(descriptorFor("reddy_szeri")) {
        medium_ = bubbleMedium(configuration_);
        radius0_ = resolvedParameter(configuration_, "equilibrium_radius_m", 2.0e-6);
        ambient_ = resolvedParameter(configuration_, "ambient_pressure_pa", 100.0e3);
        liquidSoundSpeed_ = resolvedParameter(configuration_, "liquid_sound_speed_m_s", 0.0);
        gasSoundSpeed_ = resolvedParameter(configuration_, "gas_sound_speed_m_s", 0.0);
        translationEnabled_ = resolvedParameter(configuration_, "translation_enabled", 1.0) != 0.0;
        acousticAmplitude_ = resolvedParameter(configuration_, "acoustic_pressure_amplitude", 0.0, "Pa");
        driveFrequency_ = resolvedParameter(configuration_, "drive_frequency_hz", 850.75e3);
        const auto scaling = bubbleScaling(medium_, radius0_, ambient_);
        configuration_.coordinateScaling = {scaling.timeScale, ambient_, true, true,
            "tau", "p_control/p0", "Pa"};
        groups_ = groups();
        configuration_.dimensionlessGroups = {
            {"liquid_mach", groups_.liquidMach}, {"gas_mach", groups_.gasMach},
            {"euler", groups_.euler}, {"capillarity", groups_.capillarity},
            {"viscosity", groups_.viscosity}, {"gas_pressure_ratio", groups_.gasPressure},
            {"pressure_amplitude", groups_.pressureAmplitude},
            {"pressure_derivative_amplitude", groups_.pressureDerivativeAmplitude},
            {"angular_frequency", groups_.omega}, {"wave_number", groups_.waveNumber},
            {"natural_frequency_hz", scaling.naturalFrequency},
            {"characteristic_time_s", scaling.timeScale}};
        configuration_.dimensionlessGroupUnits = {
            {"natural_frequency_hz", "Hz"}, {"characteristic_time_s", "s"}};
        if (translationEnabled_) descriptor_.observableTimeDerivatives["translation"]={"translation_velocity",false};
        const std::vector<double> standard{0.524736182055959, 0.0, 0.0, 0.1,
            0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0};
        const std::vector<double> control{1.0, 0.0, 0.0, 0.01,
            0.01, 0.0, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0};
        initial_ = selectedInitial(configuration_,
            configuration_.preset == "matlab_cleanup_control" ? control : standard);
        if (initial_.size() != 12) throw std::invalid_argument("Reddy-Szeri initial state must have size 12");
    }

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
    void rhsWithInput(double time, const double* x, double input, double* out) const override {
        const double radius = x[0], velocity = x[1];
        std::fill(out, out + 12, 0.0);
        const double radialAcceleration = volumeAcceleration(x);
        const auto control = inputMatrix(x);
        out[0] = velocity;
        out[1] = radialAcceleration + control[1] * input;
        double relative = 0.0, relativePrime = 0.0;
        if (translationEnabled_) {
            relative = x[3] + groups_.pressureAmplitude * std::sin(groups_.omega * time);
            const double acceleration = -(3.0 * velocity / radius
                + 4.5 * groups_.euler * groups_.viscosity / (radius * radius)) * relative
                + 3.0 * groups_.pressureDerivativeAmplitude * std::cos(groups_.omega * time);
            relativePrime = -(groups_.pressureDerivativeAmplitude
                * std::cos(groups_.omega * time) + acceleration);
            out[2] = x[3];
            out[3] = acceleration;
        }
        const double c = 3.0 * velocity / radius;
        auto k = [&](double n) {
            return 0.5 * groups_.euler * groups_.capillarity
                * (n - 1.0) * (n + 1.0) * (n + 2.0) / std::pow(radius, 3.0)
                - (n - 1.0) * out[1] / radius
                - 4.5 * n * std::pow(n + 1.0, 2.0)
                    / ((2.0*n - 1.0) * (2.0*n + 3.0)) * std::pow(relative / radius, 2.0);
        };
        auto d = [&](double n) {
            return 9.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)
                / (4.0*(2.0*n - 3.0)*(2.0*n - 1.0)) * std::pow(relative / radius, 2.0);
        };
        auto f = [&](double n) {
            return 3.0*(n + 1.0)*(2.0*(n - 1.0)*relative*velocity + n*radius*relativePrime)
                / (2.0*(2.0*n - 1.0)*radius*radius);
        };
        auto h = [&](double n) {
            return 3.0*std::pow(n + 1.0, 2.0)*(-2.0*relative*velocity + radius*relativePrime)
                / (2.0*(2.0*n + 3.0)*radius*radius);
        };
        auto j = [&](double n) {
            return 9.0*std::pow(n + 1.0, 2.0)*(n + 2.0)*(n + 3.0)
                / (4.0*(2.0*n + 3.0)*(2.0*n + 5.0)) * std::pow(relative / radius, 2.0);
        };
        auto q = [&](double n) { return 1.5*(n + 1.0)*relative / radius; };
        out[4] = x[5];
        out[5] = -k(2.0)*x[4] - c*x[5] + h(2.0)*x[6] + j(2.0)*x[8]
            - 2.25*relative*relative/radius + control[5]*input;
        out[6] = x[7];
        out[7] = f(3.0)*x[4] + q(3.0)*x[5] - k(3.0)*x[6] - c*x[7]
            + h(3.0)*x[8] + j(3.0)*x[10] + control[7]*input;
        out[8] = x[9];
        out[9] = d(4.0)*x[4] + f(4.0)*x[6] + q(4.0)*x[7] - k(4.0)*x[8]
            - c*x[9] + h(4.0)*x[10] + control[9]*input;
        out[10] = x[11];
        out[11] = d(5.0)*x[6] + f(5.0)*x[8] + q(5.0)*x[9] - k(5.0)*x[10]
            - c*x[11] + control[11]*input;
    }
    std::vector<double> jacobian(double time, const double* state) const override {
        std::vector<double> matrix(144), plus(state, state + 12), minus(plus), fp(12), fm(12);
        for (std::size_t column = 0; column < 12; ++column) {
            const double step = 1.0e-6 * std::max(1.0, std::abs(state[column]));
            plus[column] += step; minus[column] -= step;
            rhs(time, plus.data(), fp.data()); rhs(time, minus.data(), fm.data());
            for (std::size_t row = 0; row < 12; ++row)
                matrix[row * 12 + column] = (fp[row] - fm[row]) / (2.0 * step);
            plus[column] = state[column]; minus[column] = state[column];
        }
        return matrix;
    }
    std::vector<double> inputMatrix(const double* x) const override {
        const double radius = x[0];
        const double kappa = radius > 1.0 ? 1.0 : medium_.kappa;
        const double denominator = radius*radius
            + groups_.euler*groups_.liquidMach*groups_.viscosity*radius
            + 3.0*kappa*groups_.euler*groups_.gasPressure*groups_.gasMach*groups_.liquidMach
                * std::pow(radius, -(3.0*kappa - 1.0)/2.0);
        const double radialDenominator = denominator / radius;
        std::vector<double> matrix(12, 0.0);
        matrix[1] = -groups_.euler / radialDenominator;
        matrix[5] = -groups_.euler*x[4] / denominator;
        matrix[7] = -2.0*groups_.euler*x[6] / denominator;
        matrix[9] = -3.0*groups_.euler*x[8] / denominator;
        matrix[11] = -4.0*groups_.euler*x[10] / denominator;
        return matrix;
    }
    bool invalid(double, const double* state) const override {
        return genericInvalid(state, 12, 1.0e10) || state[0] <= 0.01;
    }
    std::vector<Observable> modelObservables() const override {
        static const char* names[]{"compressibility_loss","radial_acceleration","c",
            "k2","k3","k4","k5","g1"};
        static const char* labels[]{"Compressibility loss","Radial acceleration","c",
            "k2","k3","k4","k5","g1 = r rdot"};
        std::vector<Observable> result;
        for (std::size_t index=0;index<8;++index)
            result.emplace_back(ObservableDescriptor{names[index],labels[index],"dimensionless",
                "dimensionless",{},false},12,1,[this,index](const double* state) {
                    return Observable::Values{observableTerms(state)[index]};
                });
        return result;
    }
    std::map<std::string,std::string> observableAliases() const override {
        return {{"r","radius"},{"R","radius"},{"rdot","radial_velocity"},
            {"Rdot","radial_velocity"},{"x","translation"},
            {"xdot","translation_velocity"},{"Rdotdot","radial_acceleration"}};
    }
    std::map<std::string,std::vector<std::string>>
    observableCoordinatePresets() const override {
        return {{"radial_state",{"radius","radial_velocity"}},
                {"r_rdot",{"radius","radial_velocity"}},
                {"c_k2",{"c","k2"}},{"c_k5",{"c","k5"}}};
    }

private:
    ReddyGroups groups() const {
        const auto scale = bubbleScaling(medium_, radius0_, ambient_);
        const double liquidMach = liquidSoundSpeed_ == 0.0 ? 0.0
            : radius0_ / (liquidSoundSpeed_ * scale.timeScale);
        const double gasMach = gasSoundSpeed_ == 0.0 ? 0.0
            : radius0_ / (3.0 * gasSoundSpeed_ * scale.timeScale);
        const double euler = ambient_ * scale.timeScale * scale.timeScale
            / (medium_.density * radius0_ * radius0_);
        constexpr double twoPi = 6.283185307179586476925286766559;
        const double omega = twoPi * driveFrequency_ * scale.timeScale;
        const double pressure = -euler * liquidMach * acousticAmplitude_;
        const double wavelength = liquidSoundSpeed_ == 0.0 ? 0.0 : liquidSoundSpeed_ / driveFrequency_;
        return {liquidMach, gasMach, euler,
            2.0*medium_.surfaceTension/(ambient_*radius0_),
            4.0*medium_.viscosity/(ambient_*scale.timeScale),
            scale.gasPressure/ambient_, pressure, -omega*euler*liquidMach*acousticAmplitude_,
            omega, wavelength == 0.0 ? 0.0 : twoPi*radius0_/wavelength};
    }
    double volumeAcceleration(const double* x) const {
        const double radius = x[0], velocity = x[1];
        const double kappa = radius > 1.0 ? 1.0 : medium_.kappa;
        const double t1 = -1.5 * velocity * velocity;
        const double t2 = groups_.euler * groups_.gasPressure
            * (1.0 - 3.0*kappa*groups_.liquidMach*velocity) * std::pow(radius, -3.0*kappa);
        const double t3 = 3.0*kappa*groups_.euler*groups_.gasPressure*groups_.gasMach*velocity
            * (1.5*(kappa + 1.0)*groups_.liquidMach*velocity - 1.0)
            * std::pow(radius, -1.5*(kappa + 1.0));
        const double laplace = groups_.euler*groups_.capillarity
            * (groups_.liquidMach*velocity - 1.0) / radius;
        const double viscous = groups_.euler*groups_.viscosity*(velocity/radius)
            * (1.0 - groups_.liquidMach*velocity);
        const double denominator = radius + groups_.euler*groups_.liquidMach*groups_.viscosity
            + 3.0*kappa*groups_.euler*groups_.gasPressure*groups_.gasMach*groups_.liquidMach
                * std::pow(radius, -(3.0*kappa + 1.0)/2.0);
        return (t1 + t2 + t3 - groups_.euler - viscous + laplace) / denominator;
    }
    std::vector<double> observableTerms(const double* x) const {
        const double radius=x[0], velocity=x[1], kappa=medium_.kappa;
        const double t1=-1.5*velocity*velocity;
        const double t2=groups_.euler*groups_.gasPressure
            *(1.0-3.0*kappa*groups_.liquidMach*velocity)*std::pow(radius,-3.0*kappa);
        const double t3=3.0*kappa*groups_.euler*groups_.gasPressure*groups_.gasMach*velocity
            *(1.5*(kappa+1.0)*groups_.liquidMach*velocity-1.0)
            *std::pow(radius,-1.5*(kappa+1.0));
        const double laplace=groups_.euler*groups_.capillarity
            *(groups_.liquidMach*velocity-1.0)/radius;
        const double viscous=groups_.euler*groups_.viscosity*(velocity/radius)
            *(1.0-groups_.liquidMach*velocity);
        const double denominator=radius+groups_.euler*groups_.liquidMach*groups_.viscosity
            +3.0*kappa*groups_.euler*groups_.gasPressure*groups_.gasMach*groups_.liquidMach
            *std::pow(radius,-(3.0*kappa+1.0)/2.0);
        const double acceleration=(t1+t2+t3-groups_.euler-viscous+laplace)/denominator;
        const double c=3.0*velocity/radius;
        auto k=[&](double n) { return 0.5*groups_.euler*groups_.capillarity
            *(n-1.0)*(n+1.0)*(n+2.0)/std::pow(radius,3.0)-(n-1.0)*acceleration/radius; };
        return {t3/radius,acceleration,c,k(2.0),k(3.0),k(4.0),k(5.0),radius*velocity};
    }

    SystemConfiguration configuration_;
    SystemDescriptor descriptor_;
    BubbleMedium medium_;
    ReddyGroups groups_;
    std::vector<double> initial_;
    double radius0_{2.0e-6}, ambient_{100.0e3}, liquidSoundSpeed_{}, gasSoundSpeed_{};
    bool translationEnabled_{true};
    double acousticAmplitude_{}, driveFrequency_{850.75e3};
};

}  // namespace

std::unique_ptr<ODESystem> makeReddySzeri(const SystemConfiguration& configuration) {
    return std::make_unique<ReddySzeriNative>(configuration);
}

}  // namespace dynsyslab::core::detail
