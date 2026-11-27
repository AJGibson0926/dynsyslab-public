#include "BubbleSupport.hpp"

#include <cmath>

namespace dynsyslab::core::detail {

std::unique_ptr<ODESystem> makeRPE(SystemConfiguration configuration) {
    const auto medium = bubbleMedium(configuration);
    const double radius0 = resolvedParameter(configuration, "equilibrium_radius_m", 10.0e-6);
    const double ambient = resolvedParameter(configuration, "ambient_pressure_pa", 100000.0);
    const double inputSign = resolvedParameter(configuration, "pressure_input_sign", -1.0);
    const auto scaling = bubbleScaling(medium, radius0, ambient);
    const double gas = scaling.gasPressure / scaling.pressureScale;
    const double vapor = medium.vaporPressure / scaling.pressureScale;
    const double external = ambient / scaling.pressureScale;
    const double capillarity = 2.0 * medium.surfaceTension
        / (radius0 * scaling.pressureScale);
    const double viscosity = 4.0 * medium.viscosity
        / (scaling.timeScale * scaling.pressureScale);
    configuration.dimensionlessGroups = {
        {"gas_pressure_ratio", gas}, {"vapor_pressure_ratio", vapor},
        {"ambient_pressure_ratio", external}, {"capillarity", capillarity},
        {"viscosity", viscosity}, {"natural_frequency_hz", scaling.naturalFrequency},
        {"characteristic_time_s", scaling.timeScale},
        {"pressure_scale_pa", scaling.pressureScale}};
    configuration.dimensionlessGroupUnits = {
        {"natural_frequency_hz", "Hz"}, {"characteristic_time_s", "s"},
        {"pressure_scale_pa", "Pa"}};
    configuration.coordinateScaling = {scaling.timeScale, ambient, true, true,
        "tau", "p_control/p0", "Pa"};
    configuration.coordinateScaling.physicalStatePerModelState = {scaling.radius, scaling.radius / scaling.timeScale};
    configuration.coordinateScaling.physicalStateUnits = {"m", "m/s"};
    auto rhs = [medium, gas, vapor, external, capillarity, viscosity, inputSign]
        (double, const double* x, double input, double* out) {
        const double radius = x[0], velocity = x[1];
        const double numerator = gas * std::pow(radius, -3.0 * medium.kappa) + vapor
            - external - capillarity / radius - viscosity * velocity / radius
            - 1.5 * velocity * velocity + inputSign * external * input;
        out[0] = velocity;
        out[1] = numerator / radius;
    };
    auto jacobian = [medium, gas, vapor, external, capillarity, viscosity, inputSign,
                     driver = configuration.driver, coordinateScaling=configuration.coordinateScaling]
        (double time, const double* x) {
        const double radius = x[0], velocity = x[1];
        const double physicalTime = coordinateScaling.physicalTimeAvailable
            ? coordinateScaling.physicalTime(time) : time;
        const double physicalInput = driver.value(physicalTime);
        const double modelInput = coordinateScaling.physicalInputAvailable
            ? coordinateScaling.modelInput(physicalInput) : physicalInput;
        const double pressure = inputSign * external * modelInput;
        const double balance = gas * std::pow(radius, -3.0 * medium.kappa) + vapor
            - external + pressure - capillarity / radius - viscosity * velocity / radius;
        const double derivative = -3.0 * medium.kappa * gas
            * std::pow(radius, -3.0 * medium.kappa - 1.0)
            + capillarity / (radius * radius)
            + viscosity * velocity / (radius * radius);
        return std::vector<double>{0.0, 1.0,
            (radius * derivative - balance + 1.5 * velocity * velocity) / (radius * radius),
            (-viscosity / radius - 3.0 * velocity) / radius};
    };
    auto input = [inputSign, external](const double* x) {
        return std::vector<double>{0.0, inputSign * external / x[0]};
    };
    const auto r=AnalyticExpression::state(2,0), v=AnalyticExpression::state(2,1);
    const auto inertial=1.5*v.power(2);
    const auto internal=gas*r.power(-3.0*medium.kappa);
    const auto laplace=capillarity/r, viscous=viscosity*v/r;
    const auto momentum=internal-inertial+vapor-external-laplace-viscous;
    // Distance-scaled spherical-monopole pressure: rho*(2*R*Rdot^2+R^2*Rddot).
    // In native coordinates momentum = r*r'', not r''. The acceleration
    // contribution must therefore be r*momentum, preserving both radius factors.
    const auto echo=radius0*scaling.pressureScale*r*(2.0*v.power(2)+momentum);
    const std::vector<AnalyticExpression> terms{inertial,internal,laplace,viscous,momentum,echo};
    const std::vector<AnalyticExpression> drift{v,momentum/r};
    const char* names[]{"inertial_kinetic_energy","internal_gas_pressure","laplace_pressure",
        "viscous_damping","momentum_balance","far_field_acoustic_echo"};
    const char* labels[]{"Inertial kinetic energy","Internal gas pressure","Laplace pressure",
        "Viscous damping","Momentum balance","Far-field acoustic echo"};
    std::vector<Observable> observables;
    // First integral from the RPE integrating factor r^3 (the elliptic-solution
    // references are recorded in AnalyticRPE.cpp). For arbitrary kappa, use the
    // logarithmic gas-energy limit at kappa=1. Subtract H(1,0) so equilibrium is 0.
    const auto gasEnergy=medium.kappa==1.0 ? -gas*r.unary("log")
        : gas/(3.0*(medium.kappa-1.0))*(r.power(3.0*(1.0-medium.kappa))-1.0);
    const auto energy=0.5*r.power(3)*v.power(2)+gasEnergy
        +(external-vapor)/3.0*(r.power(3)-1.0)+capillarity/2.0*(r.power(2)-1.0);
    const auto loss=-viscosity*r*v.power(2);
    const auto inputPower=inputSign*external*r.power(2)*v;
    for(const auto& item:std::vector<std::pair<ObservableDescriptor,AnalyticExpression>>{
        {{"mechanical_energy","Mechanical energy above equilibrium","dimensionless","dimensionless",{},false},energy},
        {{"energy_dissipation","Viscous energy rate","dimensionless","dimensionless",{},false},loss},
        {{"control_power_coefficient","Energy rate per control input","dimensionless","dimensionless",{},false},inputPower}})
        observables.push_back(Observable::fromExpressions(item.first,{item.second}));
    for (std::size_t index=0;index<6;++index) {
        const std::string units=index==5 ? "Pa m" : "dimensionless";
        observables.push_back(Observable::fromExpressions(
            ObservableDescriptor{names[index],labels[index],units,"dimensionless",{},false},
            {terms[index]}));
    }
    const auto echoDot=echo.flowDerivative(drift);
    const auto echoDdot=echoDot.flowDerivative(drift);
    observables.push_back(Observable::fromExpressions(
        {"far_field_acoustic_echo_dot","\\frac{de_0}{d\\tau}","Pa m","dimensionless",{},false},
        {echoDot}));
    observables.push_back(Observable::fromExpressions(
        {"far_field_acoustic_echo_ddot","\\frac{d^{2}e_0}{d\\tau^{2}}","Pa m","dimensionless",{},false},
        {echoDdot}));
    observables.push_back(Observable::fromExpressions(
        {"far_field_acoustic_echo_dot_dr","\\frac{\\partial \\dot e_0}{\\partial r}","Pa m","dimensionless",{},false},
        {echoDot.derivative(0)}));
    observables.push_back(Observable::fromExpressions(
        {"far_field_acoustic_echo_dot_dv","\\frac{\\partial \\dot e_0}{\\partial v}","Pa m","dimensionless",{},false},
        {echoDot.derivative(1)}));
    return std::make_unique<FunctionalSystem>(descriptorFor("rpe"), configuration,
        selectedInitial(configuration, {1.0, 0.0}), rhs, jacobian, input,
        [](double, const double* x) {
            return genericInvalid(x, 2) || x[0] <= 0.01;
        }, std::vector<double>{1.0, 0.0}, std::move(observables),
        std::map<std::string,std::string>{{"r","radius"},{"R","radius"},
        {"rdot","radial_velocity"},{"Rdot","radial_velocity"},{"echo","far_field_acoustic_echo"},{"energy","mechanical_energy"},
        {"echo_dot","far_field_acoustic_echo_dot"},{"echo_ddot","far_field_acoustic_echo_ddot"},
        {"echo_dot_dr","far_field_acoustic_echo_dot_dr"},{"echo_dot_dv","far_field_acoustic_echo_dot_dv"}},
        std::map<std::string,std::vector<std::string>>{{"echo",{"far_field_acoustic_echo"}},
        {"echo_jet2",{"far_field_acoustic_echo","far_field_acoustic_echo_dot"}}},drift,
        std::vector<AnalyticExpression>{AnalyticExpression::constant(0.,2),inputSign*external/r});
}

}  // namespace dynsyslab::core::detail
