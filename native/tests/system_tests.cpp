#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message) + ": got " + std::to_string(actual)
            + ", expected " + std::to_string(expected));
}

std::unique_ptr<dynsyslab::core::ODESystem> system(const char* identifier,
                                                   const char* preset = "") {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = identifier;
    configuration.preset = preset;
    return dynsyslab::core::createSystem(configuration);
}

void testCatalogAndBasicSystems() {
    std::set<std::string> identifiers;
    for(const auto& descriptor:dynsyslab::core::systemCatalog()) {
        require(!descriptor.identifier.empty()&&descriptor.stateDimension>0,"invalid catalog descriptor");
        require(identifiers.insert(descriptor.identifier).second,"duplicate catalog identifier");
    }
    for(const auto* required:{"rpe","duffing","nonlinear_pendulum","van_der_pol","kuramoto",
        "morse_oscillator","isotonic_oscillator","poschl_teller","cosh_oscillator","radial_twist",
        "kepler_planar","euler_top","coupled_duffing_hamiltonian","toda_lattice"})
        require(identifiers.count(required)==1,"required system missing from catalog");
    auto duffing = system("duffing", "double_well_reference");
    const double state[]{0.5, -0.25}; double derivative[2]{};
    duffing->rhsWithInput(0.0, state, 0.2, derivative);
    close(derivative[0], -0.25, 1.0e-14, "Duffing position derivative");
    close(derivative[1], 0.625, 1.0e-14, "Duffing acceleration");
    require(duffing->jacobian(0.0, state).size() == 4, "Duffing Jacobian shape");
    require(duffing->initialState()[0] == 0.1, "Duffing free default IC");
    require(duffing->controlledInitialState()[0] == 1.0, "Duffing controlled equilibrium IC");
    const auto observables=duffing->observableCatalog();
    close(observables.at("x").evaluateReal(state)[0],0.5,1.0e-14,"observable alias");
    close(observables.at("total_energy").evaluateReal(state)[0],
          0.5*0.25*0.25-0.5*0.5*0.5+0.25*std::pow(0.5,4),1.0e-14,
          "Duffing total-energy observable");
    require(observables.at("state").hasAnalyticJacobian(),"state observable Jacobian");

    auto quadratic = system("quadratic_saddle_node", "exact_fold");
    const double scalarState[]{1.0}; double scalarDerivative[1]{};
    quadratic->rhs(0.0, scalarState, scalarDerivative);
    close(scalarDerivative[0], -1.0, 1.0e-14, "quadratic fold vector field");
    close(quadratic->jacobian(0.0, scalarState)[0], -2.0, 1.0e-14,
          "quadratic fold Jacobian");

    auto cubic = system("cubic_bifurcation", "symmetric_pitchfork");
    cubic->rhs(0.0, scalarState, scalarDerivative);
    close(scalarDerivative[0], -1.0, 1.0e-14, "cubic pitchfork vector field");
    close(cubic->jacobian(0.0, scalarState)[0], -3.0, 1.0e-14,
          "cubic pitchfork Jacobian");

    auto radial = system("radial_limit_cycle", "unit_cycle");
    const double cycleState[]{1.0, 0.0}; double cycleDerivative[2]{};
    radial->rhs(0.0, cycleState, cycleDerivative);
    close(cycleDerivative[0], 0.0, 1.0e-14, "radial cycle radial derivative");
    close(cycleDerivative[1], 1.0, 1.0e-14, "radial cycle angular derivative");
    const auto radialJacobian = radial->jacobian(0.0, cycleState);
    close(radialJacobian[0], -2.0, 1.0e-14, "radial cycle transverse Jacobian");
    close(radialJacobian[1], -1.0, 1.0e-14, "radial cycle rotation Jacobian");
}

void testMechanicalSystems() {
    auto pendulum = system("nonlinear_pendulum", "earth_small_angle");
    const double state[]{0.2, 0.1}; double derivative[2]{};
    pendulum->rhsWithInput(0.0, state, 0.3, derivative);
    close(derivative[0], 0.1, 1.0e-14, "pendulum velocity");
    close(derivative[1], -9.80665*std::sin(0.2) - 0.002 + 0.3, 1.0e-12,
          "pendulum acceleration");
    auto fput = system("fput", "alpha_fput_32");
    require(fput->descriptor().stateDimension == 64, "FPUT preset dimension");
    dynsyslab::core::SystemConfiguration expanded;
    expanded.identifier = "fput"; expanded.parameters["sites"] = 32;
    fput = dynsyslab::core::createSystem(expanded);
    require(fput->descriptor().stateDimension == 64, "FPUT configured dimension");
    auto kuramoto = system("kuramoto");
    require(kuramoto->initialState().size() == 16, "Kuramoto initial shape");
}

void testBiologicalSystems() {
    auto neuron = system("hindmarsh_rose", "canonical_bursting");
    const auto initial = neuron->initialState(); std::vector<double> derivative(3);
    neuron->rhsWithInput(0.0, initial.data(), 0.0, derivative.data());
    close(derivative[0], 0.125, 1.0e-12, "Hindmarsh-Rose x derivative");
    close(derivative[1], 0.75, 1.0e-12, "Hindmarsh-Rose y derivative");
    auto hodgkin = system("hodgkin_huxley", "squid_axon_rest");
    const auto hhInitial = hodgkin->initialState(); std::vector<double> hhDerivative(4);
    hodgkin->rhsWithInput(0.0, hhInitial.data(), 0.0, hhDerivative.data());
    for (double value : hhDerivative) require(std::isfinite(value), "non-finite HH derivative");
    require(hodgkin->jacobian(0.0, hhInitial.data()).size() == 16, "HH Jacobian shape");
}

void testBubbleSystems() {
    auto rpe = system("rpe", "air_water_20c");
    const double equilibrium[]{1.0, 0.0}; double derivative[2]{};
    rpe->rhsWithInput(0.0, equilibrium, 0.0, derivative);
    close(derivative[0], 0.0, 1.0e-14, "RPE velocity");
    close(derivative[1], 0.0, 1.0e-10, "RPE equilibrium acceleration");
    require(rpe->inputMatrix(equilibrium)[1] < 0.0, "RPE pressure sign");
    require(rpe->controlledInitialState() == std::vector<double>({1.0, 0.0}),
            "RPE controlled dimensionless equilibrium IC");
    require(rpe->configuration().parameters.at("equilibrium_radius_m") == 1.0e-5,
            "RPE resolved radius");
    require(rpe->configuration().parameterUnits.at("ambient_pressure_pa") == "Pa",
            "RPE pressure units");
    require(rpe->configuration().dimensionlessGroups.count("characteristic_time_s") == 1,
            "RPE characteristic time");
    const auto scaling = rpe->coordinates();
    require(scaling.physicalTimeAvailable && scaling.physicalInputAvailable,
            "RPE physical/model scaling contract");
    close(scaling.physicalInputPerModelInput, 100000.0, 1.0e-9,
          "RPE p_control/p0 scale");
    close(scaling.physicalTime(1.0),
          rpe->configuration().dimensionlessGroups.at("characteristic_time_s"),
          1.0e-18, "RPE model-time conversion");
    dynsyslab::core::SystemConfiguration drivenConfiguration;
    drivenConfiguration.identifier = "rpe";
    drivenConfiguration.preset = "air_water_20c";
    drivenConfiguration.driver.amplitude = 100000.0;
    drivenConfiguration.driver.frequencyHz =
        1.0 / scaling.secondsPerModelTime;
    auto drivenRpe = dynsyslab::core::createSystem(drivenConfiguration);
    double drivenDerivative[2]{}, explicitDerivative[2]{};
    drivenRpe->rhs(0.0, equilibrium, drivenDerivative);
    drivenRpe->rhsWithInput(0.0, equilibrium, 1.0, explicitDerivative);
    close(drivenDerivative[1], explicitDerivative[1], 1.0e-10,
          "RPE built-in physical driver conversion");
    require(!dynsyslab::core::equationsFingerprint(rpe->configuration()).empty(),
            "RPE equation fingerprint");
    const auto rpeObservables=rpe->observableCatalog();
    require(rpeObservables.contains("far_field_acoustic_echo"),"RPE echo observable");
    require(std::isfinite(rpeObservables.at("far_field_acoustic_echo")
        .evaluateReal(equilibrium)[0]),"finite RPE echo observable");
    for(double kappa:{1.,4./3.,1.5})for(double mu:{0.,.001}) {
        dynsyslab::core::SystemConfiguration c;c.identifier="rpe";c.parameters={{"kappa",kappa},{"mu_pa_s",mu}};
        const auto bubble=dynsyslab::core::createSystem(c);const auto cat=bubble->observableCatalog();
        close(cat.at("mechanical_energy").evaluateReal(equilibrium)[0],0,1e-12,"RPE equilibrium energy origin");
        for(double r:{.3,1.,2.4})for(double v:{-3.,.2})for(double u:{-.1,.05}){
            const double state[]{r,v};double f[2];bubble->rhsWithInput(0,state,u,f);
            const auto gradient=cat.at("mechanical_energy").jacobian(state);
            const double rate=cat.at("energy_dissipation").evaluateReal(state)[0]+u*cat.at("control_power_coefficient").evaluateReal(state)[0];
            close((gradient[0]*f[0]+gradient[1]*f[1]).real(),rate,1e-10,"RPE exact energy balance");
        }
    }
    auto marmottant = system("marmottant", "matlab_cleanup_linear");
    marmottant->rhs(0.0, equilibrium, derivative);
    require(std::isfinite(derivative[1]), "Marmottant derivative");
    require(marmottant->coordinates().modelInputUnits == "p_control/p0",
            "Marmottant control scaling");
    auto reddy = system("reddy_szeri", "matlab_cleanup_reference");
    const auto initial = reddy->initialState(); std::vector<double> redDerivative(12);
    reddy->rhsWithInput(0.0, initial.data(), 0.0, redDerivative.data());
    for (double value : redDerivative) require(std::isfinite(value), "non-finite Reddy derivative");
    require(reddy->inputMatrix(initial.data()).size() == 12, "Reddy input shape");
    require(reddy->coordinates().physicalTimeAvailable
            && reddy->coordinates().physicalInputAvailable,
            "Reddy physical/model scaling contract");
    const auto reddyObservables=reddy->observableCatalog();
    const auto coordinates=reddyObservables.coordinatePreset("c_k5");
    require(coordinates.outputDimension()==2,"Reddy observable coordinates");
    require(coordinates.jacobian(initial.data()).size()==24,"Reddy coordinate Jacobian");
}

}  // namespace

int main() {
    try {
        testCatalogAndBasicSystems();
        testMechanicalSystems();
        testBiologicalSystems();
        testBubbleSystems();
        std::cout << "dynsyslab_system_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_system_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
