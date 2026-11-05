#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/SystemLabels.hpp>
#include <dynsyslab/core/Extension.hpp>

#include "SystemSupport.hpp"
#include "DynamicalKnowledgeInternal.hpp"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <mutex>

namespace dynsyslab::core {
namespace {

std::mutex& providerMutex() {
    static std::mutex value;
    return value;
}

std::vector<SystemCatalogProvider>& catalogProviders() {
    static std::vector<SystemCatalogProvider> value;
    return value;
}

}  // namespace

namespace {
void requirePositiveScale(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument(std::string(name) + " must be finite and positive");
}
}  // namespace

double CoordinateScaling::physicalTime(double modelTimeValue) const {
    requirePositiveScale(secondsPerModelTime, "secondsPerModelTime");
    return modelTimeValue * secondsPerModelTime;
}

double CoordinateScaling::modelTime(double physicalTimeValue) const {
    requirePositiveScale(secondsPerModelTime, "secondsPerModelTime");
    return physicalTimeValue / secondsPerModelTime;
}

double CoordinateScaling::physicalInput(double modelInputValue) const {
    requirePositiveScale(physicalInputPerModelInput, "physicalInputPerModelInput");
    return modelInputValue * physicalInputPerModelInput;
}

double CoordinateScaling::modelInput(double physicalInputValue) const {
    requirePositiveScale(physicalInputPerModelInput, "physicalInputPerModelInput");
    return physicalInputValue / physicalInputPerModelInput;
}

CoordinateScaling ODESystem::coordinates() const {
    CoordinateScaling result = configuration().coordinateScaling;
    if (descriptor().preferredCoordinates == "dimensionless") {
        if (result.modelTimeUnits == "system time") result.modelTimeUnits = "dimensionless time";
        if (result.modelInputUnits == "system input")
            result.modelInputUnits = "dimensionless system input";
    } else {
        if (result.modelTimeUnits == "system time") result.modelTimeUnits = "s";
        result.physicalTimeAvailable = true;
    }
    requirePositiveScale(result.secondsPerModelTime, "secondsPerModelTime");
    requirePositiveScale(result.physicalInputPerModelInput,
                         "physicalInputPerModelInput");
    return result;
}

ObservableCatalog ODESystem::observableCatalog() const {
    const auto& description = descriptor();
    std::vector<std::string> names = description.stateNames;
    std::vector<std::string> labels = description.stateLabels;
    if (names.size() != description.stateDimension) {
        names.resize(description.stateDimension);
        for (std::size_t index = 0; index < names.size(); ++index)
            if (names[index].empty()) names[index] = "state_" + std::to_string(index);
    }
    if (labels.size() != description.stateDimension) labels = names;
    std::vector<Observable> result;
    result.reserve(description.stateDimension + 1);
    for (std::size_t index = 0; index < description.stateDimension; ++index)
        result.push_back(Observable::stateComponent(description.stateDimension, index,
            names[index], labels[index], description.preferredCoordinates));
    result.push_back(Observable::identity(description.stateDimension, names, labels,
                                          description.preferredCoordinates));
    auto extras = modelObservables();
    result.insert(result.end(), std::make_move_iterator(extras.begin()),
                  std::make_move_iterator(extras.end()));
    auto aliases = observableAliases();
    for (std::size_t index = 0; index < names.size(); ++index) {
        aliases.emplace("x[" + std::to_string(index) + "]", names[index]);
        aliases.emplace("state_" + std::to_string(index), names[index]);
    }
    return {std::move(result), std::move(aliases), observableCoordinatePresets()};
}

double SinusoidalDriver::value(double time) const {
    constexpr double pi = 3.141592653589793238462643383279502884;
    return amplitude * std::cos(2.0 * pi * frequencyHz * time + phaseRadians);
}

const std::vector<SystemDescriptor>& systemCatalog() {
    static const std::vector<SystemDescriptor> systems = [] {
        std::vector<SystemDescriptor> catalog{
        {"harmonic_oscillator", "Harmonic oscillator", 2, "dimensional", {"position", "velocity"},
         {"x", "dx/dt"}, {{"undamped", "Undamped"}, {"damped", "Damped"}}},
        {"rpe", "Rayleigh-Plesset", 2, "dimensionless", {"radius", "radial_velocity"},
         {"r = R/R0", "dr/dtau"}, {{"air_water_0c", "Air in water at 0 C"},
         {"air_water_20c", "Air in water at 20 C"}, {"air_water_37c", "Air in water at 37 C"},
         {"air_saline_37c", "Air in saline at 37 C"}, {"air_blood_37c", "Air in blood at 37 C"},
         {"air_glycerol_20c", "Air in glycerol at 20 C"},
         {"air_glycerol_water_50pct_20c", "Air in 50% glycerol-water at 20 C"},
         {"air_silicone_oil_100cst_25c", "Air in 100 cSt silicone oil at 25 C"},
         {"co2_water_20c", "CO2 in water at 20 C"}, {"co2_honey_25c", "CO2 in honey at 25 C"}}},
        {"duffing", "Duffing", 2, "dimensionless", {"position", "velocity"},
         {"x", "dx/dt"}, {{"double_well_reference", "Canonical double well"},
         {"hardening_single_well", "Hardening single well"},
         {"gibson_undamped_klqr_2022", "Undamped KLQR benchmark"},
         {"gibson_undamped_vdp_tracking_2022", "Undamped tracking benchmark"},
         {"gibson_damped_stabilization_2022", "Damped stabilization benchmark"},
         {"korda_mezic_optimal_2020", "Koopman MPC benchmark"}}},
        {"marmottant", "Marmottant", 2, "dimensionless", {"radius", "radial_velocity"},
         {"r", "dr/dtau"}, {{"matlab_cleanup_linear", "Cleanup-suite linear shell"},
         {"matlab_cleanup_long", "Cleanup-suite long response"}}},
        {"reddy_szeri", "Reddy-Szeri", 12, "dimensionless",
         {"radius","radial_velocity","translation","translation_velocity",
          "shape_2","shape_2_velocity","shape_3","shape_3_velocity",
          "shape_4","shape_4_velocity","shape_5","shape_5_velocity"},
         {"r = R/R0","dr/dtau","x/R0","dx/dtau","a2","da2/dtau",
          "a3","da3/dtau","a4","da4/dtau","a5","da5/dtau"},
         {{"matlab_cleanup_reference", "Cleanup-suite reference"},
         {"matlab_cleanup_control", "Cleanup-suite control case"}}},
        {"van_der_pol", "Van der Pol", 2, "dimensional", {"position", "velocity"},
         {"x", "dx/dt"}, {{"weakly_nonlinear", "Weakly nonlinear"},
         {"relaxation_oscillation", "Relaxation oscillation"}}},
        {"quadratic_saddle_node", "Quadratic saddle-node unfolding", 1, "dimensionless",
         {"state"}, {"x"}, {{"exact_fold", "Exact quadratic fold"}}},
        {"cubic_bifurcation", "Cubic bifurcation unfolding", 1, "dimensionless",
         {"state"}, {"x"}, {{"symmetric_pitchfork", "Reflection-symmetric pitchfork"},
         {"imperfect_pitchfork", "Imperfect pitchfork / cusp slice"}}},
        {"radial_limit_cycle", "Radial attracting limit cycle", 2, "dimensionless",
         {"x", "y"}, {"x", "y"}, {{"unit_cycle", "Exact unit attracting cycle"}}},
        {"rossler", "Rossler", 3, "dimensional", {"x", "y", "z"}, {"x", "y", "z"},
         {{"classic_chaos", "Classic chaos"}, {"common_alternative", "Common alternative"},
         {"brunton_shavok_2021", "HAVOK reference"}, {"funnel_chaos_2021", "Funnel chaos"},
         {"b045_chaos_2013", "b=0.45 reference"}}},
        {"lorenz", "Lorenz", 3, "dimensional", {"x", "y", "z"}, {"x", "y", "z"},
         {{"lorenz_1963", "Lorenz 1963"}, {"arbabi_mezic_hdmd_2017", "HDMD reference"},
         {"near_classic_transition", "Near transition"}}},
        {"chua", "Chua circuit", 3, "dimensional", {"x", "y", "z"}, {"x", "y", "z"},
         {{"canonical_double_scroll", "Canonical double scroll"},
         {"double_scroll_seed_left", "Left-lobe seed"}}},
        {"henon_heiles", "Henon-Heiles", 4, "dimensional", {}, {},
         {{"low_energy", "Low energy"}, {"near_escape_energy", "Near escape energy"}}},
        {"hindmarsh_rose", "Hindmarsh-Rose", 3, "dimensional", {"x", "y", "z"},
         {"x", "y", "z"}, {{"canonical_bursting", "Canonical bursting"},
         {"standard_spiking", "Standard spiking"}}},
        {"hodgkin_huxley", "Hodgkin-Huxley", 4, "dimensional", {"V", "m", "h", "n"},
         {"V", "m", "h", "n"}, {{"squid_axon_rest", "Squid axon rest"},
         {"tonic_firing_10", "Tonic firing at input 10"}}},
        {"fput", "FPUT", 32, "dimensional", {}, {}, {{"alpha_fput_32", "Alpha FPUT"},
         {"beta_fput_32", "Beta FPUT"}}},
        {"kuramoto", "Kuramoto", 16, "dimensional", {}, {}, {{"weak_coupling", "Weak coupling"},
         {"synchronizing", "Synchronizing coupling"}}},
        {"nonlinear_pendulum", "Nonlinear pendulum", 2, "dimensional", {"angle", "velocity"},
         {R"(\theta)", R"(\frac{d\theta}{dt})"}, {{"earth_small_angle", "Earth small angle"},
         {"large_angle_low_damping", "Large angle, low damping"},
         {"driven_damped", "Driven damped"}, {"holmes_chaotic", "Holmes chaotic"}}},
        {"double_pendulum", "Double pendulum", 4, "dimensional", {},
         {R"(\theta_1)", R"(\frac{d\theta_1}{dt})", R"(\theta_2)", R"(\frac{d\theta_2}{dt})"},
         {{"equal_arm_demo", "Equal-arm demonstration"},
         {"classic_chaotic_release", "Classic chaotic release"},
         {"unequal_lab_rig", "Unequal laboratory rig"}}},
        {"spring_coupled_pendulums", "Spring-coupled pendulums", 4, "dimensional", {},
         {R"(\theta_1)", R"(\frac{d\theta_1}{dt})", R"(\theta_2)", R"(\frac{d\theta_2}{dt})"},
         {{"symmetric_weak_coupling", "Symmetric weak coupling"},
         {"strong_nonlinear_coupling", "Strong nonlinear coupling"},
         {"detuned_energy_exchange", "Detuned energy exchange"}}},
        };
        for (auto& descriptor : catalog) {
            descriptor.research.papers = detail::builtInResearchPapers(descriptor.identifier);
            if (descriptor.identifier == "kuramoto")
                descriptor.stateLabels = systemStateLabels(descriptor, descriptor.stateDimension);
        }
        return catalog;
    }();
    return systems;
}

std::vector<SystemDescriptor> availableSystems() {
    std::vector<SystemDescriptor> result = systemCatalog();
    {
        std::lock_guard<std::mutex> lock(providerMutex());
        for (const auto provider : catalogProviders()) {
            auto supplied = provider();
            result.insert(result.end(), std::make_move_iterator(supplied.begin()),
                          std::make_move_iterator(supplied.end()));
        }
    }
    auto extensionSystems = ExtensionRegistry::instance().systems();
    result.insert(result.end(), std::make_move_iterator(extensionSystems.begin()),
                  std::make_move_iterator(extensionSystems.end()));
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        return left.identifier < right.identifier;
    });
    const auto duplicate = std::adjacent_find(result.begin(), result.end(),
        [](const auto& left, const auto& right) {
            return left.identifier == right.identifier;
        });
    if (duplicate != result.end())
        throw std::runtime_error("duplicate system identifier: " + duplicate->identifier);
    return result;
}

void registerSystemCatalogProvider(SystemCatalogProvider provider) {
    if (!provider) throw std::invalid_argument("system catalog provider is null");
    std::lock_guard<std::mutex> lock(providerMutex());
    if (std::find(catalogProviders().begin(), catalogProviders().end(), provider)
        == catalogProviders().end())
        catalogProviders().push_back(provider);
}

std::unique_ptr<ODESystem> createSystem(const SystemConfiguration& configuration) {
    SystemConfiguration resolved = configuration;
    if (resolved.equationsVersion.empty())
        resolved.equationsVersion = "dynsyslab." + resolved.identifier + "/1";
    using Factory = std::unique_ptr<ODESystem> (*)(const SystemConfiguration&);
    for (const Factory factory : {detail::makeBasicSystem, detail::makeMechanicalSystem,
                                  detail::makeBiologicalSystem, detail::makeBubbleSystem}) {
        if (auto system = factory(resolved)) return system;
    }
    if (auto system = ExtensionRegistry::instance().createSystem(resolved)) return system;
    throw std::invalid_argument("unsupported native system: " + configuration.identifier);
}

std::string equationsFingerprint(const SystemConfiguration& configuration) {
    std::ostringstream canonical;
    canonical << configuration.identifier << '\n' << configuration.preset << '\n'
              << configuration.equationsVersion << '\n';
    for (const auto& item : configuration.parameters) {
        canonical << item.first << '=' << std::setprecision(17) << item.second << ':';
        const auto unit = configuration.parameterUnits.find(item.first);
        if (unit != configuration.parameterUnits.end()) canonical << unit->second;
        canonical << '\n';
    }
    std::uint64_t hash = 1469598103934665603ull;
    for (const unsigned char item : canonical.str()) {
        hash ^= item;
        hash *= 1099511628211ull;
    }
    std::ostringstream result;
    result << std::hex << std::setw(16) << std::setfill('0') << hash;
    return result.str();
}

}  // namespace dynsyslab::core
