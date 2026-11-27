#include "SystemViewDefaults.hpp"

#include <algorithm>
#include <cctype>
#include <map>

namespace dsl {
namespace {

struct Defaults {
    std::vector<std::string> views;
    bool oblique = false;  // centered oblique camera that shows both lobes of a strange attractor
};

const std::map<std::string, Defaults>& table() {
    static const std::map<std::string, Defaults> defaults{
        {"rpe", {{"time", "phase", "fourier"}}},
        {"rayleigh-plesset", {{"time", "phase", "fourier"}}},
        {"duffing", {{"time", "phase", "fourier"}}},
        {"van_der_pol", {{"time", "phase"}}},
        {"nonlinear_pendulum", {{"time", "phase"}}},
        {"henon_heiles", {{"phase", "time"}}},
        {"double_pendulum", {{"time", "phase"}}},
        {"spring_coupled_pendulums", {{"time", "phase"}}},
        {"lorenz", {{"three_d", "time", "koopman_spectrum"}, true}},
        {"rossler", {{"three_d", "time", "koopman_spectrum"}, true}},
        {"chua", {{"three_d", "time", "koopman_spectrum"}, true}},
        {"hindmarsh_rose", {{"three_d", "time"}, true}},
        {"hodgkin_huxley", {{"time", "phase"}}},
        {"fput", {{"time", "fourier"}}},
        {"kuramoto", {{"time", "fourier"}}},
        {"ising", {{"time"}}},
        {"ising_3d", {{"time"}}},
        {"potts", {{"time"}}},
        {"mackey_glass", {{"time", "phase"}}},
    };
    return defaults;
}

const Defaults* find(std::string system) {
    std::transform(system.begin(), system.end(), system.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    for (auto& c : system) if (c == '-' && system != "rayleigh-plesset") c = '_';
    const auto found = table().find(system);
    return found == table().end() ? nullptr : &found->second;
}

}  // namespace

void orderViewsForSystem(const std::string& system, std::vector<std::string>& views) {
    const auto* defaults = find(system);
    if (!defaults) return;
    std::vector<std::string> ordered;
    for (const auto& view : defaults->views)
        if (std::find(views.begin(), views.end(), view) != views.end()) ordered.push_back(view);
    for (const auto& view : views)
        if (std::find(ordered.begin(), ordered.end(), view) == ordered.end()) ordered.push_back(view);
    views = std::move(ordered);
}

std::optional<std::pair<double, double>> systemCamera(const std::string& system) {
    const auto* defaults = find(system);
    // Eye direction (1.55, -1.65, 0.82) about the box centre: azimuth -46.8 deg, elevation 19.9 deg.
    if (!defaults || !defaults->oblique) return std::nullopt;
    return std::make_pair(-46.8, 19.9);
}

}  // namespace dsl
