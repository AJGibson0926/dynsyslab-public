#include "ArchiveSystem.hpp"
#include "DisplayCoordinates.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
namespace dsl {
dynsyslab::core::SystemConfiguration archiveSystemConfiguration(const json::Value &metadata) {
    dynsyslab::core::SystemConfiguration result;
    result.identifier =
        metadata.at("system").string(metadata.at("system_configuration").at("identifier").string());
    const auto folded = [](const std::string &value) {
        std::string token;
        for (const auto c : value)
            if (std::isalnum(static_cast<unsigned char>(c)))
                token.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        return token;
    };
    for (const auto &descriptor : dynsyslab::core::availableSystems())
        if (folded(result.identifier) == folded(descriptor.identifier) ||
            folded(result.identifier) == folded(descriptor.label)) {
            result.identifier = descriptor.identifier;
            break;
        }
    const auto &legacy = metadata.at("system_parameters");
    result.preset =
        metadata.at("system_preset")
            .string(
                metadata.at("system_configuration")
                    .at("preset")
                    .string(legacy.at("preset").at("key").string(legacy.at("preset").string())));
    for (const auto *parameters :
         {&legacy.at("preset"), &legacy, &metadata.at("system_configuration").at("parameters"),
          &metadata.at("resolved_system_parameters")})
        for (const auto &entry : parameters->object())
            if (std::holds_alternative<double>(entry.second.storage()))
                result.parameters[entry.first] = entry.second.number();
    const auto &driver = !metadata.at("system_driver").isNull() ? metadata.at("system_driver")
                         : metadata.at("system_configuration").at("driver");
    result.driver = {driver.at("amplitude").number(), driver.at("frequency_hz").number(),
                     driver.at("phase_radians").number()};
    if (result.identifier.empty())
        throw std::invalid_argument("Viewer requires system metadata");
    return result;
}

void appendVerifiedEquilibria(Figure &figure, const json::Value &metadata) {
    if (figure.kind != "phase")
        return;
    try {
        auto system = dynsyslab::core::createSystem(archiveSystemConfiguration(metadata));
        const auto n = system->descriptor().stateDimension;
        if (n != 2)
            return;
        const auto fingerprint = metadata.at("equations_fingerprint").string();
        if (!fingerprint.empty() &&
            fingerprint != dynsyslab::core::equationsFingerprint(system->configuration()))
            return;
        const auto knowledge = system->dynamicalKnowledge();
        if (!knowledge.equilibriumGuides)
            return;
        const auto native =
            displayCoordinateTransform(metadata, n, system->descriptor().preferredCoordinates);
        const auto display = displayCoordinateTransform(metadata, n, figure.coordinates);
        if (native.coordinates != system->descriptor().preferredCoordinates)
            return;
        const double sx = display.stateMultipliers[0] / native.stateMultipliers[0];
        const double sy = display.stateMultipliers[1] / native.stateMultipliers[1];
        const auto bounds = automaticAxisLimits(figure);
        const std::vector<double> low{bounds.xMin / sx, bounds.yMin / sy},
            high{bounds.xMax / sx, bounds.yMax / sy};
        for (const auto &guide : knowledge.equilibriumGuides(low, high)) {
            if (guide.state.size() != n || system->invalid(0, guide.state.data()))
                continue;
            std::vector<double> f(n);
            system->rhsWithInput(0, guide.state.data(), 0, f.data());
            if (std::any_of(f.begin(), f.end(),
                            [](double v) { return !std::isfinite(v) || std::abs(v) > 1e-8; }))
                continue;
            const double x = guide.state[0] * sx, y = guide.state[1] * sy;
            if (x < bounds.xMin || x > bounds.xMax || y < bounds.yMin || y > bounds.yMax)
                continue;
            if (std::none_of(figure.references.begin(), figure.references.end(),
                             [&](const auto &p) { return p.equilibrium && p.x == x && p.y == y; }))
                figure.references.push_back({x, y, "EC", true});
        }
    } catch (const std::exception
                 &) { /* Historical archives without a resolvable model remain viewable. */
    }
}
} // namespace dsl
