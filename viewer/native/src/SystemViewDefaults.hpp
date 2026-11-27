#pragma once
// Per-system presentation defaults: which views a trajectory archive opens with, and the initial
// camera of its three-state view. Presentation only; never changes archived data.
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace dsl {

// Moves the system's recommended views (in order) to the front of `views`; others keep their order.
void orderViewsForSystem(const std::string& system, std::vector<std::string>& views);
// Initial (azimuth, elevation) in degrees for the three-state view, when the system has one.
[[nodiscard]] std::optional<std::pair<double, double>> systemCamera(const std::string& system);

}  // namespace dsl
