#include <dynsyslab/core/Lattice.hpp>
#include "LatticeResearch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::size_t checkedPower(std::size_t size, std::size_t dimensions) {
    if (size < 2 || (dimensions != 2 && dimensions != 3))
        throw std::invalid_argument("lattice size must be at least two in two or three dimensions");
    std::size_t count = 1;
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
        if (count > (std::numeric_limits<std::size_t>::max)() / size)
            throw ExecutionLimitExceeded("lattice site count is not addressable");
        count *= size;
    }
    return count;
}

void validateThermodynamics(double coupling, double field, double temperature) {
    if (!std::isfinite(coupling) || !std::isfinite(field)
        || !std::isfinite(temperature) || temperature <= 0.0)
        throw std::invalid_argument(
            "lattice coupling/field must be finite and temperature positive");
}

void requireBuffers(std::size_t sites, std::size_t sweeps,
                    const ExecutionContext* execution) {
    if (sweeps > (std::numeric_limits<std::size_t>::max)() / 2)
        throw ExecutionLimitExceeded("lattice observation count is not addressable");
    if (!execution) return;
    const long double bytes = static_cast<long double>(sites) * 3.0L
        * sizeof(std::int64_t) + static_cast<long double>(sweeps) * 2.0L
        * sizeof(double);
    if (bytes > static_cast<long double>((std::numeric_limits<std::size_t>::max)()))
        throw ExecutionLimitExceeded("lattice buffers are not addressable");
    execution->requireMemory(static_cast<std::size_t>(bytes), "lattice buffers");
}

void report(const ExecutionContext* execution, const char* model,
            std::size_t sweep, std::size_t sweeps) {
    if (!execution || sweeps == 0) return;
    const std::size_t stride = std::max<std::size_t>(1, sweeps / 100);
    if ((sweep + 1) % stride == 0 || sweep + 1 == sweeps)
        execution->report({"lattice", std::string("Advanced ") + model,
            static_cast<double>(sweep + 1) / sweeps, sweep + 1, sweeps});
}

std::vector<SystemDescriptor> latticeCatalog() {
    return {ising2DDescriptor(), ising3DDescriptor(), potts2DDescriptor()};
}

const bool registered = [] {
    registerSystemCatalogProvider(&latticeCatalog);
    return true;
}();

double requestedParameter(const SystemConfiguration& configuration,
                          const char* name, double fallback) {
    const auto found = configuration.parameters.find(name);
    return found == configuration.parameters.end() ? fallback : found->second;
}

std::array<double, 2> ising2DObservation(const std::vector<std::int64_t>& spins,
                                         std::size_t size, double coupling,
                                         double field) {
    double magnetization = 0.0, energy = 0.0;
    const auto at = [&spins, size](std::size_t row, std::size_t column) {
        return spins[row * size + column];
    };
    for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column < size; ++column) {
            const double spin = static_cast<double>(at(row, column));
            magnetization += spin;
            energy -= coupling * spin * static_cast<double>(
                at((row + 1) % size, column) + at(row, (column + 1) % size));
            energy -= field * spin;
        }
    const double count = static_cast<double>(size * size);
    return {magnetization / count, energy / count};
}

std::array<double, 2> ising3DObservation(const std::vector<std::int64_t>& spins,
                                         std::size_t size, double coupling,
                                         double field) {
    const auto at = [&spins, size](std::size_t x, std::size_t y, std::size_t z) {
        return spins[(x * size + y) * size + z];
    };
    double magnetization = 0.0, energy = 0.0;
    for (std::size_t x = 0; x < size; ++x)
        for (std::size_t y = 0; y < size; ++y)
            for (std::size_t z = 0; z < size; ++z) {
                const double spin = static_cast<double>(at(x, y, z));
                magnetization += spin;
                energy -= coupling * spin * static_cast<double>(
                    at((x + 1) % size, y, z) + at(x, (y + 1) % size, z)
                    + at(x, y, (z + 1) % size));
                energy -= field * spin;
            }
    const double count = static_cast<double>(size * size * size);
    return {magnetization / count, energy / count};
}

std::array<double, 2> potts2DObservation(const std::vector<std::int64_t>& spins,
                                         std::size_t size, std::size_t states,
                                         double coupling, double field) {
    const auto at = [&spins, size](std::size_t row, std::size_t column) {
        return spins[row * size + column];
    };
    std::vector<std::size_t> populations(states, 0);
    double energy = 0.0;
    for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column < size; ++column) {
            const auto spin = at(row, column);
            ++populations[static_cast<std::size_t>(spin)];
            energy -= coupling * static_cast<double>(spin == at((row + 1) % size, column));
            energy -= coupling * static_cast<double>(spin == at(row, (column + 1) % size));
            energy -= field * static_cast<double>(spin == 0);
        }
    const auto largest = *std::max_element(populations.begin(), populations.end());
    const double count = static_cast<double>(size * size);
    return {(states * static_cast<double>(largest) / count - 1.0)
        / static_cast<double>(states - 1), energy / count};
}

}  // namespace

IsingResult runIsing2D(const IsingRequest& request,
                       const ExecutionContext* execution) {
    if (request.dimensions != 2)
        throw std::invalid_argument("2-D Ising request has the wrong dimensionality");
    const std::size_t count = checkedPower(request.size, 2);
    validateThermodynamics(request.coupling, request.field, request.temperature);
    requireBuffers(count, request.sweeps, execution);
    std::uniform_int_distribution<std::size_t> site(0, request.size - 1);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::mt19937_64 generator(request.seed);
    IsingResult result;
    result.spins = request.spins;
    if (result.spins.empty()) {
        result.spins.resize(count);
        std::uniform_int_distribution<int> initial(0, 1);
        for (auto& spin : result.spins) spin = initial(generator) ? 1 : -1;
    }
    if (result.spins.size() != count)
        throw std::invalid_argument("2-D Ising spin count does not match lattice size");
    for (const auto spin : result.spins)
        if (spin != -1 && spin != 1)
            throw std::invalid_argument("Ising spins must contain only -1 and 1");
    result.initialSpins = result.spins;
    const auto initialObservation = ising2DObservation(
        result.spins, request.size, request.coupling, request.field);
    result.initialObservation.assign(initialObservation.begin(), initialObservation.end());
    result.observations.resize(request.sweeps * 2);
    result.sweeps = request.sweeps;
    auto at = [&result, &request](std::size_t row, std::size_t column) -> std::int64_t& {
        return result.spins[row * request.size + column];
    };
    for (std::size_t sweep = 0; sweep < request.sweeps; ++sweep) {
        if (execution) execution->checkpoint();
        for (std::size_t attempt = 0; attempt < count; ++attempt) {
            const auto row = site(generator);
            const auto column = site(generator);
            const auto spin = at(row, column);
            const auto neighbours = at((row + request.size - 1) % request.size, column) +
                                    at((row + 1) % request.size, column) +
                                    at(row, (column + request.size - 1) % request.size) +
                                    at(row, (column + 1) % request.size);
            const double delta = 2.0 * static_cast<double>(spin) *
                                 (request.coupling * static_cast<double>(neighbours) + request.field);
            if (delta <= 0.0 || uniform(generator) < std::exp(-delta / request.temperature)) {
                at(row, column) = -spin;
                ++result.acceptedMoves;
            }
        }
        const auto observation = ising2DObservation(
            result.spins, request.size, request.coupling, request.field);
        result.observations[sweep * 2] = observation[0];
        result.observations[sweep * 2 + 1] = observation[1];
        report(execution, "2-D Ising lattice", sweep, request.sweeps);
    }
    return result;
}

IsingResult runIsing3D(const IsingRequest& request,
                       const ExecutionContext* execution) {
    if (request.dimensions != 3)
        throw std::invalid_argument("3-D Ising request has the wrong dimensionality");
    const std::size_t count = checkedPower(request.size, 3);
    validateThermodynamics(request.coupling, request.field, request.temperature);
    requireBuffers(count, request.sweeps, execution);
    std::mt19937_64 generator(request.seed);
    std::uniform_int_distribution<std::size_t> site(0, request.size - 1);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    IsingResult result;
    result.spins = request.spins;
    if (result.spins.empty()) {
        result.spins.resize(count);
        std::uniform_int_distribution<int> initial(0, 1);
        for (auto& spin : result.spins) spin = initial(generator) ? 1 : -1;
    }
    if (result.spins.size() != count)
        throw std::invalid_argument("3-D Ising spin count does not match lattice size");
    for (const auto spin : result.spins)
        if (spin != -1 && spin != 1)
            throw std::invalid_argument("Ising spins must contain only -1 and 1");
    result.initialSpins = result.spins;
    const auto initialObservation = ising3DObservation(
        result.spins, request.size, request.coupling, request.field);
    result.initialObservation.assign(initialObservation.begin(), initialObservation.end());
    result.observations.resize(request.sweeps * 2);
    result.sweeps = request.sweeps;
    const auto index = [&request](std::size_t x, std::size_t y, std::size_t z) {
        return (x * request.size + y) * request.size + z;
    };
    for (std::size_t sweep = 0; sweep < request.sweeps; ++sweep) {
        if (execution) execution->checkpoint();
        for (std::size_t attempt = 0; attempt < count; ++attempt) {
            const auto x = site(generator), y = site(generator), z = site(generator);
            const auto active = index(x, y, z);
            const auto spin = result.spins[active];
            const auto neighbours =
                result.spins[index((x + request.size - 1) % request.size, y, z)]
                + result.spins[index((x + 1) % request.size, y, z)]
                + result.spins[index(x, (y + request.size - 1) % request.size, z)]
                + result.spins[index(x, (y + 1) % request.size, z)]
                + result.spins[index(x, y, (z + request.size - 1) % request.size)]
                + result.spins[index(x, y, (z + 1) % request.size)];
            const double delta = 2.0 * static_cast<double>(spin)
                * (request.coupling * static_cast<double>(neighbours) + request.field);
            if (delta <= 0.0 || uniform(generator) < std::exp(-delta / request.temperature)) {
                result.spins[active] = -spin;
                ++result.acceptedMoves;
            }
        }
        const auto observation = ising3DObservation(
            result.spins, request.size, request.coupling, request.field);
        result.observations[sweep * 2] = observation[0];
        result.observations[sweep * 2 + 1] = observation[1];
        report(execution, "3-D Ising lattice", sweep, request.sweeps);
    }
    return result;
}

PottsResult runPotts2D(const PottsRequest& request,
                       const ExecutionContext* execution) {
    const std::size_t count = checkedPower(request.size, 2);
    if (request.states < 2
        || request.states > static_cast<std::size_t>((std::numeric_limits<std::int64_t>::max)()))
        throw std::invalid_argument("Potts model needs an addressable state count of at least two");
    validateThermodynamics(request.coupling, request.field, request.temperature);
    requireBuffers(count, request.sweeps, execution);
    std::mt19937_64 generator(request.seed);
    std::uniform_int_distribution<std::size_t> site(0, request.size - 1);
    std::uniform_int_distribution<std::int64_t> initial(
        0, static_cast<std::int64_t>(request.states - 1));
    std::uniform_int_distribution<std::int64_t> proposal(
        0, static_cast<std::int64_t>(request.states - 2));
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    PottsResult result;
    result.spins = request.spins;
    if (result.spins.empty()) {
        result.spins.resize(count);
        for (auto& spin : result.spins) spin = initial(generator);
    }
    if (result.spins.size() != count)
        throw std::invalid_argument("Potts spin count does not match lattice size");
    for (const auto spin : result.spins)
        if (spin < 0 || spin >= static_cast<std::int64_t>(request.states))
            throw std::invalid_argument("Potts spins are outside the configured state range");
    result.initialSpins = result.spins;
    const auto initialObservation = potts2DObservation(
        result.spins, request.size, request.states,
        request.coupling, request.field);
    result.initialObservation.assign(initialObservation.begin(), initialObservation.end());
    result.observations.resize(request.sweeps * 2);
    result.sweeps = request.sweeps;
    const auto index = [&request](std::size_t row, std::size_t column) {
        return row * request.size + column;
    };
    for (std::size_t sweep = 0; sweep < request.sweeps; ++sweep) {
        if (execution) execution->checkpoint();
        for (std::size_t attempt = 0; attempt < count; ++attempt) {
            const auto row = site(generator), column = site(generator);
            const auto active = index(row, column);
            const auto old = result.spins[active];
            auto next = proposal(generator);
            if (next >= old) ++next;
            const std::int64_t neighbours[] = {
                result.spins[index((row + request.size - 1) % request.size, column)],
                result.spins[index((row + 1) % request.size, column)],
                result.spins[index(row, (column + request.size - 1) % request.size)],
                result.spins[index(row, (column + 1) % request.size)]};
            std::size_t oldMatches = 0, newMatches = 0;
            for (const auto neighbour : neighbours) {
                oldMatches += neighbour == old;
                newMatches += neighbour == next;
            }
            const double fieldChange = static_cast<double>(next == 0)
                - static_cast<double>(old == 0);
            const double delta = -request.coupling
                * (static_cast<double>(newMatches) - static_cast<double>(oldMatches))
                - request.field * fieldChange;
            if (delta <= 0.0 || uniform(generator) < std::exp(-delta / request.temperature)) {
                result.spins[active] = next;
                ++result.acceptedMoves;
            }
        }
        const auto observation = potts2DObservation(
            result.spins, request.size, request.states,
            request.coupling, request.field);
        result.observations[sweep * 2] = observation[0];
        result.observations[sweep * 2 + 1] = observation[1];
        report(execution, "2-D Potts lattice", sweep, request.sweeps);
    }
    return result;
}

const SystemDescriptor& ising2DDescriptor() {
    static const SystemDescriptor value = [] {
      SystemDescriptor result{"ising", "2-D Ising model", 2, "discrete",
        {"magnetization", "energy_density"}, {"magnetization", "energy density"},
        {{"critical", "Square lattice near the critical temperature"}}, "lattice"};
      result.research = detail::latticeResearch("ising");
      return result;
    }();
    return value;
}

const SystemDescriptor& ising3DDescriptor() {
    static const SystemDescriptor value = [] {
      SystemDescriptor result{"ising_3d", "3-D Ising model", 2, "discrete",
        {"magnetization", "energy_density"}, {"magnetization", "energy density"},
        {{"critical", "Cubic lattice near the critical temperature"}}, "lattice"};
      result.research = detail::latticeResearch("ising_3d");
      return result;
    }();
    return value;
}

const SystemDescriptor& potts2DDescriptor() {
    static const SystemDescriptor value = [] {
      SystemDescriptor result{"potts", "2-D Potts model", 2, "discrete",
        {"order_parameter", "energy_density"}, {"order parameter", "energy density"},
        {{"q3_reference", "Three-state square-lattice reference"}}, "lattice"};
      result.research = detail::latticeResearch("potts");
      return result;
    }();
    return value;
}

bool isLatticeSystemIdentifier(const std::string& identifier) {
    return identifier == "ising" || identifier == "ising_2d"
        || identifier == "ising-2d" || identifier == "ising_3d"
        || identifier == "ising-3d" || identifier == "potts"
        || identifier == "potts_2d" || identifier == "potts-2d";
}

SystemConfiguration resolveLatticeSystem(const SystemConfiguration& requested) {
    static_cast<void>(registered);
    SystemConfiguration result = requested;
    const std::string id = requested.identifier;
    if (!isLatticeSystemIdentifier(id))
        throw std::invalid_argument("unknown native lattice system: " + id);
    if (id == "ising" || id == "ising_2d" || id == "ising-2d") {
        result.identifier = "ising";
        if (result.preset.empty()) result.preset = "critical";
        result.parameters["size"] = requestedParameter(result, "size", 32.0);
        result.parameters["coupling"] = requestedParameter(result, "coupling", 1.0);
        result.parameters["field"] = requestedParameter(result, "field", 0.0);
        result.parameters["temperature"] = requestedParameter(result, "temperature", 2.269);
        result.equationsVersion = "dynsyslab.ising_2d/1";
    } else if (id == "ising_3d" || id == "ising-3d") {
        result.identifier = "ising_3d";
        if (result.preset.empty()) result.preset = "critical";
        result.parameters["size"] = requestedParameter(result, "size", 16.0);
        result.parameters["coupling"] = requestedParameter(result, "coupling", 1.0);
        result.parameters["field"] = requestedParameter(result, "field", 0.0);
        result.parameters["temperature"] = requestedParameter(result, "temperature", 4.5115);
        result.equationsVersion = "dynsyslab.ising_3d/1";
    } else {
        result.identifier = "potts";
        if (result.preset.empty()) result.preset = "q3_reference";
        result.parameters["size"] = requestedParameter(result, "size", 32.0);
        result.parameters["states"] = requestedParameter(result, "states", 3.0);
        result.parameters["coupling"] = requestedParameter(result, "coupling", 1.0);
        result.parameters["field"] = requestedParameter(result, "field", 0.0);
        result.parameters["temperature"] = requestedParameter(result, "temperature", 1.0);
        result.equationsVersion = "dynsyslab.potts_2d/1";
    }
    result.parameterUnits = {{"size", "sites/edge"}, {"coupling", "energy"},
        {"field", "energy"}, {"temperature", "energy"}};
    if (result.identifier == "potts") result.parameterUnits["states"] = "count";
    result.initialState = {0.0, 0.0};
    return result;
}

}  // namespace dynsyslab::core
