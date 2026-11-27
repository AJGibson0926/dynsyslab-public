#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Called every `interval` sweeps (before the last) with the sweeps done, the observations so far
// (sweeps x 2 row-major: energy, magnetisation) and the current spins: durable mid-run history.
struct LatticeCheckpoint {
    std::size_t interval{0};
    std::function<void(std::size_t sweeps, const std::vector<double>& observations,
                       const std::vector<std::int64_t>& spins)> callback;
};

struct IsingRequest {
    std::vector<std::int64_t> spins;
    std::size_t size{};
    std::size_t dimensions{2};
    double coupling{1.0};
    double field{0.0};
    double temperature{1.0};
    std::uint64_t seed{};
    std::size_t sweeps{};
    LatticeCheckpoint checkpoint;
};

struct IsingResult {
    std::vector<double> observations;
    std::vector<double> initialObservation;
    std::vector<std::int64_t> initialSpins;
    std::vector<std::int64_t> spins;
    std::size_t sweeps{};
    std::size_t acceptedMoves{};
};

struct PottsRequest {
    std::vector<std::int64_t> spins;
    std::size_t size{};
    std::size_t states{3};
    double coupling{1.0};
    double field{0.0};
    double temperature{1.0};
    std::uint64_t seed{};
    std::size_t sweeps{};
    LatticeCheckpoint checkpoint;
};

struct PottsResult {
    std::vector<double> observations;
    std::vector<double> initialObservation;
    std::vector<std::int64_t> initialSpins;
    std::vector<std::int64_t> spins;
    std::size_t sweeps{};
    std::size_t acceptedMoves{};
};

[[nodiscard]] IsingResult runIsing2D(
    const IsingRequest& request, const ExecutionContext* execution = nullptr);
[[nodiscard]] IsingResult runIsing3D(
    const IsingRequest& request, const ExecutionContext* execution = nullptr);
[[nodiscard]] PottsResult runPotts2D(
    const PottsRequest& request, const ExecutionContext* execution = nullptr);

[[nodiscard]] const SystemDescriptor& ising2DDescriptor();
[[nodiscard]] const SystemDescriptor& ising3DDescriptor();
[[nodiscard]] const SystemDescriptor& potts2DDescriptor();
[[nodiscard]] bool isLatticeSystemIdentifier(const std::string& identifier);
[[nodiscard]] SystemConfiguration resolveLatticeSystem(
    const SystemConfiguration& requested);

}  // namespace dynsyslab::core
