#pragma once
#include <cstddef>
#include <optional>
#include <memory>
#include <string>
#include <vector>
namespace dynsyslab::core {
struct PooledHankelBasis {
    std::vector<double> leftSingularVectors; // delay x rank, row-major
    std::vector<double> singularValues;
    std::vector<double> cumulativeEnergy;
    double observationScale{1.0};
    double totalHankelEnergy{};
    std::size_t delays{};
    std::size_t rank{};
    std::size_t trainingTrajectories{};
    std::size_t trainingWindows{};
};

// Immutable fitted transform. h is ordered oldest-to-newest and the retained
// coordinates are (U^T h)/sigma, preserving the established native convention.
// observationScale records fitting normalization; it is NOT applied a second time.
class FrozenHistoryTransform {
public:
    FrozenHistoryTransform(PooledHankelBasis basis, double sampleInterval,
        std::vector<std::size_t> coordinates = {}, std::string measurement = {},
        std::string timeUnits = "model-time", std::string observationUnits = "unspecified");
    const PooledHankelBasis& basis() const noexcept { return *basis_; }
    double sampleInterval() const noexcept { return sampleInterval_; }
    const std::vector<std::size_t>& retainedCoordinates() const noexcept { return coordinates_; }
    const std::string& measurement() const noexcept { return measurement_; }
    const std::string& timeUnits() const noexcept { return timeUnits_; }
    const std::string& observationUnits() const noexcept { return observationUnits_; }
    const std::string& identity() const noexcept { return identity_; }
private:
    std::shared_ptr<const PooledHankelBasis> basis_;
    double sampleInterval_{};
    std::vector<std::size_t> coordinates_;
    std::string measurement_,timeUnits_,observationUnits_,identity_;
};

struct HistoryProjection {
    std::vector<double> coordinates; // window x rank, row-major
    std::size_t windows{};
    std::size_t rank{};
};

// Projects a complete observation record through a previously fitted basis.
// The coordinate at row j belongs to the history ending at sample
// j + basis.delays - 1.
[[nodiscard]] HistoryProjection projectHistory(
    const std::vector<double>& observations,
    const PooledHankelBasis& basis);

// Online counterpart to projectHistory. The fitted basis is immutable; only
// the causal observation ring buffer changes while measurements arrive.
class CausalHistoryProjector {
public:
    explicit CausalHistoryProjector(PooledHankelBasis basis);
    explicit CausalHistoryProjector(FrozenHistoryTransform transform);
    [[nodiscard]] std::optional<std::vector<double>> push(double observation);
    [[nodiscard]] std::optional<std::vector<double>> push(double time, double observation);
    // Multi-rate acquisition: append every sensor tick, project only when needed.
    void append(double observation);
    void append(double time, double observation);
    [[nodiscard]] std::optional<std::vector<double>> projectLatest() const;
    void reset() noexcept;
    [[nodiscard]] std::size_t samplesSeen() const noexcept { return samplesSeen_; }
    [[nodiscard]] const PooledHankelBasis& basis() const noexcept { return transform_.basis(); }
    [[nodiscard]] const FrozenHistoryTransform& transform() const noexcept { return transform_; }

private:
    void appendValue(double observation);
    FrozenHistoryTransform transform_;
    std::vector<double> ring_;
    std::size_t writeIndex_{};
    std::size_t samplesSeen_{};
    std::optional<double> lastTime_;
};

} // namespace dynsyslab::core
