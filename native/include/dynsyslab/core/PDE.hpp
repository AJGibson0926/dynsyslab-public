#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class PDEIntegrationMethod { euler, rk4 };
enum class PDEBoundary { edge, periodic, zeroValue };

struct StructuredMesh {
    std::vector<std::size_t> shape;
    std::vector<double> spacing;
    std::vector<double> origin;
    std::vector<PDEBoundary> boundary;

    [[nodiscard]] std::size_t cellCount() const;
    void validate() const;
};

struct PDEIntegrationOptions {
    double initialTime{};
    double finalTime{1.0};
    double outputStep{0.01};
    double internalStep{0.01};
    PDEIntegrationMethod method{PDEIntegrationMethod::rk4};
};

using PDERightHandSide = std::function<void(
    double time, const StructuredMesh& mesh, const double* field,
    std::size_t valueCount, double* derivative)>;

struct PDEIntegrationResult {
    std::vector<double> time;
    std::vector<double> values;
    std::vector<std::size_t> fieldShape;
    std::size_t sampleCount{};
    std::size_t fieldValueCount{};
    std::size_t acceptedSteps{};
};

[[nodiscard]] PDEIntegrationResult integratePDE(
    const StructuredMesh& mesh, const std::vector<std::size_t>& fieldShape,
    const std::vector<double>& initialField, const PDEIntegrationOptions& options,
    const PDERightHandSide& rhs,
    const ExecutionContext* execution = nullptr);

[[nodiscard]] std::vector<AnalyticSolutionFamily> heatEquationAnalyticSolutions(
    const SystemConfiguration& configuration);
[[nodiscard]] DynamicalKnowledge heatEquationKnowledge();
[[nodiscard]] SystemDescriptor heatEquationDescriptor(
    std::size_t points = 64);
[[nodiscard]] SystemConfiguration resolveHeatEquation(
    const SystemConfiguration& requested);
[[nodiscard]] StructuredMesh heatEquationMesh(
    const SystemConfiguration& configuration);
[[nodiscard]] std::vector<double> heatEquationCoordinates(
    const SystemConfiguration& configuration);
[[nodiscard]] std::vector<double> heatEquationInitialField(
    const SystemConfiguration& configuration);
[[nodiscard]] PDERightHandSide heatEquationRightHandSide(
    const SystemConfiguration& configuration);
[[nodiscard]] double heatEquationStableStep(
    const SystemConfiguration& configuration);

}  // namespace dynsyslab::core
