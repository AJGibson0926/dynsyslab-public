#include <dynsyslab/core/PDE.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

std::size_t checkedProduct(const std::vector<std::size_t>& shape) {
    std::size_t result = 1;
    for (const auto extent : shape) {
        if (!extent || result > (std::numeric_limits<std::size_t>::max)() / extent)
            throw std::invalid_argument("PDE field shape is empty or too large");
        result *= extent;
    }
    return result;
}

double parameter(const SystemConfiguration& configuration,
                 const char* name, double fallback) {
    const auto found = configuration.parameters.find(name);
    return found == configuration.parameters.end() ? fallback : found->second;
}

std::size_t pointCount(const SystemConfiguration& configuration) {
    const double value = parameter(configuration, "points", 64.0);
    if (!std::isfinite(value) || value < 3.0 || std::floor(value) != value
        || value > static_cast<double>((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument("heat-equation points must be an integer of at least 3");
    return static_cast<std::size_t>(value);
}

int boundaryCode(const SystemConfiguration& configuration) {
    const double value = parameter(configuration, "boundary_code", 0.0);
    if (!std::isfinite(value) || std::floor(value) != value || value < 0.0 || value > 2.0)
        throw std::invalid_argument(
            "heat-equation boundary_code must be 0 (periodic), 1 (zero flux), or 2 (zero value)");
    return static_cast<int>(value);
}

std::vector<SystemDescriptor> pdeCatalog() {
    return {heatEquationDescriptor()};
}

struct CatalogRegistration {
    CatalogRegistration() { registerSystemCatalogProvider(&pdeCatalog); }
} registration;

void requireFinite(const std::vector<double>& values) {
    if (std::any_of(values.begin(), values.end(),
                    [](double value) { return !std::isfinite(value); }))
        throw std::runtime_error("PDE integration produced a non-finite field");
}

}  // namespace

std::size_t StructuredMesh::cellCount() const {
    return checkedProduct(shape);
}

void StructuredMesh::validate() const {
    if (shape.empty() || spacing.size() != shape.size())
        throw std::invalid_argument("structured PDE mesh shape and spacing must agree");
    if (!origin.empty() && origin.size() != shape.size())
        throw std::invalid_argument("structured PDE mesh origin must match its dimension");
    if (!boundary.empty() && boundary.size() != shape.size())
        throw std::invalid_argument("structured PDE mesh boundary must match its dimension");
    for (std::size_t axis = 0; axis < shape.size(); ++axis)
        if (shape[axis] < 2 || !std::isfinite(spacing[axis]) || spacing[axis] <= 0.0)
            throw std::invalid_argument("structured PDE mesh axes require at least two cells and positive spacing");
    (void)cellCount();
}

PDEIntegrationResult integratePDE(
    const StructuredMesh& mesh, const std::vector<std::size_t>& fieldShape,
    const std::vector<double>& initialField, const PDEIntegrationOptions& options,
    const PDERightHandSide& rhs, const ExecutionContext* execution) {
    mesh.validate();
    if (!rhs) throw std::invalid_argument("PDE right-hand side is empty");
    if (options.method == PDEIntegrationMethod::adaptive)
        return integratePDEAdaptive(mesh, fieldShape, initialField, options, rhs, execution);
    if (!std::isfinite(options.initialTime) || !std::isfinite(options.finalTime)
        || !std::isfinite(options.outputStep) || !std::isfinite(options.internalStep)
        || !(options.initialTime < options.finalTime) || options.outputStep <= 0.0
        || options.internalStep <= 0.0)
        throw std::invalid_argument("PDE time span and steps are invalid");
    const std::size_t values = checkedProduct(fieldShape);
    if (fieldShape.size() < mesh.shape.size()
        || !std::equal(mesh.shape.begin(), mesh.shape.end(), fieldShape.begin())
        || values != initialField.size() || values % mesh.cellCount() != 0)
        throw std::invalid_argument("PDE initial field does not match its mesh/field shape");
    const double duration = options.finalTime - options.initialTime;
    const double outputIntervals = std::ceil(duration / options.outputStep);
    if (!std::isfinite(outputIntervals) || outputIntervals < 1.0
        || outputIntervals >= static_cast<double>((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument("PDE output grid is too large");
    const std::size_t outputs = static_cast<std::size_t>(outputIntervals) + 1;
    if (outputs > (std::numeric_limits<std::size_t>::max)() / values)
        throw std::invalid_argument("PDE solution allocation is too large");
    const std::size_t stored = outputs * values;
    if (stored > (std::numeric_limits<std::size_t>::max)() / sizeof(double)
        || values > (std::numeric_limits<std::size_t>::max)() / (6 * sizeof(double)))
        throw std::invalid_argument("PDE solution allocation is too large");
    if (execution) {
        const std::size_t bytes = stored * sizeof(double);
        const std::size_t workBytes = 6 * values * sizeof(double);
        if (bytes > (std::numeric_limits<std::size_t>::max)() - workBytes)
            throw std::invalid_argument("PDE solution allocation is too large");
        execution->requireMemory(bytes + workBytes,
                                 "PDE solution and RK work arrays");
        execution->report({"pde", "Preparing PDE integration", 0.0, 0, outputs});
    }
    PDEIntegrationResult result;
    result.fieldShape = fieldShape;
    result.fieldValueCount = values;
    result.time.reserve(outputs);
    result.values.reserve(stored);
    std::vector<double> field = initialField;
    std::vector<double> k1(values), k2(values), k3(values), k4(values), work(values);
    result.time.push_back(options.initialTime);
    result.values.insert(result.values.end(), field.begin(), field.end());
    double current = options.initialTime;
    std::size_t output = 1;
    int reportedPercent = 0;
    while (current < options.finalTime) {
        const double target = std::min(
            options.finalTime, options.initialTime + output * options.outputStep);
        while (current < target) {
            if (execution) execution->checkpoint();
            const double step = std::min(options.internalStep, target - current);
            if (!(step > 0.0) || current + step == current)
                throw std::runtime_error("PDE step is below floating-point time resolution");
            if (options.method == PDEIntegrationMethod::euler) {
                rhs(current, mesh, field.data(), values, k1.data());
                for (std::size_t index = 0; index < values; ++index)
                    field[index] += step * k1[index];
            } else {
                rhs(current, mesh, field.data(), values, k1.data());
                for (std::size_t i = 0; i < values; ++i) work[i] = field[i] + 0.5 * step * k1[i];
                rhs(current + 0.5 * step, mesh, work.data(), values, k2.data());
                for (std::size_t i = 0; i < values; ++i) work[i] = field[i] + 0.5 * step * k2[i];
                rhs(current + 0.5 * step, mesh, work.data(), values, k3.data());
                for (std::size_t i = 0; i < values; ++i) work[i] = field[i] + step * k3[i];
                rhs(current + step, mesh, work.data(), values, k4.data());
                for (std::size_t i = 0; i < values; ++i)
                    field[i] += step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
            }
            current += step;
            ++result.acceptedSteps;
            requireFinite(field);
        }
        current = target;
        result.time.push_back(current);
        result.values.insert(result.values.end(), field.begin(), field.end());
        const double fraction = (current - options.initialTime) / duration;
        const int percent = static_cast<int>(100.0 * fraction);
        if (execution && (percent > reportedPercent || target == options.finalTime)) {
            reportedPercent = percent;
            execution->report({"pde", "Integrating PDE field", fraction,
                output, outputs - 1});
        }
        ++output;
    }
    result.sampleCount = result.time.size();
    return result;
}

SystemDescriptor heatEquationDescriptor(std::size_t points) {
    SystemDescriptor result{"heat_equation_1d", "1-D heat equation", points,
        "dimensional", {}, {}, {{"periodic_diffusion", "Periodic diffusion"},
        {"periodic_decay", "Periodic diffusion with linear decay"}}, "pde"};
    result.research = heatEquationKnowledge();
    result.analyticSolutions = heatEquationAnalyticSolutions;
    result.stateNames.reserve(points);
    result.stateLabels.reserve(points);
    for (std::size_t index = 0; index < points; ++index) {
        result.stateNames.push_back("u_" + std::to_string(index));
        result.stateLabels.push_back("u(x_" + std::to_string(index) + ")");
    }
    return result;
}

SystemConfiguration resolveHeatEquation(const SystemConfiguration& requested) {
    if (!requested.identifier.empty() && requested.identifier != "heat_equation_1d")
        throw std::invalid_argument("unsupported native PDE system: " + requested.identifier);
    SystemConfiguration result = requested;
    result.identifier = "heat_equation_1d";
    if (result.preset.empty()) result.preset = "periodic_diffusion";
    std::map<std::string, double> defaults{{"points", 64.0}, {"length", 1.0},
        {"diffusivity", 0.01}, {"decay_rate", 0.0}, {"amplitude", 1.0},
        {"offset", 0.0}, {"wave_number", 1.0}, {"boundary_code", 0.0}};
    if (result.preset == "periodic_decay") defaults["decay_rate"] = 0.1;
    else if (result.preset != "periodic_diffusion")
        throw std::invalid_argument("unknown heat-equation preset: " + result.preset);
    for (const auto& item : result.parameters) defaults[item.first] = item.second;
    result.parameters = std::move(defaults);
    result.parameterUnits = {{"points", "grid points"}, {"length", "m"},
        {"diffusivity", "m^2/s"}, {"decay_rate", "1/s"},
        {"amplitude", "field units"}, {"offset", "field units"},
        {"wave_number", "dimensionless"}, {"boundary_code", "enum"}};
    result.equationsVersion = "dynsyslab.heat_equation_1d/1";
    const double length = parameter(result, "length", 1.0);
    const double diffusivity = parameter(result, "diffusivity", 0.01);
    const double decay = parameter(result, "decay_rate", 0.0);
    const double wave = parameter(result, "wave_number", 1.0);
    const double amplitude = parameter(result, "amplitude", 1.0);
    const double offset = parameter(result, "offset", 0.0);
    if (!std::isfinite(length) || length <= 0.0 || !std::isfinite(diffusivity)
        || diffusivity < 0.0 || !std::isfinite(decay) || decay < 0.0
        || !std::isfinite(wave) || std::floor(wave) != wave || wave < 0.0
        || !std::isfinite(amplitude) || !std::isfinite(offset))
        throw std::invalid_argument("heat-equation physical parameters are invalid");
    (void)pointCount(result);
    (void)boundaryCode(result);
    return result;
}

StructuredMesh heatEquationMesh(const SystemConfiguration& configuration) {
    const auto resolved = resolveHeatEquation(configuration);
    const std::size_t points = pointCount(resolved);
    const double length = parameter(resolved, "length", 1.0);
    const int boundary = boundaryCode(resolved);
    const PDEBoundary mode = boundary == 0 ? PDEBoundary::periodic
        : (boundary == 1 ? PDEBoundary::edge : PDEBoundary::zeroValue);
    return {{points}, {length / (boundary == 0 ? points : points - 1)},
            {0.0}, {mode}};
}

std::vector<double> heatEquationCoordinates(const SystemConfiguration& configuration) {
    const auto resolved = resolveHeatEquation(configuration);
    const auto mesh = heatEquationMesh(resolved);
    std::vector<double> result(mesh.shape[0]);
    for (std::size_t index = 0; index < result.size(); ++index)
        result[index] = mesh.origin[0] + index * mesh.spacing[0];
    return result;
}

std::vector<double> heatEquationInitialField(const SystemConfiguration& configuration) {
    const auto resolved = resolveHeatEquation(configuration);
    const auto coordinates = heatEquationCoordinates(resolved);
    const double length = parameter(resolved, "length", 1.0);
    const double amplitude = parameter(resolved, "amplitude", 1.0);
    const double offset = parameter(resolved, "offset", 0.0);
    const double wave = parameter(resolved, "wave_number", 1.0);
    std::vector<double> result(coordinates.size());
    for (std::size_t index = 0; index < result.size(); ++index)
        result[index] = offset + amplitude * std::sin(2.0 * kPi * wave * coordinates[index] / length);
    if (boundaryCode(resolved) == 2) result.front() = result.back() = 0.0;
    return result;
}

PDERightHandSide heatEquationRightHandSide(const SystemConfiguration& configuration) {
    const auto resolved = resolveHeatEquation(configuration);
    const double diffusivity = parameter(resolved, "diffusivity", 0.01);
    const double decay = parameter(resolved, "decay_rate", 0.0);
    const int boundary = boundaryCode(resolved);
    return [diffusivity, decay, boundary](double, const StructuredMesh& mesh,
        const double* field, std::size_t count, double* derivative) {
        if (mesh.shape.size() != 1 || count != mesh.shape[0])
            throw std::invalid_argument("1-D heat equation received an incompatible field");
        const std::size_t points = mesh.shape[0];
        const double scale = diffusivity / (mesh.spacing[0] * mesh.spacing[0]);
        for (std::size_t index = 0; index < points; ++index) {
            if (boundary == 2 && (index == 0 || index + 1 == points)) {
                derivative[index] = 0.0;
                continue;
            }
            const std::size_t left = index ? index - 1 : (boundary == 0 ? points - 1 : 0);
            const std::size_t right = index + 1 < points ? index + 1
                                                         : (boundary == 0 ? 0 : points - 1);
            derivative[index] = scale * (field[left] - 2.0 * field[index] + field[right])
                - decay * field[index];
        }
    };
}

double heatEquationStableStep(const SystemConfiguration& configuration) {
    const auto resolved = resolveHeatEquation(configuration);
    const double diffusivity = parameter(resolved, "diffusivity", 0.01);
    const double decay = parameter(resolved, "decay_rate", 0.0);
    const double spacing = heatEquationMesh(resolved).spacing[0];
    const double fastestRate = 4.0 * diffusivity / (spacing * spacing) + decay;
    return fastestRate == 0.0
        ? (std::numeric_limits<double>::infinity)() : 1.8 / fastestRate;
}

}  // namespace dynsyslab::core
