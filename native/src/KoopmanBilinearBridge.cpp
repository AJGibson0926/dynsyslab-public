// KoopmanModel -> discrete target-centred BilinearKrom. Method and citations:
// core/KoopmanBilinearBridge.hpp (Van Loan 1978; Goswami & Paley 2022; Bruder et al.
// 2021; Strasser, Worthmann et al., SafEDMD, arXiv:2402.03145).
#include <dynsyslab/core/KoopmanBilinearBridge.hpp>

#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <variant>

namespace dynsyslab::core {
namespace {

std::vector<double> multiply(const std::vector<double>& M, const std::vector<double>& x, std::size_t rows) {
    const auto columns = x.size();
    std::vector<double> y(rows, 0.0);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < columns; ++j) y[i] += M[i * columns + j] * x[j];
    return y;
}

double norm2(const std::vector<double>& x) {
    double sum = 0.0;
    for (double v : x) sum += v * v;
    return std::sqrt(sum);
}

// Upper-right block of expm([[G, N], [0, G]] dt) = int_0^dt exp(G(dt-s)) N exp(G s) ds.
std::vector<double> vanLoanIntegral(const std::vector<double>& G, const std::vector<double>& N,
                                    std::size_t n, double dt) {
    const std::size_t N2 = 2 * n;
    std::vector<double> block(N2 * N2, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            block[i * N2 + j] = G[i * n + j];
            block[i * N2 + n + j] = N[i * n + j];
            block[(n + i) * N2 + n + j] = G[i * n + j];
        }
    const auto exponential = numerics::matrixExponential(block, N2, dt);
    std::vector<double> out(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) out[i * n + j] = exponential[i * N2 + n + j];
    return out;
}

}  // namespace

BilinearKrom continuousBilinearSplit(const KoopmanModel& model, const std::string& representation) {
    if (!model.spectrum().hasContinuousRates())
        throw std::invalid_argument("Bilinear bridge needs finite continuous rates");
    if (representation == "exact" || !model.controlMaps().count(representation))
        throw std::invalid_argument("Bilinear bridge needs a stored control map named " + representation);
    if (!std::holds_alternative<BilinearControlMap>(model.storedControlMap(representation)))
        throw std::invalid_argument("Control map " + representation + " is not bilinear");
    const auto n = model.splitDimension(), m = model.inputDimension();
    BilinearKrom out;
    out.stateDimension = n; out.inputDimension = m;
    out.A = model.generatorSplit();
    const auto mapAt = [&](const std::vector<double>& w) {
        return model.realControlMap(model.controlMapFromLifted(model.valuesFromSplit(w), representation));
    };
    // The map is real-affine in the split coordinates, so basis evaluations are exact.
    std::vector<double> basis(n, 0.0);
    out.B0 = mapAt(basis);
    out.B.assign(m, std::vector<double>(n * n, 0.0));
    for (std::size_t column = 0; column < n; ++column) {
        basis[column] = 1.0;
        const auto map = mapAt(basis);
        basis[column] = 0.0;
        for (std::size_t input = 0; input < m; ++input)
            for (std::size_t row = 0; row < n; ++row)
                out.B[input][row * n + column] = map[row * m + input] - out.B0[row * m + input];
    }
    return out;
}

KoopmanBilinearBridge bridgeKoopmanBilinear(const KoopmanModel& model, const KoopmanBilinearBridgeOptions& o) {
    if (!(o.sampleTime > 0.0) || !std::isfinite(o.sampleTime))
        throw std::invalid_argument("Bilinear bridge needs a positive sample time");
    if (!(o.equilibriumTolerance >= 0.0) || !std::isfinite(o.equilibriumTolerance))
        throw std::invalid_argument("Bilinear bridge equilibrium tolerance must be non-negative");
    KoopmanBilinearBridge out;
    out.continuous = continuousBilinearSplit(model, o.representation);
    const auto n = out.continuous.stateDimension, m = out.continuous.inputDimension;
    out.targetLifted = o.targetLifted.empty() ? std::vector<double>(n, 0.0) : o.targetLifted;
    out.targetInput = o.targetInput.empty() ? std::vector<double>(m, 0.0) : o.targetInput;
    if (out.targetLifted.size() != n || out.targetInput.size() != m)
        throw std::invalid_argument("Bilinear bridge target has the wrong dimension");
    for (double v : out.targetLifted) if (!std::isfinite(v)) throw std::invalid_argument("Target lifted value must be finite");
    for (double v : out.targetInput) if (!std::isfinite(v)) throw std::invalid_argument("Target input must be finite");

    out.frozenTransition = numerics::matrixExponential(out.continuous.A, n, o.sampleTime);
    out.discreteOffsets = numerics::zeroOrderHold(out.continuous.A, out.continuous.B0, n, m, o.sampleTime).response;
    for (const auto& N : out.continuous.B)
        out.discreteOperators.push_back(vanLoanIntegral(out.continuous.A, N, n, o.sampleTime));

    out.heldTargetInput = std::any_of(out.targetInput.begin(), out.targetInput.end(), [](double v) { return v != 0.0; });
    auto effective = out.frozenTransition;  // bit-for-bit A when u* = 0
    if (out.heldTargetInput)
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t k = 0; k < n * n; ++k) effective[k] += out.targetInput[i] * out.discreteOperators[i][k];
    auto inputMatrix = out.discreteOffsets;
    for (std::size_t i = 0; i < m; ++i) {
        const auto column = multiply(out.discreteOperators[i], out.targetLifted, n);
        for (std::size_t row = 0; row < n; ++row) inputMatrix[row * m + i] += column[row];
    }
    // c = A_eff w* + B0 u* - w*
    out.constantTerm = multiply(effective, out.targetLifted, n);
    const auto forced = multiply(out.discreteOffsets, out.targetInput, n);
    for (std::size_t row = 0; row < n; ++row) out.constantTerm[row] += forced[row] - out.targetLifted[row];
    out.constantNorm = norm2(out.constantTerm);
    out.constantTolerance = o.equilibriumTolerance * std::max(1.0, norm2(out.targetLifted));
    out.targetIsEquilibrium = out.constantNorm <= out.constantTolerance;

    out.model.stateDimension = n; out.model.inputDimension = m;
    out.model.sampleTime = o.sampleTime;
    out.model.A = std::move(effective);
    out.model.B0 = std::move(inputMatrix);
    out.model.B = out.discreteOperators;
    out.model.provenance = "frozen-generator bridge from K-ROM control map '" + o.representation
        + "': A = exp(G dt); B0 = ZOH(b0); B_i = Van Loan first-order-in-u integral";
    validateBilinearKrom(out.model);
    return out;
}

BilinearFit refitBridgeOperators(const KoopmanBilinearBridge& bridge, const BilinearDataset& centred,
    const BilinearFitOptions& options) {
    if (bridge.heldTargetInput)
        throw std::invalid_argument("Operator refit is not available with a nonzero held target input: "
                                    "A_eff itself depends on B_i");
    auto settings = options;
    settings.sampleTime = bridge.model.sampleTime;
    if (settings.provenance.empty())
        settings.provenance = "frozen-A regression of B0 and B_i on centred held-input pairs";
    return fitBilinearKromFrozenA(centred, bridge.model.A, bridge.model.stateDimension,
                                  bridge.model.inputDimension, settings);
}

BilinearConversionDefect conversionDefect(const KoopmanBilinearBridge& bridge, const std::vector<double>& bounds) {
    const auto n = bridge.continuous.stateDimension, m = bridge.continuous.inputDimension;
    if (bounds.size() != m) throw std::invalid_argument("Conversion defect needs one bound per input");
    if (m > 12) throw std::invalid_argument("Conversion defect enumerates at most 12 inputs");
    const double dt = bridge.model.sampleTime;
    BilinearConversionDefect out;
    for (std::size_t mask = 0; mask < (std::size_t{1} << m); ++mask) {
        std::vector<double> theta(m);
        for (std::size_t i = 0; i < m; ++i)
            theta[i] = bridge.targetInput[i] + (((mask >> i) & 1U) ? bounds[i] : -bounds[i]);
        auto generator = bridge.continuous.A;
        auto approximate = bridge.frozenTransition;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t k = 0; k < n * n; ++k) {
                generator[k] += theta[i] * bridge.continuous.B[i][k];
                approximate[k] += theta[i] * bridge.discreteOperators[i][k];
            }
        const auto held = numerics::zeroOrderHold(generator, bridge.continuous.B0, n, m, dt);
        double transition = 0.0;
        for (std::size_t k = 0; k < n * n; ++k) {
            const double d = held.transition[k] - approximate[k];
            transition += d * d;
        }
        out.transition = std::max(out.transition, std::sqrt(transition));
        auto difference = held.response;
        for (std::size_t k = 0; k < n * m; ++k) difference[k] -= bridge.discreteOffsets[k];
        out.input = std::max(out.input, norm2(multiply(difference, theta, n)));
    }
    return out;
}

}  // namespace dynsyslab::core
