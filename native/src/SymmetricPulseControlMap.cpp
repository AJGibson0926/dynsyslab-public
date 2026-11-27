#include <dynsyslab/core/SymmetricPulseControlMap.hpp>

#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>

#include <cmath>
#include <map>
#include <stdexcept>

namespace dynsyslab::core {

using Complex = std::complex<double>;

SymmetricPulseMapSamples estimateSymmetricPulseControlMap(const KoopmanModel& model,
                                                          const std::vector<SymmetricPulseRun>& runs) {
    if (model.inputDimension() != 1) throw std::invalid_argument("symmetric-pulse estimation supports one input");
    std::map<std::size_t, std::map<int, const SymmetricPulseRun*>> pairs;
    for (const auto& run : runs) {
        if (run.polarity != -1 && run.polarity != 1) throw std::invalid_argument("symmetric-pulse polarity must be -1 or +1");
        if (!pairs[run.pairId].emplace(run.polarity, &run).second)
            throw std::invalid_argument("each symmetric-pulse pair needs unique -1 and +1 runs");
    }
    SymmetricPulseMapSamples samples;
    samples.modes = model.liftedDimension();
    samples.stateDimension = model.stateDimension();
    for (const auto& [id, pair] : pairs) {
        if (pair.size() != 2) throw std::invalid_argument("symmetric-pulse pair " + std::to_string(id) + " is incomplete");
        const auto& negative = *pair.at(-1);
        const auto& positive = *pair.at(1);
        if (negative.anchor != positive.anchor)
            throw std::invalid_argument("symmetric-pulse pair " + std::to_string(id) + " does not share an exact anchor");
        if (negative.sampleTime != positive.sampleTime || !(positive.sampleTime > 0.0))
            throw std::invalid_argument("symmetric-pulse pair " + std::to_string(id) + " has inconsistent time steps");
        const double difference = positive.input - negative.input;
        if (difference == 0.0) throw std::invalid_argument("symmetric-pulse pair " + std::to_string(id) + " has no opposing input");
        if (positive.anchor.size() != samples.stateDimension) throw std::invalid_argument("symmetric-pulse state dimension mismatch");
        const auto plus = model.evaluate(positive.response.data());
        const auto minus = model.evaluate(negative.response.data());
        double norm = 0.0;
        samples.anchors.insert(samples.anchors.end(), positive.anchor.begin(), positive.anchor.end());
        for (std::size_t j = 0; j < samples.modes; ++j) {
            samples.maps.push_back((plus[j] - minus[j]) / (difference * positive.sampleTime));
            norm += std::norm(plus[j] - minus[j]);
        }
        samples.responseNorms.push_back(std::sqrt(norm));
        ++samples.pairs;
    }
    if (!samples.pairs) throw std::invalid_argument("symmetric-pulse estimation needs at least one pair");
    return samples;
}

SampledMapFit fitConstantControlMapSamples(const KoopmanModel& model, const SymmetricPulseMapSamples& samples) {
    if (samples.modes != model.liftedDimension()) throw std::invalid_argument("control-map samples do not align with the K-ROM");
    ConstantControlMap map;
    map.modes = samples.modes; map.inputs = 1;
    map.matrix.assign(samples.modes, 0.0);
    for (std::size_t p = 0; p < samples.pairs; ++p)
        for (std::size_t j = 0; j < samples.modes; ++j) map.matrix[j] += samples.maps[p * samples.modes + j] / double(samples.pairs);
    double residual = 0.0, total = 0.0;
    for (std::size_t p = 0; p < samples.pairs; ++p)
        for (std::size_t j = 0; j < samples.modes; ++j) {
            residual += std::norm(samples.maps[p * samples.modes + j] - map.matrix[j]);
            total += std::norm(samples.maps[p * samples.modes + j]);
        }
    map.provenance = "sampled-pointwise-constant-control-map-regression/v1 (symmetric pulses)";
    map.modelIdentity = koopmanModelIdentity(model);
    SampledMapFit fit{map, total > 0 ? std::sqrt(residual / total) : 0.0,
                      std::sqrt(residual / double(samples.pairs * samples.modes))};
    return fit;
}

SampledMapFit fitBilinearControlMapSamples(const KoopmanModel& model, const SymmetricPulseMapSamples& samples, double ridge) {
    const std::size_t m = samples.modes, columns = 1 + m;
    if (samples.pairs < columns) throw std::invalid_argument("a bilinear map fit needs at least modes + 1 pulse pairs");
    std::vector<Complex> design(samples.pairs * columns);
    for (std::size_t p = 0; p < samples.pairs; ++p) {
        const auto psi = model.evaluate(samples.anchors.data() + p * samples.stateDimension);
        design[p * columns] = 1.0;
        for (std::size_t i = 0; i < m; ++i) design[p * columns + 1 + i] = psi[i];
    }
    numerics::LeastSquaresOptions options; options.ridge = ridge;
    const auto beta = numerics::solveLeastSquares(design, samples.pairs, columns, samples.maps, m, options);
    BilinearControlMap map;
    map.modes = m; map.inputs = 1;
    map.offsets.assign(m, 0.0); map.operators.assign(m * m, 0.0);
    for (std::size_t j = 0; j < m; ++j) {
        map.offsets[j] = beta[j];
        for (std::size_t i = 0; i < m; ++i) map.operators[j * m + i] = beta[(1 + i) * m + j];
    }
    double residual = 0.0, total = 0.0;
    for (std::size_t p = 0; p < samples.pairs; ++p)
        for (std::size_t j = 0; j < m; ++j) {
            Complex fit = map.offsets[j];
            for (std::size_t i = 0; i < m; ++i) fit += map.operators[j * m + i] * design[p * columns + 1 + i];
            residual += std::norm(samples.maps[p * m + j] - fit);
            total += std::norm(samples.maps[p * m + j]);
        }
    map.provenance = "sampled-pointwise-bilinear-control-map-regression/v1 (symmetric pulses)";
    map.modelIdentity = koopmanModelIdentity(model);
    return {map, total > 0 ? std::sqrt(residual / total) : 0.0, std::sqrt(residual / double(samples.pairs * m))};
}

}  // namespace dynsyslab::core
