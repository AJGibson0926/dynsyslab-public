#include <dynsyslab/core/HankelOperator.hpp>
#include "HistoryCorrelation.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
using detail::CorrelationPlan;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Vector = Eigen::VectorXd;

// Squared Frobenius norm of the one-channel Hankel matrix of signal.
double channelEnergy(const std::vector<double>& signal, std::size_t delays) {
    const std::size_t columns = signal.size() - delays + 1;
    double result = 0.0;
    for (std::size_t sample = 0; sample < signal.size(); ++sample) {
        const std::size_t firstDelay = sample >= columns ? sample - columns + 1 : 0;
        const std::size_t lastDelay = std::min(sample, delays - 1);
        const std::size_t count = lastDelay >= firstDelay ? lastDelay - firstDelay + 1 : 0;
        result += static_cast<double>(count) * signal[sample] * signal[sample];
    }
    return result;
}

Matrix orthonormalize(const Matrix& values) {
    Eigen::HouseholderQR<Matrix> qr(values);
    return qr.householderQ() * Matrix::Identity(values.rows(), values.cols());
}
}  // namespace

// One (segment, channel) block of X is the plain Hankel matrix of that
// channel's trimmed signal t[i] = s[warmup - lag + i] / scale, so X^T w and
// X v are FFT correlations of t with the block's kernel.
struct HankelOperator::Plans {
    struct Block { CorrelationPlan adjoint, forward; };
    std::vector<std::vector<Block>> blocks; // [segment][channel]
};

HankelOperator::~HankelOperator() = default;
HankelOperator::HankelOperator(HankelOperator&&) noexcept = default;
HankelOperator& HankelOperator::operator=(HankelOperator&&) noexcept = default;

HankelOperator::HankelOperator(const std::vector<HankelSegment>& segments,
                               std::vector<HankelChannel> channels, std::size_t delays)
    : delays_(delays), channels_(std::move(channels)), plans_(std::make_unique<Plans>()) {
    if (segments.empty() || channels_.empty() || delays_ < 1)
        throw std::invalid_argument("Hankel operator needs segments, channels and delays");
    std::size_t warmup = 0;
    for (const auto& channel : channels_) {
        if (!(channel.scale > 0.) || !std::isfinite(channel.scale))
            throw std::invalid_argument("invalid Hankel channel scale");
        warmup = std::max(warmup, channel.lag);
    }
    for (const auto& segment : segments) {
        if (!segment.inputs || segment.values.size() != segment.samples * segment.inputs)
            throw std::invalid_argument("invalid Hankel segment shape");
        if (segment.samples < warmup + delays_)
            throw std::invalid_argument("Hankel segment is shorter than its delay window");
        const std::size_t columns = segment.samples - warmup - delays_ + 1;
        const std::size_t span = columns + delays_ - 1;
        std::vector<Plans::Block> blocks;
        blocks.reserve(channels_.size());
        for (const auto& channel : channels_) {
            if (channel.input >= segment.inputs)
                throw std::invalid_argument("Hankel channel index exceeds the segment width");
            std::vector<double> trimmed(span);
            const std::size_t first = warmup - channel.lag;
            for (std::size_t i = 0; i < span; ++i) {
                const double value = segment.values[(first + i) * segment.inputs + channel.input] / channel.scale;
                if (!std::isfinite(value)) throw std::invalid_argument("nonfinite Hankel measurement");
                trimmed[i] = value;
            }
            energy_ += channelEnergy(trimmed, delays_);
            blocks.push_back({CorrelationPlan(trimmed, delays_), CorrelationPlan(trimmed, columns)});
        }
        plans_->blocks.push_back(std::move(blocks));
        segmentColumns_.push_back(columns);
        columns_ += columns;
    }
}

std::vector<double> HankelOperator::applyAdjoint(const std::vector<double>& w, std::size_t k) const {
    const std::size_t C = channels_.size();
    if (w.size() != rows() * k) throw std::invalid_argument("Hankel adjoint operand has the wrong shape");
    std::vector<double> result(columns_ * k, 0.0);
    Vector kernel(static_cast<Eigen::Index>(delays_));
    std::size_t offset = 0;
    for (std::size_t s = 0; s < plans_->blocks.size(); ++s) {
        for (std::size_t column = 0; column < k; ++column)
            for (std::size_t c = 0; c < C; ++c) {
                for (std::size_t d = 0; d < delays_; ++d)
                    kernel[static_cast<Eigen::Index>(d)] = w[(d * C + c) * k + column];
                const Vector part = plans_->blocks[s][c].adjoint.correlate(kernel);
                for (std::size_t j = 0; j < segmentColumns_[s]; ++j)
                    result[(offset + j) * k + column] += part[static_cast<Eigen::Index>(j)];
            }
        offset += segmentColumns_[s];
    }
    return result;
}

std::vector<double> HankelOperator::apply(const std::vector<double>& v, std::size_t k) const {
    const std::size_t C = channels_.size();
    if (v.size() != columns_ * k) throw std::invalid_argument("Hankel operand has the wrong shape");
    std::vector<double> result(rows() * k, 0.0);
    std::size_t offset = 0;
    for (std::size_t s = 0; s < plans_->blocks.size(); ++s) {
        Vector kernel(static_cast<Eigen::Index>(segmentColumns_[s]));
        for (std::size_t column = 0; column < k; ++column) {
            for (std::size_t j = 0; j < segmentColumns_[s]; ++j)
                kernel[static_cast<Eigen::Index>(j)] = v[(offset + j) * k + column];
            for (std::size_t c = 0; c < C; ++c) {
                const Vector part = plans_->blocks[s][c].forward.correlate(kernel);
                for (std::size_t d = 0; d < delays_; ++d)
                    result[(d * C + c) * k + column] += part[static_cast<Eigen::Index>(d)];
            }
        }
        offset += segmentColumns_[s];
    }
    return result;
}

HankelLeftBasis fitHankelLeftBasis(const HankelOperator& hankel, const HankelBasisOptions& options) {
    const std::size_t rows = hankel.rows();
    if (!options.rank || !options.subspaceIterations)
        throw std::invalid_argument("invalid Hankel basis options");
    const std::size_t subspace = std::min(rows, options.rank + options.oversampling);
    if (subspace < options.rank) throw std::invalid_argument("Hankel basis rank exceeds the row dimension");
    if (!(hankel.energy() > 0.)) throw std::invalid_argument("Hankel matrix has zero energy");
    std::mt19937_64 generator(options.seed);
    std::normal_distribution<double> normal;
    Matrix q(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(subspace));
    for (Eigen::Index row = 0; row < q.rows(); ++row)
        for (Eigen::Index column = 0; column < q.cols(); ++column) q(row, column) = normal(generator);
    q = orthonormalize(q);
    // X X^T q, accumulated segment by segment (X_s X_s^T q summed in order).
    auto gram = [&](const Matrix& values) {
        const std::size_t k = static_cast<std::size_t>(values.cols());
        const std::vector<double> flat(values.data(), values.data() + values.size());
        const auto adjoint = hankel.applyAdjoint(flat, k);
        const auto product = hankel.apply(adjoint, k);
        return Matrix(Eigen::Map<const Matrix>(product.data(), values.rows(), values.cols()));
    };
    for (std::size_t iteration = 0; iteration < options.subspaceIterations; ++iteration)
        q = orthonormalize(gram(q));
    const Matrix rayleigh = q.transpose() * gram(q);
    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(rayleigh);
    if (eigensolver.info() != Eigen::Success)
        throw std::runtime_error("Hankel Rayleigh-Ritz solve did not converge");
    Matrix ritz(static_cast<Eigen::Index>(subspace), static_cast<Eigen::Index>(options.rank));
    Vector selected(static_cast<Eigen::Index>(options.rank));
    for (std::size_t mode = 0; mode < options.rank; ++mode) {
        const auto source = static_cast<Eigen::Index>(subspace - 1 - mode);
        ritz.col(static_cast<Eigen::Index>(mode)) = eigensolver.eigenvectors().col(source);
        selected[static_cast<Eigen::Index>(mode)] = eigensolver.eigenvalues()[source];
    }
    Matrix left = q * ritz;
    HankelLeftBasis result;
    result.rows = rows; result.rank = options.rank; result.columns = hankel.columns();
    result.totalEnergy = hankel.energy();
    result.leftSingularVectors.resize(rows * options.rank);
    result.singularValues.resize(options.rank);
    result.cumulativeEnergy.resize(options.rank);
    double captured = 0.0;
    for (std::size_t mode = 0; mode < options.rank; ++mode) {
        Eigen::Index pivot = 0;
        left.col(static_cast<Eigen::Index>(mode)).cwiseAbs().maxCoeff(&pivot);
        if (left(pivot, static_cast<Eigen::Index>(mode)) < 0.0) left.col(static_cast<Eigen::Index>(mode)) *= -1.0;
        const double singular = std::sqrt(std::max(selected[static_cast<Eigen::Index>(mode)], 0.0));
        result.singularValues[mode] = singular;
        captured += singular * singular;
        result.cumulativeEnergy[mode] = captured / result.totalEnergy;
    }
    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t mode = 0; mode < options.rank; ++mode)
            result.leftSingularVectors[row * options.rank + mode]
                = left(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(mode));
    return result;
}

} // namespace dynsyslab::core
