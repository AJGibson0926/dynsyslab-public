#include "CoordinateInverse.hpp"
#include "ObjectPublishing.hpp"
#include "CoordinateChannels.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using J = json::Value;
constexpr std::size_t kNeighbours = 8, kMaximumPairs = 120000;

// Minimal k-d tree over scaled points for k-nearest-neighbour queries.
class KdTree {
public:
    KdTree(std::vector<double> points, std::size_t dimension) : p_(std::move(points)), d_(dimension) {
        const std::size_t n = p_.size() / d_;
        order_.resize(n); std::iota(order_.begin(), order_.end(), std::size_t{0});
        build(0, n);
    }
    // Indices and distances of the k nearest points, nearest first.
    std::vector<std::pair<double, std::size_t>> nearest(const double* q, std::size_t k) const {
        std::priority_queue<std::pair<double, std::size_t>> best;  // max-heap of squared distances
        search(0, order_.size(), q, k, best);
        std::vector<std::pair<double, std::size_t>> out;
        while (!best.empty()) { out.push_back({std::sqrt(best.top().first), best.top().second}); best.pop(); }
        std::reverse(out.begin(), out.end());
        return out;
    }
private:
    struct Node { std::size_t axis; double split; };
    void build(std::size_t lo, std::size_t hi) {
        if (hi - lo <= 16) return;
        std::size_t axis = 0; double spread = -1;
        for (std::size_t a = 0; a < d_; ++a) {
            double mn = INFINITY, mx = -INFINITY;
            for (std::size_t i = lo; i < hi; ++i) { const double v = p_[order_[i] * d_ + a]; mn = std::min(mn, v); mx = std::max(mx, v); }
            if (mx - mn > spread) { spread = mx - mn; axis = a; }
        }
        const std::size_t mid = (lo + hi) / 2;
        std::nth_element(order_.begin() + std::ptrdiff_t(lo), order_.begin() + std::ptrdiff_t(mid), order_.begin() + std::ptrdiff_t(hi),
            [&](std::size_t x, std::size_t y) { return p_[x * d_ + axis] < p_[y * d_ + axis]; });
        nodes_[key(lo, hi)] = {axis, p_[order_[mid] * d_ + axis]};
        build(lo, mid); build(mid, hi);
    }
    static std::size_t key(std::size_t lo, std::size_t hi) { return lo * 1000003u + hi; }
    void search(std::size_t lo, std::size_t hi, const double* q, std::size_t k,
                std::priority_queue<std::pair<double, std::size_t>>& best) const {
        if (hi - lo <= 16) {
            for (std::size_t i = lo; i < hi; ++i) {
                double s = 0; for (std::size_t a = 0; a < d_; ++a) { const double t = p_[order_[i] * d_ + a] - q[a]; s += t * t; }
                if (best.size() < k) best.push({s, order_[i]});
                else if (s < best.top().first) { best.pop(); best.push({s, order_[i]}); }
            }
            return;
        }
        const auto& node = nodes_.at(key(lo, hi));
        const std::size_t mid = (lo + hi) / 2;
        const double delta = q[node.axis] - node.split;
        const bool left = delta < 0;
        search(left ? lo : mid, left ? mid : hi, q, k, best);
        if (best.size() < k || delta * delta < best.top().first) search(left ? mid : lo, left ? hi : mid, q, k, best);
    }
    std::vector<double> p_;
    std::size_t d_;
    std::vector<std::size_t> order_;
    std::unordered_map<std::size_t, Node> nodes_;
};

std::vector<std::string> strings(const J& array) {
    std::vector<std::string> out; for (const auto& item : array.array()) out.push_back(item.string()); return out;
}
}  // namespace

struct CoordinateInverse::Readout {
    std::vector<double> scale, targets;
    std::unique_ptr<KdTree> tree;
};

CoordinateInverse::CoordinateInverse() = default;
CoordinateInverse::~CoordinateInverse() = default;
CoordinateInverse::CoordinateInverse(CoordinateInverse&&) noexcept = default;
CoordinateInverse& CoordinateInverse::operator=(CoordinateInverse&&) noexcept = default;

std::vector<double> CoordinateInverse::apply(const std::vector<double>& values, std::size_t rows) const {
    if (values.size() != rows * inputs_) throw std::invalid_argument("Coordinate inverse input has the wrong width");
    if (kind_ == "exact") return values;
    std::vector<double> out(rows * outputs_, 0.0), q(inputs_);
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t a = 0; a < inputs_; ++a) q[a] = values[r * inputs_ + a] / readout_->scale[a];
        const auto found = readout_->tree->nearest(q.data(), kNeighbours);
        double total = 0;
        for (const auto& [distance, index] : found) {
            if (distance == 0) {  // an exact training pair
                std::copy_n(readout_->targets.begin() + std::ptrdiff_t(index * outputs_), outputs_, out.begin() + std::ptrdiff_t(r * outputs_));
                total = -1; break;
            }
            const double w = 1.0 / distance; total += w;
            for (std::size_t j = 0; j < outputs_; ++j) out[r * outputs_ + j] += w * readout_->targets[index * outputs_ + j];
        }
        if (total > 0) for (std::size_t j = 0; j < outputs_; ++j) out[r * outputs_ + j] /= total;
    }
    return out;
}

c::TimeSeries CoordinateInverse::apply(const c::TimeSeries& series) const {
    const std::vector<double> values(series.valuesData(), series.valuesData() + series.sampleCount() * series.componentCount());
    return {std::vector<double>(series.timeData(), series.timeData() + series.sampleCount()),
            apply(values, series.sampleCount()), outputs_};
}

CoordinateInverse loadCoordinateInverse(const std::filesystem::path& path) {
    const auto record = c::readNpz(path);
    const auto meta = record.metadataJson().empty() ? J{} : json::parse(record.metadataJson());
    const auto& info = meta.at("coordinate_inverse");
    const auto kind = info.at("kind").string();
    if (kind != "exact" && kind != "readout")
        throw std::invalid_argument("This coordinate change retains no inverse: " + info.at("reason").string("not recorded"));
    CoordinateInverse inverse;
    inverse.kind_ = kind; inverse.description_ = info;
    inverse.names_ = strings(info.at("target_names")); inverse.labels_ = strings(info.at("target_labels"));
    inverse.inputs_ = std::size_t(info.at("inputs").number()); inverse.outputs_ = std::size_t(info.at("outputs").number());
    if (kind == "readout") {
        auto readout = std::make_unique<CoordinateInverse::Readout>();
        const auto& inputs = record.at("inverse_inputs");
        readout->targets = record.at("inverse_targets").doubles();
        readout->scale = record.at("inverse_scale").doubles();
        auto points = inputs.doubles();
        for (std::size_t r = 0; r < inputs.shape[0]; ++r)
            for (std::size_t a = 0; a < inverse.inputs_; ++a) points[r * inverse.inputs_ + a] /= readout->scale[a];
        readout->tree = std::make_unique<KdTree>(std::move(points), inverse.inputs_);
        inverse.readout_ = std::move(readout);
    }
    return inverse;
}

void recordCoordinateInverse(const SourceData& source, const std::string& kind, const std::string& reason,
    const c::TimeSeries& z, J& meta, c::NpzRecord& archive) {
    const auto sourceMeta = source.metadataJson.empty() ? J{} : json::parse(source.metadataJson);
    const auto names = timeSeriesColumnNames(source);
    J::Array targetNames, targetLabels;
    const auto& labels = sourceMeta.at("coordinate_labels").array();
    for (std::size_t k = 0; k < names.size(); ++k) {
        targetNames.emplace_back(names[k]);
        targetLabels.emplace_back(labels.size() == names.size() ? labels[k].string() : names[k]);
    }
    J::Object info{{"kind", kind}, {"inputs", double(z.componentCount())}, {"outputs", double(source.dimension())},
                   {"target_names", targetNames}, {"target_labels", targetLabels},
                   {"target_archive", std::filesystem::absolute(source.path).generic_u8string()}};
    if (kind == "none") { info["reason"] = reason; info.erase("outputs"); meta.object()["coordinate_inverse"] = J(std::move(info)); return; }
    if (kind == "exact") { info["method"] = reason; meta.object()["coordinate_inverse"] = J(std::move(info)); return; }
    // readout: pairs (z(t), source(t)) at the coordinates' sample times, thinned evenly.
    const std::size_t n = z.sampleCount(), d = z.componentCount(), w = source.dimension();
    const std::size_t stride = std::max<std::size_t>(1, n / kMaximumPairs);
    std::vector<double> inputs, targets, mean(d, 0.0), scale(d, 0.0);
    std::size_t row = 0, pairs = 0;
    for (std::size_t k = 0; k < n; k += stride) {
        const double t = z.timeData()[k];
        while (row < source.samples() && source.time()[row] < t) ++row;
        if (row == source.samples() || source.time()[row] != t) continue;
        for (std::size_t a = 0; a < d; ++a) inputs.push_back(z.valuesData()[k * d + a]);
        for (std::size_t j = 0; j < w; ++j) targets.push_back(source.values()[row * w + j]);
        ++pairs;
    }
    if (pairs < kNeighbours) { recordCoordinateInverse(source, "none", "too few aligned training pairs for a readout", z, meta, archive); return; }
    for (std::size_t r = 0; r < pairs; ++r) for (std::size_t a = 0; a < d; ++a) mean[a] += inputs[r * d + a] / double(pairs);
    for (std::size_t r = 0; r < pairs; ++r) for (std::size_t a = 0; a < d; ++a) {
        const double t = inputs[r * d + a] - mean[a]; scale[a] += t * t / double(pairs); }
    for (auto& s : scale) s = std::sqrt(s) > 0 ? std::sqrt(s) : 1.0;
    archive.arrays["inverse_inputs"] = c::NpyArray::fromDoubles({pairs, d}, inputs);
    archive.arrays["inverse_targets"] = c::NpyArray::fromDoubles({pairs, w}, targets);
    archive.arrays["inverse_scale"] = c::NpyArray::fromDoubles({d}, scale);
    info["method"] = "inverse-distance weighting over the 8 nearest training pairs, coordinates scaled by their standard deviation";
    info["training_pairs"] = double(pairs);
    info["scope"] = "approximate; valid only where the coordinates stay within the training trajectory's range";
    meta.object()["coordinate_inverse"] = J(std::move(info));
}

void writeCoordinateInverse(const std::filesystem::path& coordinates, const std::filesystem::path& input,
                            const std::filesystem::path& output) {
    const auto inverse = loadCoordinateInverse(coordinates);
    const auto record = c::readNpz(input);
    const auto& values = record.at("values");
    if (values.shape.size() != 2 || values.shape[1] != inverse.inputs())
        throw std::invalid_argument("Input values must be samples x " + std::to_string(inverse.inputs()));
    c::NpzRecord out;
    if (record.contains("time")) out.arrays["time"] = record.at("time");
    out.arrays["values"] = c::NpyArray::fromDoubles({values.shape[0], inverse.outputs()}, inverse.apply(values.doubles(), values.shape[0]));
    J::Array names, labels;
    for (const auto& n : inverse.targetNames()) names.emplace_back(n);
    for (const auto& l : inverse.targetLabels()) labels.emplace_back(l);
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(J::Object{
        {"method", "coordinate-inverse"}, {"coordinate_archive", std::filesystem::absolute(coordinates).generic_u8string()},
        {"inverse", inverse.description()}, {"coordinate_names", names}, {"coordinate_labels", labels}})));
    stampDataObject(out, c::object_kind::timeSeries, "coordinate-inverse", {coordinates});
    c::writeNpzAtomic(output, out);
}

} // namespace dynsyslab::workflow::postprocess_detail

#include <dynsyslab/workflow/CoordinateInverse.hpp>
namespace dynsyslab::workflow {
void invertCoordinates(const std::filesystem::path& coordinateArchive,
                       const std::filesystem::path& input, const std::filesystem::path& output) {
    postprocess_detail::writeCoordinateInverse(coordinateArchive, input, output);
}
}  // namespace dynsyslab::workflow
