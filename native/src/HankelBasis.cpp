#include <dynsyslab/core/HankelRepresentations.hpp>
#include <dynsyslab/core/HankelBasis.hpp>

#include <stdexcept>

namespace dynsyslab::core {

HankelCoordinateScaling hankelCoordinateScaling(const std::string& name) {
    if (name.empty() || name == "energy") return HankelCoordinateScaling::energy;
    if (name == "whitened") return HankelCoordinateScaling::whitened;
    throw std::invalid_argument("coordinate_scaling must be energy or whitened");
}

std::string hankelCoordinateScalingName(HankelCoordinateScaling scaling) {
    return scaling == HankelCoordinateScaling::energy ? "energy" : "whitened";
}

HankelBasis::HankelBasis(HankelWindowLayout layout, std::vector<double> modes, std::vector<double> singularValues,
                         std::vector<double> mean, DataOrigin origin)
    : layout_(std::move(layout)), modes_(std::move(modes)), singular_(std::move(singularValues)),
      mean_(std::move(mean)), origin_(std::move(origin)) {
    const std::size_t rows = layout_.rows();
    if (singular_.empty() || modes_.size() != rows * singular_.size() || (!mean_.empty() && mean_.size() != rows))
        throw std::invalid_argument("Hankel basis needs rows x rank modes, rank singular values and a rows mean");
    for (const double value : singular_)
        if (!(value > 0.0)) throw std::invalid_argument("Hankel basis singular values must be positive");
}

std::vector<double> HankelBasis::coordinates(const double* window, HankelCoordinateScaling scaling) const {
    const std::size_t rows = layout_.rows(), rank = singular_.size(), channels = layout_.channels.size();
    std::vector<double> q(rank, 0.0);
    for (std::size_t r = 0; r < rows; ++r) {
        const double value = window[r] / layout_.channels[r % channels].scale - (mean_.empty() ? 0.0 : mean_[r]);
        for (std::size_t j = 0; j < rank; ++j) q[j] += modes_[r * rank + j] * value;
    }
    if (scaling == HankelCoordinateScaling::whitened)
        for (std::size_t j = 0; j < rank; ++j) q[j] /= singular_[j];
    return q;
}

std::vector<double> HankelBasis::coordinates(const TimeSeries& record, const std::vector<std::size_t>& labels,
                                             HankelCoordinateScaling scaling) const {
    std::vector<double> out;
    out.reserve(labels.size() * rank());
    for (const auto label : labels) {
        const auto window = layout_.window_at(record, label);
        const auto q = coordinates(window.data(), scaling);
        out.insert(out.end(), q.begin(), q.end());
    }
    return out;
}

json::Value HankelBasis::describe() const {
    json::Value::Array channels;
    for (std::size_t c = 0; c < layout_.channels.size(); ++c)
        channels.emplace_back(json::Value::Object{{"column", double(layout_.channels[c].input)},
            {"lag", double(layout_.channels[c].lag)}, {"scale", layout_.channels[c].scale},
            {"observable", c < layout_.channelObservables.size() ? layout_.channelObservables[c] : ""}});
    return json::Value::Object{{"schema", kHankelBasisSchema}, {"origin", dataOriginJson(origin_)},
        {"rank", double(rank())}, {"centered", !mean_.empty()},
        {"layout", json::Value::Object{{"delays", double(layout_.window.delays)},
            {"sample_stride", double(layout_.window.sampleStride)}, {"delay_stride", double(layout_.window.delayStride)},
            {"causal", layout_.alignEnd}, {"label_offset", double(layout_.labelOffset)}, {"channels", channels}}},
        {"form", "q = U^T (h - mean); xi = Sigma^-1 q"}};
}

HankelBasis hankelBasisFromArchive(const json::Value& block, std::vector<double> modes,
                                   std::vector<double> singularValues, std::vector<double> mean) {
    if (block.at("schema").string() != kHankelBasisSchema) throw std::invalid_argument("not a Hankel basis block");
    return HankelBasis(hankelWindowLayoutFromJson(block.at("layout")), std::move(modes), std::move(singularValues),
                       block.at("centered").boolean(false) ? std::move(mean) : std::vector<double>{},
                       dataOriginFromJson(block.at("origin")));
}

}  // namespace dynsyslab::core
