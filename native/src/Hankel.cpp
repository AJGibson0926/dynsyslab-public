#include <dynsyslab/core/Hankel.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
std::size_t multiply(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max()/b)
        throw std::invalid_argument("Hankel dimensions overflow");
    return a*b;
}
}

HankelWindow fixedPairHankelWindow(std::size_t samples,std::size_t maxLag,
    std::size_t requestedDelays,std::size_t pairs,std::size_t sampleStride,std::size_t snapshotStride) {
    if (!samples || !sampleStride || !snapshotStride || pairs<2 || requestedDelays<2)
        throw std::invalid_argument("invalid fixed-pair Hankel geometry");
    const auto last=(samples-1)/sampleStride;
    if (maxLag>=last || pairs-1>(last-maxLag-1)/snapshotStride)
        throw std::invalid_argument("not enough samples for requested H-DMD pairs");
    const auto delays=std::min(requestedDelays,last-maxLag-(pairs-1)*snapshotStride);
    if (delays<2) throw std::invalid_argument("not enough samples for two H-DMD delays");
    return {delays,sampleStride,1,snapshotStride,pairs};
}

HankelMatrix buildHankelMatrix(const std::vector<double>& values, std::size_t samples,
    std::size_t inputs, const std::vector<HankelChannel>& channels,
    const HankelWindow& w, bool paired) {
    if (!samples || !inputs || values.size()!=multiply(samples,inputs) || channels.empty()
        || !w.delays || !w.sampleStride || !w.delayStride || !w.snapshotStride)
        throw std::invalid_argument("invalid Hankel input shape or window");
    std::size_t warmup=0;
    for (const auto& c:channels) {
        if (c.input>=inputs || !(c.scale>0.) || !std::isfinite(c.scale))
            throw std::invalid_argument("invalid Hankel channel index or scale");
        warmup=std::max(warmup,c.lag);
    }
    const auto last=(samples-1)/w.sampleStride;
    const auto span=multiply(w.delays-1,w.delayStride);
    if (warmup>last || span>last-warmup || last-warmup-span<std::size_t(paired))
        throw std::invalid_argument("source is too short for the complete Hankel window");
    const auto available=1+(last-warmup-span-std::size_t(paired))/w.snapshotStride;
    if (w.columns>available) throw std::invalid_argument("requested Hankel columns exceed available complete windows");
    HankelMatrix h;
    h.window=w; h.channels=channels; h.sourceSamples=samples;
    h.rows=multiply(w.delays,channels.size()); h.columns=w.columns?w.columns:available;
    h.x.resize(multiply(h.rows,h.columns));
    if (paired) h.y.resize(h.x.size());
    for (std::size_t j=0;j<h.columns;++j) {
        const auto start=warmup+j*w.snapshotStride;
        h.starts.push_back(start*w.sampleStride);
        h.ends.push_back((start+span)*w.sampleStride);
    }
    for (std::size_t d=0;d<w.delays;++d) {
        for (std::size_t c=0;c<channels.size();++c) {
            const auto& channel=channels[c];
            auto tick=(warmup-channel.lag+d*w.delayStride)*w.sampleStride;
            const auto row=(d*channels.size()+c)*h.columns;
            for (std::size_t j=0;j<h.columns;++j) {
                auto value=[&](std::size_t t) {
                    const double v=values[t*inputs+channel.input]/channel.scale;
                    if (!std::isfinite(v)) throw std::invalid_argument("nonfinite Hankel measurement");
                    return v;
                };
                h.x[row+j]=value(tick);
                if (paired) h.y[row+j]=value(tick+w.sampleStride);
                if (j+1<h.columns) tick+=w.snapshotStride*w.sampleStride;
            }
        }
    }
    return h;
}

std::vector<double> HankelWindowLayout::window_at(const TimeSeries& record, std::size_t label) const {
    if (label < labelOffset) throw std::out_of_range("Hankel window starts before the record");
    const std::size_t start = label - labelOffset, s = window.sampleStride, channelsCount = channels.size();
    std::vector<double> h(rows());
    for (std::size_t d = 0; d < window.delays; ++d)
        for (std::size_t c = 0; c < channelsCount; ++c) {
            const auto& channel = channels[c];
            const auto offset = static_cast<std::ptrdiff_t>(s * d * window.delayStride)
                              - static_cast<std::ptrdiff_t>(s * channel.lag);
            const auto row = static_cast<std::ptrdiff_t>(start) + offset;
            if (row < 0 || std::size_t(row) >= record.sampleCount() || channel.input >= record.componentCount())
                throw std::out_of_range("Hankel window extends past the record");
            h[d * channelsCount + c] = record.value(std::size_t(row), channel.input);
        }
    return h;
}

bool HankelWindowLayout::complete(std::size_t samples, std::size_t label) const noexcept {
    if (label < labelOffset || !window.delays) return false;
    const auto start = static_cast<std::ptrdiff_t>(label - labelOffset), s = static_cast<std::ptrdiff_t>(window.sampleStride);
    for (const auto& channel : channels) {
        const auto first = start - s * static_cast<std::ptrdiff_t>(channel.lag);
        const auto last = first + s * static_cast<std::ptrdiff_t>((window.delays - 1) * window.delayStride);
        if (first < 0 || last >= static_cast<std::ptrdiff_t>(samples)) return false;
    }
    return true;
}

json::Value hankelWindowLayoutJson(const HankelWindowLayout& layout) {
    json::Value::Array channels;
    for (std::size_t c = 0; c < layout.channels.size(); ++c)
        channels.emplace_back(json::Value::Object{{"column", double(layout.channels[c].input)},
            {"lag", double(layout.channels[c].lag)}, {"scale", layout.channels[c].scale},
            {"observable", c < layout.channelObservables.size() ? layout.channelObservables[c] : ""}});
    return json::Value::Object{{"causal", layout.alignEnd}, {"delays", double(layout.window.delays)},
        {"sample_stride", double(layout.window.sampleStride)}, {"delay_stride", double(layout.window.delayStride)},
        {"label_offset", double(layout.labelOffset)}, {"channels", channels}, {"rows", double(layout.rows())}};
}

HankelWindowLayout hankelWindowLayoutFromJson(const json::Value& value) {
    HankelWindowLayout layout;
    layout.window.delays = std::size_t(value.at("delays").number());
    layout.window.sampleStride = std::size_t(value.at("sample_stride").number(1));
    layout.window.delayStride = std::size_t(value.at("delay_stride").number(1));
    layout.alignEnd = value.at("causal").boolean(false);
    layout.labelOffset = std::size_t(value.at("label_offset").number(0));
    for (const auto& channel : value.at("channels").array()) {
        layout.channels.push_back({std::size_t(channel.at("column").number()), std::size_t(channel.at("lag").number()),
                                   channel.at("scale").number(1.0)});
        layout.channelObservables.push_back(channel.at("observable").string());
    }
    if (!layout.window.delays || layout.channels.empty())
        throw std::invalid_argument("window layout needs delays and channels");
    return layout;
}
} // namespace dynsyslab::core
