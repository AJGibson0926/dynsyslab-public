// Delay-tap stages over core/DelayDiagnostics.hpp (which cites the methods):
//   delay_pairs        two taps per lag from data available at each end: settings "channel" (0),
//                      "lags" (whole samples, [0]), "half_width" (1), "stride" (1). Writes a
//                      time_series whose values are (midpoint, secant) per lag at t_e - h dt, plus
//                      raw_taps (older, newer), curvature and the latency h dt.
//   delay_consistency  a forecast of taps, columns (value, rate) per tap with tap 0 at lag 0 and
//                      "lags" in time units, against the cubic Hermite shift of its current tap;
//                      "reference" (default: the source itself) sets each error scale; "stride".
#include "DelayDiagnosticsPostprocess.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/DelayDiagnostics.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using J = json::Value;

std::size_t whole(const J& value, double fallback, const char* name) {
    const double v = value.number(fallback);
    if (!(v >= 0) || v != std::floor(v) || v > 1e9) throw std::invalid_argument(std::string(name) + " must be a whole number");
    return std::size_t(v);
}
J list(const std::vector<double>& values) {
    J::Array a;
    for (const auto v : values) a.emplace_back(v);
    return J(std::move(a));
}
}  // namespace

std::filesystem::path runDelayPairs(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                    const std::string& label) {
    const auto source = loadSource(path);
    const auto channel = whole(settings.at("channel"), 0, "channel");
    const auto h = whole(settings.at("half_width"), 1, "half_width"), stride = whole(settings.at("stride"), 1, "stride");
    if (channel >= source.dimension() || h < 1 || stride < 1 || source.samples() < 3)
        throw std::invalid_argument("delay_pairs needs a channel of a record, a positive half_width and stride");
    std::vector<std::size_t> lags;
    for (const auto& l : settings.at("lags").array()) lags.push_back(whole(l, 0, "lags"));
    if (lags.empty()) lags.push_back(0);
    std::size_t longest = 0;
    for (const auto l : lags) longest = std::max(longest, l);
    std::vector<double> signal(source.samples());
    for (std::size_t k = 0; k < signal.size(); ++k) signal[k] = source.values()[k * source.dimension() + channel];
    const double dt = source.time()[1] - source.time()[0];
    std::vector<std::size_t> ends;
    for (std::size_t e = longest + 2 * h; e < signal.size(); e += stride) ends.push_back(e);
    if (ends.empty()) throw std::invalid_argument("delay_pairs: the record is shorter than the longest lag");
    const auto pairs = c::scalarDelayPairs(signal, ends, lags, dt, h);
    std::vector<double> time;
    for (const auto e : ends) time.push_back(source.time()[e - h]);
    c::NpzRecord archive;
    const auto n = ends.size(), m = lags.size();
    archive.arrays["time"] = c::NpyArray::fromDoubles({n}, time);
    archive.arrays["values"] = c::NpyArray::fromDoubles({n, 2 * m}, pairs.state);
    archive.arrays["raw_taps"] = c::NpyArray::fromDoubles({n, 2 * m}, pairs.raw);
    archive.arrays["curvature"] = c::NpyArray::fromDoubles({n, m}, pairs.curvature);
    std::vector<double> lagValues(lags.begin(), lags.end());
    auto meta = json::parse(metadata(source, "delay-pairs", label, "")).object();
    meta["delay_pairs"] = J::Object{{"channel", double(channel)}, {"lags_samples", list(lagValues)}, {"half_width", double(h)},
        {"latency", pairs.latency}, {"columns", "(midpoint, secant) per lag"}, {"estimate_time", "t_end - half_width dt"}};
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::timeSeries, "delay-pairs", std::move(archive));
}

std::filesystem::path runDelayConsistency(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                          const std::string& label) {
    const auto source = loadSource(path);
    const auto referencePath = settings.at("reference").string();
    const auto reference = referencePath.empty() ? source : loadSource(std::filesystem::u8path(referencePath));
    std::vector<double> lags;
    for (const auto& l : settings.at("lags").array()) lags.push_back(l.number());
    if (source.dimension() != 2 * lags.size() || reference.dimension() != source.dimension() || reference.samples() != source.samples())
        throw std::invalid_argument("delay_consistency needs (value, rate) per lag in aligned forecast and reference records");
    const auto result = c::delayConsistency(source.time(), source.values(), reference.values(), lags,
                                            whole(settings.at("stride"), 1, "stride"));
    c::NpzRecord archive;
    const auto taps = result.lags.size();
    archive.arrays["delayed_lags"] = c::NpyArray::fromDoubles({taps}, result.lags);
    archive.arrays["value_relative_error"] = c::NpyArray::fromDoubles({taps}, result.value);
    archive.arrays["rate_relative_error"] = c::NpyArray::fromDoubles({taps}, result.rate);
    auto meta = json::parse(metadata(source, "delay-consistency", label, "")).object();
    J::Array errors;
    for (std::size_t j = 0; j < taps; ++j) {
        const auto defined = [](double v) { return std::isfinite(v) ? J(v) : J(); };
        errors.emplace_back(J::Object{{"lag", result.lags[j]}, {"value", defined(result.value[j])}, {"rate", defined(result.rate[j])}});
    }
    meta["delay_consistency"] = J::Object{{"lags", list(lags)}, {"errors", J(std::move(errors))},
        {"reference", referencePath.empty() ? J("source") : J(referencePath)}};
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::report, "delay-consistency", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
