// JSON result documents as Viewer payloads (DocumentPayload.hpp).
#include "DocumentPayload.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dsl {
namespace {
namespace json = dynsyslab::core::json;

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.isNull() ? std::nan("") : item.number());
    return out;
}

Trace curve(std::string label, std::vector<double> x, std::vector<double> y, std::string color,
            std::string style = "solid", double width = 1.6) {
    Trace t; t.label = std::move(label); t.x = std::move(x); t.y = std::move(y); t.color = std::move(color);
    t.lineStyle = std::move(style); t.lineWidthPoints = width;
    return t;
}

// Port of the periodic comparison pages: values are never clipped, renormalized or smoothed.
Payload periodicComparison(const json::Value& comparison, const std::string& view, const std::string& source) {
    const bool byFrequency = view != "parameter";
    const auto xs = numbers(comparison.at(byFrequency ? "validation_frequencies" : "validation_parameters"));
    if (std::any_of(xs.begin(), xs.end(), [](double x) { return std::isnan(x); }))
        throw std::invalid_argument("the comparison does not supply the requested abscissa");
    std::vector<std::string> families;
    for (const auto& family : comparison.at("validation_families").array()) families.push_back(family.string());
    const std::size_t n = xs.size();
    const auto angles = numbers(comparison.at("angles"));
    const auto centers = numbers(comparison.at("packet_centers"));
    constexpr const char* colors[] = {"#B57B00", "#0072B2", "#009E73", "#E69F00"};
    constexpr const char* styles[] = {"dash-dot", "dash", "dot", "solid"};
    Payload payload{source, "periodic_comparison", comparison.at("title").string("Periodic observable spectra"), 2};
    std::vector<std::string> order;
    for (const auto& family : families) if (std::find(order.begin(), order.end(), family) == order.end()) order.push_back(family);
    for (const auto& reference : comparison.at("references").array()) {
        std::vector<const json::Value*> matches;
        for (const auto& row : comparison.at("comparisons").array())
            if (row.at("channel").string() == reference.at("channel").string()) matches.push_back(&row);
        for (const auto& family : order) {
            std::vector<std::size_t> ids;
            for (std::size_t j = 0; j < n; ++j) if (families[j] == family) ids.push_back(j);
            std::sort(ids.begin(), ids.end(), [&](std::size_t a, std::size_t b) { return xs[a] < xs[b]; });
            std::vector<double> x; for (const auto j : ids) x.push_back(xs[j]);
            for (std::size_t k = 0; k < centers.size(); ++k) {
                const int page = static_cast<int>(payload.pageLabels.size());
                const auto name = reference.at("name").string();
                payload.pageLabels.push_back(name + " | " + family + " | angle=" + std::to_string(centers[k]));
                const auto packets = [&](const json::Value& row) {
                    const auto rms = numbers(row.at("packet_rms")); std::vector<double> y;
                    for (const auto j : ids) y.push_back(rms.at(k * n + j));
                    return y;
                };
                Figure density{"fourier", name + ": spectral measure (training bank)", "Spectral angle (rad/sample)",
                               "Raw spectral density", false};
                density.valueScale = "symlog"; density.dashboardPage = page;
                density.traces.push_back(curve("Fourier reference", angles, numbers(reference.at("density")), "#C52B24", "solid", 2.1));
                Figure amplitude{"phase", name + ": packet at angle " + std::to_string(centers[k]) + " | " + family,
                                 byFrequency ? "Orbit angular frequency (rad/model time)" : "Orbit parameter",
                                 "RMS packet amplitude", false};
                amplitude.valueScale = "linear"; amplitude.dashboardPage = page;
                amplitude.traces.push_back(curve("Fourier reference", x, packets(reference), "#C52B24", "solid", 2.1));
                for (std::size_t i = 0; i < matches.size(); ++i) {
                    const auto label = matches[i]->at("model").string() + " (rank "
                        + json::dump(matches[i]->at("rank")) + ")";
                    density.traces.push_back(curve(label, angles, numbers(matches[i]->at("density")), colors[i % 4], styles[i % 4]));
                    amplitude.traces.push_back(curve(label, x, packets(*matches[i]), colors[i % 4], styles[i % 4]));
                }
                payload.figures.push_back(std::move(density));
                payload.figures.push_back(std::move(amplitude));
            }
        }
    }
    return payload;
}

// Phase-cover tower: level j covers the base circle p^j times; colours show mapped base phase.
// display.elapsed_base_periods (a list) animates the exact clock: frame t has base phase
// base + t (mod 1) and the tower residue advanced by the floor(base + t) completed base turns
// (carries into the higher digits), as the solenoid flow does.
void phaseCoverPage(const json::Value& geometry, double base, int digit0, int digit1, int page, Payload& payload) {
    const auto prime = static_cast<int>(geometry.at("prime").number());
    const auto samples = static_cast<std::size_t>(geometry.at("samples").number());
    const auto& orbits = geometry.at("orbits").array();
    for (std::size_t j = 0; j < orbits.size(); ++j) {
        const auto values = numbers(orbits[j].at("values"));
        const auto dimension = static_cast<std::size_t>(orbits[j].at("state_dimension").number(2));
        if (dimension != 2 || values.size() != 2 * samples) throw std::invalid_argument("phase covers need planar sampled orbits");
        const int count = static_cast<int>(std::lround(std::pow(prime, static_cast<double>(j))));
        Trace path; path.label = "orbit"; path.colorValues.reserve(samples);
        for (std::size_t i = 0; i <= samples; ++i) {
            const auto k = i % samples;
            path.x.push_back(values[2 * k]); path.y.push_back(values[2 * k + 1]);
            path.colorValues.push_back(std::fmod(count * (static_cast<double>(k) + 0.5) / samples, 1.0));
        }
        path.encoding = "time";  // colour by value (here the mapped base phase)
        Figure figure{"phase", "omega_0/" + std::to_string(count) + ": " + std::to_string(count) + " compatible phase point(s)",
                      "x_1", "x_2", true};
        figure.colorLabel = "mapped base phase (turns)"; figure.dashboardPage = page;
        figure.traces.push_back(std::move(path));
        const int selected = j == 0 ? 0 : digit0 + (j == 2 ? prime * digit1 : 0);
        for (int residue = 0; residue < count; ++residue) {
            const double turn = std::fmod((base + residue) / count, 1.0), index = turn * samples;
            const auto low = static_cast<std::size_t>(index) % samples;
            const double fraction = index - std::floor(index);
            Trace point; point.marker = true; point.label = "phase point " + std::to_string(residue);
            point.markerShape = residue == selected ? "" : "open-circle"; point.color = "#17232E";
            point.x = {(1 - fraction) * values[2 * low] + fraction * values[2 * ((low + 1) % samples)]};
            point.y = {(1 - fraction) * values[2 * low + 1] + fraction * values[2 * ((low + 1) % samples) + 1]};
            figure.traces.push_back(std::move(point));
        }
        payload.figures.push_back(std::move(figure));
    }
}

Payload phaseCover(const json::Value& geometry, const std::string& source) {
    const auto prime = static_cast<int>(geometry.at("prime").number());
    const auto& display = geometry.at("display");
    const double base = display.at("base_phase_turns").number(1.0 / 6.0);
    std::vector<int> digits;
    for (const auto& d : display.at("digits").array()) digits.push_back(static_cast<int>(d.number()));
    while (digits.size() < 2) digits.push_back(0);
    for (const auto d : digits) if (d < 0 || d >= prime) throw std::invalid_argument("phase-cover digits must lie in [0, prime)");
    Payload payload{source, "phase_cover", "Phase-cover tower, p=" + std::to_string(prime), 3};
    const auto& elapsed = display.at("elapsed_base_periods").array();
    if (elapsed.empty()) { phaseCoverPage(geometry, base, digits[0], digits[1], -1, payload); return payload; }
    const int square = prime * prime;
    for (std::size_t frame = 0; frame < elapsed.size(); ++frame) {
        const double total = base + elapsed[frame].number();
        const auto turns = static_cast<long long>(std::floor(total));
        const auto residue = static_cast<int>(((digits[0] + prime * digits[1] + turns) % square + square) % square);
        payload.pageLabels.push_back("t/T0 = " + std::to_string(elapsed[frame].number()) + "; residues "
                                     + std::to_string(residue % prime) + ", " + std::to_string(residue));
        phaseCoverPage(geometry, total - static_cast<double>(turns), residue % prime, residue / prime,
                       static_cast<int>(frame), payload);
    }
    payload.animation = Payload::Animation{};
    return payload;
}

}  // namespace

bool isDocumentPayload(const json::Value& document) {
    const auto schema = document.at("schema").string();
    return document.at("method").string() == "periodic-spectral-comparison"
        || schema == "dynsyslab.phase-cover-geometry/v1" || schema == "dynsyslab.research-figure/v1";
}

bool isDocumentPayloadFile(const std::filesystem::path& path) {
    if (path.extension() != L".json") return false;
    try { return isDocumentPayload(json::load(path.wstring())); } catch (...) { return false; }
}

Payload buildDocumentPayload(const json::Value& document, const std::string& view, const std::string& source) {
    const auto& data = document.at("result").isNull() ? document : document.at("result");
    if (data.at("method").string() == "periodic-spectral-comparison") return periodicComparison(data, view, source);
    if (data.at("schema").string() == "dynsyslab.phase-cover-geometry/v1") return phaseCover(data, source);
    if (data.at("schema").string() == "dynsyslab.research-figure/v1") {
        auto payload = buildResearchFigurePayload(data);
        payload.source = source;
        return payload;
    }
    throw std::invalid_argument("this JSON document has no Viewer presentation");
}

Payload buildDocumentPayload(const std::filesystem::path& source, const std::string& view) {
    return buildDocumentPayload(json::load(source.wstring()), view, std::filesystem::absolute(source).u8string());
}

}  // namespace dsl
