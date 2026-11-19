#include "ModeWeightPayload.hpp"

#include "ComplexStemPlot.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace dsl {
namespace {

using dynsyslab::core::ArrayType;
using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;
namespace json = dynsyslab::core::json;

constexpr const char* kParentColor = "#2563eb";   // blue: every mode of the parent model
constexpr const char* kKeptColor = "#ef4444";     // red: modes the reduced model kept

bool reconstructionKey(const std::string& name) {
    return name.rfind("reconstruction__", 0) == 0;
}

std::vector<std::complex<double>> complexVector(const NpyArray& array) {
    if (array.type == ArrayType::complex128) return array.complexDoubles();
    std::vector<std::complex<double>> out;
    for (const double value : array.doubles()) out.emplace_back(value, 0.0);
    return out;
}

std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> out; std::stringstream stream(text); std::string line;
    while (std::getline(stream, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
        if (!line.empty()) out.push_back(line);
    }
    return out;
}

// Weights of one model: the eigenvalues it is drawn on and w_oj per scalar observable.
struct ModeWeights {
    std::vector<std::complex<double>> continuous;   // lambda (empty when unavailable)
    std::vector<std::complex<double>> discrete;     // mu (empty when unavailable)
    std::size_t modes = 0;
    std::vector<std::string> names;
    std::vector<std::vector<double>> weights;   // [observable][mode]
};

std::optional<ModeWeights> modeWeights(const NpzRecord& record, const json::Value& metadata) {
    ModeWeights out;
    if (record.contains("continuous_eigenvalues")) out.continuous = complexVector(record.at("continuous_eigenvalues"));
    if (record.contains("discrete_eigenvalues")) out.discrete = complexVector(record.at("discrete_eigenvalues"));
    // Derive the missing plane from the sampling step when it is recorded: mu = exp(lambda dt).
    const double dt = metadata.at("delta_t").number(0.0);
    if (dt > 0 && out.discrete.empty() && !out.continuous.empty())
        for (const auto z : out.continuous) out.discrete.push_back(std::exp(z * dt));
    if (dt > 0 && out.continuous.empty() && !out.discrete.empty())
        for (const auto z : out.discrete) out.continuous.push_back(std::log(z) / dt);
    if (!out.continuous.empty() && !out.discrete.empty() && out.continuous.size() != out.discrete.size())
        out.discrete.clear();
    out.modes = std::max(out.continuous.size(), out.discrete.size());
    const std::size_t modes = out.modes;
    if (modes == 0) return std::nullopt;

    // rms over the training samples of each eigenfunction (scale of phi_j).
    std::vector<double> scale(modes, 1.0);
    if (record.contains("eigenfunction_values")) {
        const auto& phi = record.at("eigenfunction_values");
        if (phi.shape.size() == 2 && phi.shape[1] == modes && phi.shape[0] > 0) {
            const auto values = complexVector(phi);
            std::fill(scale.begin(), scale.end(), 0.0);
            for (std::size_t i = 0; i < phi.shape[0]; ++i)
                for (std::size_t j = 0; j < modes; ++j) scale[j] += std::norm(values[i * modes + j]);
            for (auto& value : scale) value = std::sqrt(value / double(phi.shape[0]));
        }
    }

    // Scalar observables: named one-row reconstructions, else rows of reconstruction_matrix.
    std::vector<std::pair<std::string, std::vector<std::complex<double>>>> rows;
    std::vector<std::string> names;
    if (record.contains("reconstruction_names")) names = splitLines(record.at("reconstruction_names").string());
    for (const auto& [key, array] : record.arrays)
        if (reconstructionKey(key)) {
            const auto name = key.substr(std::string("reconstruction__").size());
            if (std::find(names.begin(), names.end(), name) == names.end()) names.push_back(name);
        }
    for (const auto& name : names) {
        const auto key = "reconstruction__" + name;
        if (!record.contains(key)) continue;
        const auto& C = record.at(key);
        if (C.shape.size() != 2 || C.shape[0] != 1 || C.shape[1] != modes) continue;
        rows.emplace_back(name, complexVector(C));
    }
    if (rows.empty() && record.contains("reconstruction_matrix")) {
        const auto& C = record.at("reconstruction_matrix");
        if (C.shape.size() == 2 && C.shape[1] == modes) {
            const auto values = complexVector(C);
            const auto& observables = metadata.at("observables");
            for (std::size_t r = 0; r < C.shape[0]; ++r) {
                std::string name = "observable " + std::to_string(r);
                if (!observables.isNull() && r < observables.array().size()
                    && !observables.array()[r].string().empty()) name = observables.array()[r].string();
                rows.emplace_back(name, std::vector<std::complex<double>>(
                    values.begin() + std::ptrdiff_t(r * modes), values.begin() + std::ptrdiff_t((r + 1) * modes)));
            }
        }
    }
    if (rows.empty()) return std::nullopt;

    // Observables in the training order when the metadata gives one.
    const auto& order = metadata.at("observables");
    if (!order.isNull()) {
        std::vector<std::string> wanted;
        for (const auto& item : order.array()) wanted.push_back(item.string());
        std::stable_sort(rows.begin(), rows.end(), [&](const auto& a, const auto& b) {
            const auto ia = std::find(wanted.begin(), wanted.end(), a.first) - wanted.begin();
            const auto ib = std::find(wanted.begin(), wanted.end(), b.first) - wanted.begin();
            return ia < ib;
        });
    }
    for (auto& [name, c] : rows) {
        std::vector<double> w(modes);
        for (std::size_t j = 0; j < modes; ++j) w[j] = std::abs(c[j]) * scale[j];
        out.names.push_back(name); out.weights.push_back(std::move(w));
    }
    return out;
}

std::vector<double> totalWeights(const ModeWeights& model) {
    std::vector<double> total(model.modes, 0.0);
    for (const auto& w : model.weights)
        for (std::size_t j = 0; j < total.size(); ++j) total[j] += w[j] * w[j];
    for (auto& value : total) value = std::sqrt(value);
    return total;
}

std::string readableName(std::string name) {
    std::replace(name.begin(), name.end(), '_', ' ');
    return name;
}

}  // namespace

Payload modeWeightPayload(const std::filesystem::path& source, const std::string& view,
                          const NpzRecord& record, const json::Value& metadata) {
    const auto own = modeWeights(record, metadata);
    if (!own) throw std::runtime_error("archive has no eigenvalues with a reconstruction map");

    // A reduced model extracted from a parent: draw the parent spectrum and mark the kept modes.
    std::optional<ModeWeights> parent;
    std::vector<std::size_t> kept;
    const auto parentPath = metadata.at("source_krom").string(metadata.at("source_archive").string());
    const auto& selected = metadata.at("selected_modes");
    if (!parentPath.empty() && !selected.isNull() && !selected.array().empty()) {
        std::error_code missing;
        if (std::filesystem::is_regular_file(std::filesystem::u8path(parentPath), missing)) {
            try {
                const auto parentRecord = dynsyslab::core::readNpz(std::filesystem::u8path(parentPath));
                const auto text = parentRecord.metadataJson();
                parent = modeWeights(parentRecord, text.empty() ? json::Value{} : json::parse(text));
            } catch (const std::exception&) { parent.reset(); }
        }
        if (parent) {
            for (const auto& item : selected.array()) {
                const double index = item.number(-1.0);
                if (index < 0 || index >= double(parent->modes)) { parent.reset(); break; }
                kept.push_back(std::size_t(index));
            }
            // The parent must decompose the same observables.
            if (parent && parent->names != own->names) {
                ModeWeights matched = *parent; matched.names.clear(); matched.weights.clear();
                for (std::size_t o = 0; o < parent->names.size(); ++o) {   // parent's (training) order
                    if (std::find(own->names.begin(), own->names.end(), parent->names[o]) == own->names.end()) continue;
                    matched.names.push_back(parent->names[o]);
                    matched.weights.push_back(parent->weights[o]);
                }
                if (matched.names.empty()) parent.reset(); else *parent = std::move(matched);
            }
        }
    }

    const ModeWeights& drawn = parent ? *parent : *own;
    const auto total = totalWeights(drawn);
    const double heaviest = total.empty() ? 0.0 : *std::max_element(total.begin(), total.end());
    Payload payload{source.u8string(), view, metadata.at("label").string(source.stem().u8string()), 2};

    // One page per plane. The weights belong to the modes, so every page shows the
    // same stems; only the base point of each stem (mu or lambda) changes.
    enum class Plane { discreteFull, discreteZoom, continuous };
    const auto page = [&](Plane plane, const std::string& label) {
        const auto& bases = plane == Plane::continuous ? drawn.continuous : drawn.discrete;
        if (bases.size() != drawn.modes) return;
        const bool discrete = plane != Plane::continuous;
        // View window: the extent of the modes carrying >= 0.1% of the heaviest total
        // weight, plus 10%. The full discrete page always spans the unit circle.
        double reLow = 0, reHigh = 0, imHigh = 0; bool any = false;
        for (std::size_t j = 0; j < total.size(); ++j) {
            if (!(total[j] >= 1e-3 * heaviest) || !std::isfinite(bases[j].real()) || !std::isfinite(bases[j].imag())) continue;
            const auto z = bases[j];
            if (!any) { reLow = reHigh = z.real(); any = true; }
            reLow = std::min(reLow, z.real()); reHigh = std::max(reHigh, z.real());
            imHigh = std::max(imHigh, std::abs(z.imag()));
        }
        if (plane == Plane::discreteFull) {
            reLow = std::min(reLow, -1.0); reHigh = std::max(reHigh, 1.0); imHigh = std::max(imHigh, 1.0);
        }
        if (plane == Plane::continuous) reHigh = std::max(reHigh, 0.0);   // keep the boundary Re = 0 in view
        const double rePad = 0.1 * std::max(reHigh - reLow, 1e-9), imPad = 0.1 * std::max(imHigh, 1e-9);
        const std::array<double, 4> window{reLow - rePad, reHigh + rePad, -(imHigh + imPad), imHigh + imPad};
        const auto inView = [&](std::complex<double> z) {
            return z.real() >= window[0] && z.real() <= window[1] && z.imag() >= window[2] && z.imag() <= window[3];
        };
        std::size_t outside = 0;
        for (const auto& z : bases) if (!inView(z)) ++outside;
        const int pageIndex = static_cast<int>(payload.pageLabels.size());
        payload.pageLabels.push_back(label);
        const auto panel = [&](const std::string& name, const std::vector<double>& weights) {
            ComplexStemPlot plot;
            plot.baseSymbol = discrete ? "mu" : "lambda";
            plot.unitCircle = discrete;                    // |mu| > 1 grows
            plot.imaginaryAxis = !discrete;                // Re lambda > 0 grows
            if (plane != Plane::discreteFull) plot.referenceWindow = window;
            plot.logHeights = false;
            plot.heightLabel = "|c_j| rms(phi_j)";
            for (std::size_t j = 0; j < weights.size(); ++j)
                if (inView(bases[j])) { plot.bases.push_back(bases[j]); plot.heights.push_back(weights[j]); }
            char text[160];
            if (parent) {
                plot.seriesLabel = "parent modes"; plot.color = kParentColor;
                StemSeries keptSeries; keptSeries.label = "kept modes"; keptSeries.color = kKeptColor;
                double all = 0, share = 0;
                for (const double w : weights) all += w * w;
                for (const auto j : kept) {
                    if (inView(bases[j])) { keptSeries.bases.push_back(bases[j]); keptSeries.heights.push_back(weights[j]); }
                    share += weights[j] * weights[j];
                }
                plot.overlays.push_back(std::move(keptSeries));
                std::snprintf(text, sizeof text, " (%zu of %zu modes kept, %.1f%% of energy)",
                              kept.size(), drawn.modes, all > 0 ? 100.0 * share / all : 0.0);
            } else {
                plot.seriesLabel = "modes"; plot.color = kParentColor;
                std::snprintf(text, sizeof text, " (%zu modes)", drawn.modes);
            }
            plot.title = "Mode weights: " + name + text;
            if (outside) plot.title += "; " + std::to_string(outside) + " light modes outside view";
            auto figure = complexPlaneStemFigure(plot);
            figure.dashboardPage = pageIndex;
            return figure;
        };
        payload.figures.push_back(panel("total", total));
        for (std::size_t o = 0; o < drawn.names.size(); ++o)
            payload.figures.push_back(panel(readableName(drawn.names[o]), drawn.weights[o]));
    };
    page(Plane::discreteFull, "Discrete-time mu: unit circle");
    page(Plane::discreteZoom, "Discrete-time mu: near the weighted modes");
    page(Plane::continuous, "Continuous-time lambda");
    if (!parent && !selected.isNull() && !parentPath.empty())
        payload.warnings.push_back("parent model " + parentPath + " could not be read; showing this model's own modes");
    return payload;
}

}  // namespace dsl
