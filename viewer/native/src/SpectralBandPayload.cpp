// The krom "spectral_measures" view (SpectralBandPayload.hpp).
#include "SpectralBandPayload.hpp"

#include "ComplexStemPlot.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <map>
#include <utility>

namespace dsl {
namespace {
namespace c = dynsyslab::core;
const char* bandColors[] = {"#e8745b", "#4e8bdd", "#b384e0", "#e9ba51", "#64be7a", "#3bbfb8"};
// Bands of one observable often coincide (an exact reference and its fit): each band its own dash.
const char* bandStyles[] = {"solid", "dash", "dot", "dash-dot"};

// The band parameter of a point: its angle (unit circle) or frequency (imaginary axis).
double parameterOf(const c::SpectralBand& band, std::complex<double> mu, std::complex<double> lambda) {
    return band.support == c::band_support::imaginaryAxis ? lambda.imag() : std::arg(mu);
}

std::string outputName(const std::string& observable, std::size_t output, std::size_t outputs) {
    return outputs > 1 ? observable + "[" + std::to_string(output) + "]" : observable;
}
}  // namespace

Payload spectralMeasuresPayload(const std::filesystem::path& source, const std::string& view) {
    const auto model = c::loadKoopmanModel(source);
    const auto& spectrum = model.spectrum();
    Payload payload{source.u8string(), view, source.stem().u8string(), 1};
    if (spectrum.measures().empty()) {
        payload.warnings.push_back("The spectrum holds no spectral measures"
                                   + std::string(spectrum.bands().empty() ? " and no bands." : "; its bands carry none."));
        return payload;
    }
    std::vector<std::complex<double>> mu, lambda;
    try { mu = spectrum.discreteMultipliers(); } catch (const std::exception&) {
        payload.warnings.push_back("The spectrum has no sample interval: atoms are not placed on the unit circle.");
    }
    if (spectrum.hasContinuousRates()) lambda = spectrum.continuousRates();
    // Measures by observable and output, in the spectrum's order.
    std::vector<std::pair<std::string, std::size_t>> keys;
    std::map<std::pair<std::string, std::size_t>, std::vector<const c::SpectralMeasure*>> byKey;
    std::map<std::string, std::size_t> outputs;
    for (const auto& m : spectrum.measures()) {
        const std::pair<std::string, std::size_t> key{m.observable, m.output};
        if (!byKey.count(key)) keys.push_back(key);
        byKey[key].push_back(&m);
        outputs[m.observable] = std::max(outputs[m.observable], m.output + 1);
    }
    std::map<std::string, const char*> colorOf, styleOf;
    for (const auto& b : spectrum.bands()) {
        styleOf[b.name] = bandStyles[colorOf.size() % 4];
        colorOf[b.name] = bandColors[colorOf.size() % 6];
    }
    const double gap = std::numeric_limits<double>::quiet_NaN();
    for (const auto& key : keys) {
        const auto name = outputName(key.first, key.second, outputs[key.first]);
        std::vector<c::SpectralAtom> atoms;
        try { atoms = spectrum.atoms(key.first, key.second); } catch (const std::exception&) {}
        // 3-D: atoms on the unit circle and the unit-circle bands' densities around it.
        ComplexStemPlot plot;
        plot.title = "Spectral measure of " + name + " on the unit circle";
        plot.baseSymbol = "mu"; plot.unitCircle = true;
        plot.heightLabel = "atom mass |c_j|^2, density"; plot.seriesLabel = "atoms |c_j|^2";
        for (const auto& atom : atoms)
            if (const auto at = spectrum.modeIndex(atom.mode); at && *at < mu.size()) {
                plot.bases.push_back(mu[*at]); plot.heights.push_back(atom.mass);
            }
        // 2-D: each band's density against its parameter, with the atoms as stems.
        Figure flat{"time", "Spectral measure of " + name + " by band", "", "density, atom mass", false};
        bool unitCircleBand = false, frequencyBand = false;
        for (const auto* m : byKey[key]) {
            const auto* band = spectrum.band(m->band);
            if (!band || band->grid.size() != m->density.size() || band->grid.empty()) continue;
            const bool circle = band->support == c::band_support::unitCircle;
            (circle ? unitCircleBand : frequencyBand) = true;
            char mass[48]; std::snprintf(mass, sizeof mass, " (mass %.4g)", m->totalMass);
            const std::string label = m->band + (band->method.empty() ? "" : " [" + band->method + "]") + mass;
            if (circle) {
                StemSeries curve; curve.label = label; curve.color = colorOf[m->band];
                for (std::size_t k = 0; k < band->grid.size(); ++k) {
                    curve.bases.push_back(std::polar(1.0, band->grid[k])); curve.heights.push_back(m->density[k]);
                }
                plot.curves.push_back(std::move(curve));
            }
            Trace line; line.label = label; line.color = colorOf[m->band]; line.role = "data";
            line.x = band->grid; line.y = m->density;
            if (band->grid.size() == 1) { line.marker = true; line.lineStyle = "none"; line.markerShape = "open-circle"; }
            else line.lineStyle = styleOf[m->band];
            flat.traces.push_back(std::move(line));
        }
        if (unitCircleBand || !plot.bases.empty()) payload.figures.push_back(complexPlaneStemFigure(plot));
        if (!flat.traces.empty()) {
            flat.xLabel = unitCircleBand && frequencyBand ? "band parameter (angle or frequency)"
                          : unitCircleBand ? "angle theta (rad)" : "frequency (rad / time)";
            const auto* first = spectrum.band(byKey[key].front()->band);
            Trace stems; stems.label = "atoms |c_j|^2"; stems.color = "#0f766e"; stems.role = "data";
            for (const auto& atom : atoms) {
                const auto at = spectrum.modeIndex(atom.mode);
                if (!at || !first) continue;
                const auto m = *at < mu.size() ? mu[*at] : std::complex<double>(gap, gap);
                const auto l = *at < lambda.size() ? lambda[*at] : std::complex<double>(gap, gap);
                const double p = parameterOf(*first, m, l);
                if (!std::isfinite(p)) continue;
                stems.x.insert(stems.x.end(), {p, p, gap}); stems.y.insert(stems.y.end(), {0.0, atom.mass, gap});
            }
            if (!stems.x.empty()) flat.traces.push_back(std::move(stems));
            payload.figures.push_back(std::move(flat));
        }
    }
    if (payload.figures.empty()) payload.warnings.push_back("No measure lies on a band with a grid.");
    return payload;
}

}  // namespace dsl
