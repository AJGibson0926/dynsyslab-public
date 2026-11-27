#include "dynsyslab_viewer/KromWorkspace.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "ComplexHeatmap.hpp"
#include "KromTraining.hpp"
#include <dynsyslab/core/DomainGrid.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <numeric>
#include <stdexcept>

namespace dsl {
namespace {
constexpr double pi = 3.14159265358979323846;
const char* colors[] = {"#3bbfb8", "#e8745b", "#4e8bdd", "#b384e0", "#e9ba51", "#64be7a"};
std::vector<std::size_t> sourceIndices(const dynsyslab::core::NpzRecord& a,
    const json::Value& m, std::size_t count) {
    std::vector<std::size_t> result(count); std::iota(result.begin(),result.end(),0);
    if (a.contains("source_mode_indices") && a.at("source_mode_indices").elementCount() == count) {
        const auto values = a.at("source_mode_indices").doubles();
        for (std::size_t i = 0; i < count; ++i) result[i] = static_cast<std::size_t>(values[i]);
        return result;
    }
    const auto& e = m.at("analysis_metadata").at("extraction");
    for (const auto* v : {&e.at("derived_retained_source_indices"), &e.at("selected_source_indices"),
        &m.at("analysis_metadata").at("krom").at("retained_indices")})
        if (v->array().size() == count) {
            for (std::size_t i = 0; i < count; ++i) result[i] = static_cast<std::size_t>(v->array()[i].number());
            break;
        }
    return result;
}
std::vector<std::complex<double>> spectrum(const dynsyslab::core::NpzRecord& a,
    const json::Value& m, const dynsyslab::core::KoopmanModel& model) {
    if (model.spectrum().sampleInterval() > 0) return model.spectrum().discreteMultipliers();
    (void)a;
    const double dt = kromSampleTime(m);
    if (!(dt > 0)) throw std::runtime_error("Discrete KROM spectrum requires archived discrete_eigenvalues or a positive delta_t");
    std::vector<std::complex<double>> result;
    for (const auto value : model.continuousEigenvalues()) result.push_back(std::exp(dt * value));
    return result;
}
void coordinates(Figure& f, const KromTraining& training, bool map, bool fixed = false) {
    const auto& current = training.coordinates;
    f.coordinates = current.coordinates;
    if (current.coordinates == "causal-history") return;
    const auto other = displayCoordinateTransform(training.metadata,training.states.shape[1],
        current.coordinates == "dimensionless" ? "dimensional" : "dimensionless");
    if (other.coordinates == current.coordinates) return;
    Figure::CoordinateAxes axes; axes.name = other.coordinates;
    axes.xLabel = fixed ? f.xLabel : map ? other.stateLabels[0] : training.sampleIndex ? "sample index" : other.timeLabel;
    axes.yLabel = fixed || !map ? f.yLabel : other.stateLabels[1];
    if (!fixed) {
        axes.xScale = map ? other.stateMultipliers[0] : training.sampleIndex ? 1 : other.timeMultiplier;
        axes.yScale = map ? other.stateMultipliers[1] : 1;
    }
    axes.colorLabel = f.colorLabel; f.alternateCoordinates = axes;
}
// The pseudospectrum as a heat map of log10 of its minimum residual under the eigenvalues, when its grid
// is a full rectangular lattice of the mu plane (Trefethen & Embree 2005; Colbrook & Townsend 2024).
// Other grids are not drawn.
std::optional<Trace> pseudospectrumTrace(const dynsyslab::core::Pseudospectrum& p) {
    std::vector<double> xs, ys;
    for (const auto z : p.grid) { xs.push_back(z.real()); ys.push_back(z.imag()); }
    std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    std::sort(ys.begin(), ys.end()); ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    if (xs.size() < 2 || ys.size() < 2 || xs.size() * ys.size() != p.grid.size()) return std::nullopt;
    Trace t; t.label = "pseudospectrum"; t.heatmap = true; t.encoding = "scalar"; t.x = xs; t.y = ys;
    t.zShape = {static_cast<int>(ys.size()), static_cast<int>(xs.size())};
    t.z.assign(p.grid.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t k = 0; k < p.grid.size(); ++k) {
        const auto i = std::size_t(std::lower_bound(xs.begin(), xs.end(), p.grid[k].real()) - xs.begin());
        const auto j = std::size_t(std::lower_bound(ys.begin(), ys.end(), p.grid[k].imag()) - ys.begin());
        t.z[j * xs.size() + i] = std::log10(std::max(p.minimum[k], 1e-300));
    }
    return t;
}
std::string residualLabel(const std::optional<dynsyslab::core::EigenpairResidual>& residual) {
    if (!residual) return {};
    char text[64];
    std::snprintf(text, sizeof text, ", r = %.3g (%s)", residual->value, residual->definition.c_str());
    return text;
}
double component(std::complex<double> value, const std::string& display, unsigned part) {
    return complexComponent(value, display, part);
}
std::string componentLabel(const std::string& psi, const std::string& display, unsigned part) {
    return complexComponentLabel(psi, display, part);
}
}
Payload buildKromWorkspace(const std::filesystem::path& source, const KromWorkspaceOptions& options) {
    const auto& display = options.complexDisplay;
    if (!validComplexDisplay(display))
        throw std::invalid_argument("complex-display must be real-imaginary, magnitude-phase, or domain-coloring");
    const auto archive = composedKrom(source, dynsyslab::core::readNpz(source));
    const auto metadata = json::parse(archive.metadataJson());
    const auto model = loadKromEigenfunctions(source,archive,metadata);
    const auto discrete = spectrum(archive,metadata,model);
    const auto ids = sourceIndices(archive,metadata,discrete.size());
    auto modes = options.modes.value_or(discrete.empty() ? std::vector<std::size_t>{} : std::vector<std::size_t>{0});
    if (options.allModes) { modes.resize(discrete.size()); std::iota(modes.begin(),modes.end(),0); }
    std::sort(modes.begin(),modes.end()); modes.erase(std::unique(modes.begin(),modes.end()),modes.end());
    for (const auto mode : modes) if (mode >= discrete.size()) throw std::out_of_range("KROM mode index is outside the retained spectrum");
    const auto training = readKromTraining(source,archive,metadata,model);
    Payload result{source.u8string(),"krom_dashboard",metadata.at("label").string(source.stem().u8string()),2};
    result.layout = "krom"; result.complexDisplay = display; result.selectedModes = modes;
    result.warnings = training.warnings;
    Figure spec{"krom_spectrum","KROM eigenvalues","Re(\\mu)","Im(\\mu)",true};
    spec.dashboardPage = -2;
    Trace circle; circle.label = "unit circle"; circle.color = "#777777"; circle.lineStyle = "dash";
    for (int i = 0; i <= 360; ++i) { circle.x.push_back(std::cos(2*pi*i/360)); circle.y.push_back(std::sin(2*pi*i/360)); }
    spec.traces.push_back(std::move(circle));
    if (const auto& pseudo = model.spectrum().pseudospectrum(); pseudo && !pseudo->empty()) {
        if (auto heat = pseudospectrumTrace(*pseudo)) {
            spec.colorLabel = "log10 min residual (" + pseudo->definition + ")";
            spec.traces.insert(spec.traces.begin(), std::move(*heat));
        } else result.warnings.push_back("The pseudospectrum grid is not a rectangular lattice; it is not drawn.");
    }
    for (std::size_t i = 0; i < discrete.size(); ++i) {
        // Each eigenvalue with its pair's residual, when the method computed one.
        Trace point; point.label = "\\psi_{" + std::to_string(ids[i]) + "}" + residualLabel(model.spectrum().residual(i));
        point.marker = true; point.selectionId = static_cast<int>(i);
        point.selected = std::binary_search(modes.begin(),modes.end(),i);
        point.color = point.selected ? colors[i%6] : "#9a9a9a";
        point.x = {discrete[i].real()}; point.y = {discrete[i].imag()}; spec.traces.push_back(std::move(point));
    }
    coordinates(spec,training,false,true); result.figures.push_back(std::move(spec));
    if (modes.empty()) { result.pageLabels = {"Select eigenvalues"}; return result; }
    const auto selected = model.extract(modes);
    const auto states = training.states.doubles();
    const auto dimension = training.states.shape[1];
    std::vector<std::size_t> samples, starts{0};
    for (std::size_t t = 0; t+1 < training.starts.size(); ++t) {
        const auto n = training.starts[t+1]-training.starts[t], shown = std::min<std::size_t>(n,100000);
        for (std::size_t i = 0; i < shown; ++i) samples.push_back(training.starts[t]+(shown < 2 ? 0 : i*(n-1)/(shown-1)));
        starts.push_back(samples.size());
    }
    std::vector<std::vector<std::complex<double>>> values(modes.size()), grid(modes.size());
    for (const auto row : samples) {
        const auto v = selected.requiresExternalCoordinates()
            ? selected.evaluateCoordinates(states.data()+row*dimension)
            : selected.evaluate(states.data()+row*dimension);
        for (std::size_t i = 0; i < modes.size(); ++i) values[i].push_back(v[i]);
    }
    std::vector<double> x, y;
    if (dimension == 2) {
        AxisLimits bounds;
        if (selected.requiresExternalCoordinates()) {
            bounds.xMin = bounds.yMin = std::numeric_limits<double>::infinity();
            bounds.xMax = bounds.yMax = -std::numeric_limits<double>::infinity();
            for (std::size_t row = 0; row < training.states.shape[0]; ++row) {
                bounds.xMin = std::min(bounds.xMin,states[row*dimension]);
                bounds.xMax = std::max(bounds.xMax,states[row*dimension]);
                bounds.yMin = std::min(bounds.yMin,states[row*dimension+1]);
                bounds.yMax = std::max(bounds.yMax,states[row*dimension+1]);
            }
            const double xPadding = 0.04 * std::max(bounds.xMax-bounds.xMin,1.0e-12);
            const double yPadding = 0.04 * std::max(bounds.yMax-bounds.yMin,1.0e-12);
            bounds.xMin -= xPadding; bounds.xMax += xPadding;
            bounds.yMin -= yPadding; bounds.yMax += yPadding;
        } else {
            dynsyslab::core::NpzRecord record;
            record.arrays["values"] = training.states;
            record.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
                {training.time.size()},training.time);
            const auto phase = buildNativePayload(
                source,"phase",record,training.metadata).figures.front();
            bounds = automaticAxisLimits(phase);
        }
        // The plotted domain: the physical state, or the K-ROM's own coordinates q when they
        // cannot be computed from the state.
        const auto external = selected.requiresExternalCoordinates();
        auto domainGrid = dynsyslab::core::DomainGrid::box(
            external ? dynsyslab::core::FunctionDomain::supplied(2) : dynsyslab::core::FunctionDomain::physicalState(2),
            {bounds.xMin,bounds.yMin},{bounds.xMax,bounds.yMax},101);
        const auto field = dynsyslab::core::evaluateOnGrid([&](const double* point) {
            return external ? selected.evaluateCoordinates(point) : selected.evaluate(point);
        },modes.size(),domainGrid);
        x = field.grid.axisValues(0); y = field.grid.axisValues(1);
        for (std::size_t k = 0; k < field.grid.size(); ++k)
            for (std::size_t i = 0; i < modes.size(); ++i) grid[i].push_back(field.value(k,i));
        const auto invalid = field.count(dynsyslab::core::GridStatus::undefined);
        if (invalid) result.warnings.push_back(std::to_string(invalid)+" grid points could not be evaluated and are left blank.");
    } else result.warnings.push_back("Phase-space heat maps require a two-state system; training time series are shown.");
    std::vector<double> magnitudeScale(modes.size(),0);
    for (std::size_t i = 0; i < modes.size(); ++i) {
        for (const auto v : values[i]) if (std::isfinite(std::abs(v))) magnitudeScale[i] = std::max(magnitudeScale[i],std::abs(v));
        for (const auto v : grid[i]) if (std::isfinite(std::abs(v))) magnitudeScale[i] = std::max(magnitudeScale[i],std::abs(v));
        if (!(magnitudeScale[i] > 0)) magnitudeScale[i] = 1;
    }
    const auto panel = [&](bool map, std::size_t i, unsigned part) {
            const char* grays[]{"#111827","#26313d","#3f4b57","#5f6b76","#7f8993"};
            const auto psi = "\\psi_{"+std::to_string(ids[modes[i]])+"}";
            const auto label = componentLabel(psi,display,part);
            Figure f{map ? "krom_heatmap" : "krom_time",
                map && display == "domain-coloring" ? psi+" — domain coloring" : label,
                map ? training.coordinates.stateLabels[0] : training.coordinates.timeLabel,
                map ? training.coordinates.stateLabels[1] : label,false};
            f.square = map;
            f.titleColor = colors[modes[i]%6];
            f.colorLabel = map && display == "domain-coloring"
                ? "phase (rad); brightness = magnitude" : map ? label : "";
            if (map) {
                Trace t; t.label = psi; t.heatmap = true; t.x = x; t.y = y;
                t.zShape = {static_cast<int>(y.size()),static_cast<int>(x.size())};
                t.encoding = display == "domain-coloring" ? "domain-color" : "scalar";
                for (const auto v : grid[i]) {
                    t.z.push_back(display == "domain-coloring" ? std::arg(v) : component(v,display,part));
                    if (display == "domain-coloring") t.magnitudes.push_back(std::abs(v)/magnitudeScale[i]);
                }
                f.traces.push_back(std::move(t));
                for (std::size_t trajectory = 0; trajectory + 1 < starts.size(); ++trajectory) {
                    Trace curve; curve.label = training.labels[trajectory];
                    curve.color = grays[trajectory % 5]; curve.lineStyle = "dash";
                    curve.opacity = std::max(0.48, 0.88 - 0.10 * static_cast<double>(trajectory));
                    for (std::size_t j = starts[trajectory]; j < starts[trajectory + 1]; ++j) {
                        const auto row = samples[j];
                        curve.x.push_back(states[row * dimension]);
                        curve.y.push_back(states[row * dimension + 1]);
                    }
                    f.traces.push_back(std::move(curve));
                }
            } else for (std::size_t t = 0; t+1 < starts.size(); ++t) {
                Trace curve; curve.label = training.labels[t]; curve.color = grays[t % 5];
                curve.opacity = std::max(0.48, 0.88 - 0.10 * static_cast<double>(t));
                for (std::size_t j = starts[t]; j < starts[t+1]; ++j) {
                    curve.x.push_back(training.time[samples[j]]); curve.y.push_back(component(values[i][j],display,part));
                }
                f.traces.push_back(std::move(curve));
            }
            coordinates(f,training,map);
            f.dashboardPage = 0;
            return f;
    };
    for (std::size_t i = 0; i < modes.size(); ++i) {
        if (display == "domain-coloring") {
            if (!x.empty()) {
                auto map = panel(true,i,0); map.rowSpan = 2;
                result.figures.push_back(std::move(map));
            }
            result.figures.push_back(panel(false,i,0));
            result.figures.push_back(panel(false,i,1));
        } else for (unsigned part = 0; part < 2; ++part) {
            if (!x.empty()) result.figures.push_back(panel(true,i,part));
            result.figures.push_back(panel(false,i,part));
        }
    }
    return result;
}
}
