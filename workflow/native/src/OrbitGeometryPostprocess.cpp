// Post-processing method orbit_geometry: cross-orbit injectivity and closure of coordinate records
// (core/OrbitGeometry.hpp). Pools every source of the job; each is one orbit in coordinates (for
// example a history_projection through a frozen basis), with its physical states in
// `source_states` when it carries them.
//
//   {"name": "orbit_geometry", "settings": {
//      "role": "training" | "validation",      the summary's key; select stages score validation
//      "axes": [0, 1],                         the plane the orbits are compared in
//      "dimensions": [2, 4],                   lifted separations at crossings (default [2, all])
//      "peaks": {"array": "measurement_values" | "values", "column": 0,
//                "prominence": 0.15, "min_spacing_time": 0.5},   prominence: share of the range
//      "cycle_samples": 400, "cycle": 0,       representative cycle: last complete period
//      "groups": ["inside", "inside", "outside"]}}                 one label per source
//
// Per record: peak count, peak closure (median |z(p_{k+1}) - z(p_k)| / RMS radius), the area
// enclosed by its representative cycle. Per pair of records: proper crossings of the two cycles in
// the plane (coordinates scaled by the range of all cycles), each lifted to the higher coordinates
// and to the physical states (scaled by their ranges), and the nearest separation of the cycles.
// Summary: crossing totals by group pair, closure medians by group and the state ambiguity
// (physical separation at crossings). Published as a report object; arrays cycles
// (records x samples x width) and crossings (rows: record a, record b, x, y, state separation,
// lifted separation per dimension).
#include "GeometryPostprocess.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/OrbitGeometry.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

struct Orbit {
    std::string source, group;
    std::vector<double> coordinates, states, cycle, cycleStates;
    std::size_t width{}, stateWidth{}, peaks{};
    double closure{}, area{};
};

double median(std::vector<double> v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
}

std::vector<double> ranges(const std::vector<const std::vector<double>*>& paths, std::size_t width) {
    std::vector<double> lo(width, std::numeric_limits<double>::infinity()), hi(width, -std::numeric_limits<double>::infinity());
    for (const auto* p : paths)
        for (std::size_t k = 0; k + width <= p->size(); k += width)
            for (std::size_t c = 0; c < width; ++c) { lo[c] = std::min(lo[c], (*p)[k + c]); hi[c] = std::max(hi[c], (*p)[k + c]); }
    std::vector<double> out(width);
    for (std::size_t c = 0; c < width; ++c) out[c] = hi[c] - lo[c] > 1e-300 ? hi[c] - lo[c] : 1.0;
    return out;
}

double scaledDistance(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& scale, std::size_t width) {
    double s = 0;
    for (std::size_t c = 0; c < width; ++c) s += std::pow((a[c] - b[c]) / scale[c], 2);
    return std::sqrt(s);
}

std::vector<double> pointOn(const std::vector<double>& path, std::size_t width, std::size_t segment, double f) {
    std::vector<double> out(width);
    for (std::size_t c = 0; c < width; ++c) out[c] = path[segment * width + c] + f * (path[(segment + 1) * width + c] - path[segment * width + c]);
    return out;
}

}  // namespace

std::filesystem::path runOrbitGeometry(const std::vector<std::filesystem::path>& sources, const J& settings,
                                       const std::filesystem::path& root, const std::string& label) {
    if (sources.size() < 1) throw std::invalid_argument("orbit_geometry needs sources");
    const auto role = settings.at("role").string("training");
    if (role != "training" && role != "validation") throw std::invalid_argument("orbit_geometry role is training or validation");
    const auto& axesSetting = settings.at("axes");
    const std::size_t ax = axesSetting.isNull() ? 0 : std::size_t(axesSetting.array().at(0).number());
    const std::size_t ay = axesSetting.isNull() ? 1 : std::size_t(axesSetting.array().at(1).number());
    const auto& peakSpec = settings.at("peaks");
    const auto samples = std::size_t(settings.at("cycle_samples").number(400));
    const auto which = std::size_t(settings.at("cycle").number(0));
    std::vector<Orbit> orbits;
    for (std::size_t r = 0; r < sources.size(); ++r) {
        const auto source = loadSource(sources[r]);
        const auto raw = c::readNpz(sources[r]);
        Orbit orbit;
        orbit.source = sources[r].generic_u8string();
        orbit.group = settings.at("groups").isNull() ? "all" : settings.at("groups").array().at(r).string();
        orbit.width = source.dimension();
        orbit.coordinates = source.values();
        if (std::max(ax, ay) >= orbit.width) throw std::invalid_argument("orbit_geometry axes exceed the coordinates of " + orbit.source);
        if (raw.contains("source_states") && raw.at("source_states").shape.size() == 2 && raw.at("source_states").shape[0] == source.samples()) {
            orbit.states = raw.at("source_states").doubles();
            orbit.stateWidth = raw.at("source_states").shape[1];
        }
        // The peak signal: an archive column, or the first plotted coordinate.
        const auto arrayName = peakSpec.at("array").string();
        std::vector<double> signal;
        if (arrayName.empty()) for (std::size_t k = 0; k < source.samples(); ++k) signal.push_back(orbit.coordinates[k * orbit.width + ax]);
        else {
            const auto& a = raw.at(arrayName);
            const auto columns = a.shape.size() == 2 ? a.shape[1] : 1, column = std::size_t(peakSpec.at("column").number(0));
            if (a.shape.empty() || a.shape[0] != source.samples() || column >= columns)
                throw std::invalid_argument("orbit_geometry peak array " + arrayName + " does not match the record");
            const auto all = a.doubles();
            for (std::size_t k = 0; k < source.samples(); ++k) signal.push_back(all[k * columns + column]);
        }
        const auto [lo, hi] = std::minmax_element(signal.begin(), signal.end());
        const double dt = uniformDeltaTime(source);
        const auto spacing = std::size_t(std::max(1.0, std::round(peakSpec.at("min_spacing_time").number(0.5) / dt)));
        const auto peaks = c::findPeaks(signal, spacing, peakSpec.at("prominence").number(0.15) * (*hi - *lo));
        orbit.peaks = peaks.size();
        orbit.closure = c::peakClosure(orbit.coordinates, orbit.width, peaks);
        orbit.cycle = c::representativeCycle(orbit.coordinates, orbit.width, peaks, samples, which);
        if (orbit.stateWidth) orbit.cycleStates = c::representativeCycle(orbit.states, orbit.stateWidth, peaks, samples, which);
        std::vector<double> x, y;
        for (std::size_t s = 0; s < samples; ++s) { x.push_back(orbit.cycle[s * orbit.width + ax]); y.push_back(orbit.cycle[s * orbit.width + ay]); }
        orbit.area = c::polygonArea(x, y);
        orbits.push_back(std::move(orbit));
    }
    const auto width = orbits.front().width;
    for (const auto& o : orbits) if (o.width != width) throw std::invalid_argument("orbit_geometry sources differ in coordinate dimension");
    std::vector<std::size_t> dimensions{2, width};
    if (!settings.at("dimensions").isNull()) {
        dimensions.clear();
        for (const auto& d : settings.at("dimensions").array()) dimensions.push_back(std::min(width, std::size_t(d.number())));
    }
    std::vector<const std::vector<double>*> cycles, cycleStates;
    for (const auto& o : orbits) { cycles.push_back(&o.cycle); if (o.stateWidth) cycleStates.push_back(&o.cycleStates); }
    const auto scale = ranges(cycles, width);
    const auto stateScale = cycleStates.empty() ? std::vector<double>{} : ranges(cycleStates, orbits.front().stateWidth);

    J::Array records, pairs;
    std::map<std::string, double> byGroup;
    std::map<std::string, std::vector<double>> closures;
    std::vector<double> ambiguity, table;
    double nearest = std::numeric_limits<double>::infinity();
    for (const auto& o : orbits) {
        records.push_back(J::Object{{"source", o.source}, {"group", o.group}, {"peaks", double(o.peaks)},
            {"peak_closure", std::isfinite(o.closure) ? J(o.closure) : J()}, {"cycle_area", o.area}});
        if (std::isfinite(o.closure)) closures[o.group].push_back(o.closure);
    }
    const auto scaled = [&](const Orbit& o, std::size_t axis) {
        std::vector<double> out;
        for (std::size_t s = 0; s < samples; ++s) out.push_back(o.cycle[s * width + axis] / scale[axis]);
        return out;
    };
    for (std::size_t a = 0; a < orbits.size(); ++a)
        for (std::size_t b = a + 1; b < orbits.size(); ++b) {
            const auto& A = orbits[a]; const auto& B = orbits[b];
            const auto crossings = c::findCurveIntersections(scaled(A, ax), scaled(A, ay), scaled(B, ax), scaled(B, ay), 2e-4);
            std::vector<double> sa, sb;  // cycles in scaled first-two coordinates for the nearest separation
            for (std::size_t s = 0; s < samples; ++s) {
                sa.push_back(A.cycle[s * width + ax] / scale[ax]); sa.push_back(A.cycle[s * width + ay] / scale[ay]);
                sb.push_back(B.cycle[s * width + ax] / scale[ax]); sb.push_back(B.cycle[s * width + ay] / scale[ay]);
            }
            const double separation = c::nearestSeparation(sa, sb, 2);
            nearest = std::min(nearest, separation);
            double worst = 0;
            for (const auto& x : crossings) {
                const auto la = pointOn(A.cycle, width, x.first, x.u), lb = pointOn(B.cycle, width, x.second, x.v);
                double state = std::numeric_limits<double>::quiet_NaN();
                if (A.stateWidth && A.stateWidth == B.stateWidth)
                    state = scaledDistance(pointOn(A.cycleStates, A.stateWidth, x.first, x.u),
                                           pointOn(B.cycleStates, B.stateWidth, x.second, x.v), stateScale, A.stateWidth);
                if (std::isfinite(state)) { ambiguity.push_back(state); worst = std::max(worst, state); }
                table.insert(table.end(), {double(a), double(b), x.x * scale[ax], x.y * scale[ay], state});
                for (const auto d : dimensions) table.push_back(scaledDistance(la, lb, scale, d));
            }
            auto groups = std::vector<std::string>{A.group, B.group};
            std::sort(groups.begin(), groups.end());
            byGroup[groups[0] + "_" + groups[1]] += double(crossings.size());
            pairs.push_back(J::Object{{"a", double(a)}, {"b", double(b)}, {"crossings", double(crossings.size())},
                {"nearest_separation", separation}, {"state_ambiguity_maximum", worst}});
        }
    J::Object groupCounts, groupClosure;
    double total = 0;
    for (const auto& [k, v] : byGroup) { groupCounts[k] = v; total += v; }
    std::vector<double> allClosures;
    for (const auto& [k, v] : closures) { groupClosure[k] = median(v); allClosures.insert(allClosures.end(), v.begin(), v.end()); }
    J::Array dims;
    for (const auto d : dimensions) dims.push_back(double(d));
    const auto finite = [](double v) { return std::isfinite(v) ? J(v) : J(); };
    J summary = J::Object{{"crossings", total}, {"crossings_by_group", groupCounts},
        {"peak_closure_median", finite(median(allClosures))}, {"peak_closure_median_by_group", groupClosure},
        {"state_ambiguity_median", finite(median(ambiguity))},
        {"state_ambiguity_maximum", ambiguity.empty() ? J(0.0) : J(*std::max_element(ambiguity.begin(), ambiguity.end()))},
        {"nearest_separation", finite(nearest)}, {"records", records}, {"pairs", pairs}, {"lifted_dimensions", dims}};

    const auto first = loadSource(sources.front());
    auto meta = json::parse(metadata(first, "orbit-geometry", label, "\"schema\":\"dynsyslab.orbit-geometry/v1\"")).object();
    meta["settings"] = settings;
    meta["role"] = role;
    meta[role] = summary;
    J::Array inputs;
    for (const auto& s : sources) inputs.push_back(s.generic_u8string());
    meta["sources"] = inputs;
    meta["references"] = J::Array{"Takens 1981, Lecture Notes in Math. 898:366-381",
        "Sauer, Yorke and Casdagli 1991, J. Stat. Phys. 65:579-616", "O'Rourke 1998, Computational Geometry in C, section 7.2"};
    c::NpzRecord archive;
    std::vector<double> allCycles;
    for (const auto& o : orbits) allCycles.insert(allCycles.end(), o.cycle.begin(), o.cycle.end());
    archive.arrays["cycles"] = c::NpyArray::fromDoubles({orbits.size(), samples, width}, allCycles);
    const auto columns = 5 + dimensions.size();
    archive.arrays["crossings"] = c::NpyArray::fromDoubles({table.size() / columns, columns}, table);
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, first, c::object_kind::report, "orbit-geometry", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
