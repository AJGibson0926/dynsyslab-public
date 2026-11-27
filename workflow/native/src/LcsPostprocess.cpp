// Post-processing method ftle: finite-time Lyapunov exponent fields of the source's system on a grid
// of initial conditions, their ridges and the variational LCS (core/LagrangianCoherence.hpp),
// archived as a dynsyslab.grid-field/v1 field object that Viewer draws as heat maps.
//
//   {"name": "ftle", "settings": {
//      "domain": {"lower": [...], "upper": [...], "points": 101, "axes": [0, 1], "anchor": [...]},
//                default: the source trajectory's bounds, the first two coordinates
//      "initial_time": 0, "integration_time": 10,       |T|; default the source's duration
//      "directions": ["forward", "backward"],           forward: repelling, backward: attracting
//      "maximum_step": 0.01, "auxiliary_spacing": 0,    RK4 step; 0: gradient from grid neighbours
//      "ridge_threshold": 0.3,                          share of the largest FTLE
//      "lcs": {"enabled": true, "seed_stride": 3, "step": 0.25, "minimum_length": 0,
//              "maximum_length": 0, "failure_length": 0, "separation": 0}}}   two-axis grids
//
// Outputs per direction d: ftle_d, ridge_d (0/1) and, on two-axis grids, lcs_repelling or
// lcs_attracting (0/1 at nodes the strainlines cross). The strainlines themselves are the arrays
// lcs_points (samples x 2), lcs_offsets (curves + 1), lcs_direction (+1 forward, -1 backward),
// lcs_length and lcs_mean_stretching; cauchy_green_eigenvalues_d holds lambda_1..lambda_k per node.
#include "GeometryPostprocess.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/GridFieldArchive.hpp>
#include <dynsyslab/core/LagrangianCoherence.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

std::vector<double> numbers(const J& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number());
    return out;
}

c::DomainGrid grid(const J& spec, const c::ODESystem& system, const SourceData& source) {
    const auto n = system.descriptor().stateDimension;
    auto domain = c::FunctionDomain::physicalState(n);
    domain.names = system.descriptor().stateNames;
    const auto& given = spec.at("points");
    const auto points = std::holds_alternative<double>(given.storage()) ? std::size_t(given.number()) : std::size_t(n <= 2 ? 101 : 41);
    auto out = spec.at("lower").isNull() ? c::DomainGrid::around(domain, source.values(), points)
                                         : c::DomainGrid::box(domain, numbers(spec.at("lower")), numbers(spec.at("upper")), points);
    std::vector<std::size_t> axes;
    if (!spec.at("axes").isNull()) for (const auto a : numbers(spec.at("axes"))) axes.push_back(std::size_t(a));
    else for (std::size_t j = 0; j < std::min<std::size_t>(n, 2); ++j) axes.push_back(j);
    for (std::size_t j = 0; j < n; ++j) out.anchor[j] = (out.axes[j].lower + out.axes[j].upper) / 2;
    if (!spec.at("anchor").isNull()) out.anchor = numbers(spec.at("anchor"));
    std::vector<c::GridAxis> kept;
    for (const auto a : axes) kept.push_back(out.axes.at(a));
    // Per-axis node counts: "points" may be a list.
    if (!spec.at("points").isNull() && !std::holds_alternative<double>(spec.at("points").storage())) {
        const auto counts = numbers(spec.at("points"));
        for (std::size_t a = 0; a < kept.size() && a < counts.size(); ++a) kept[a].points = std::size_t(counts[a]);
    }
    out.axes = std::move(kept);
    out.validate();
    return out;
}

void addOutput(c::GridField& field, const std::string& name, const std::vector<double>& values) {
    const auto old = field.outputs, points = field.grid.size();
    std::vector<std::complex<double>> merged(points * (old + 1));
    for (std::size_t p = 0; p < points; ++p) {
        for (std::size_t r = 0; r < old; ++r) merged[p * (old + 1) + r] = field.values[p * old + r];
        merged[p * (old + 1) + old] = values[p];
    }
    field.values = std::move(merged);
    field.outputs = old + 1;
    field.outputNames.push_back(name);
}

}  // namespace

std::filesystem::path runFtle(const std::filesystem::path& sourcePath, const J& settings, const std::filesystem::path& root,
                              const std::string& label) {
    const auto source = loadSource(sourcePath);
    const std::shared_ptr<const c::ODESystem> system = c::createSystem(systemConfiguration(source));
    const auto g = grid(settings.at("domain"), *system, source);
    const double duration = source.samples() > 1 ? source.time().back() - source.time().front() : 1;
    const double T = std::abs(settings.at("integration_time").number(duration));
    c::FlowMapGridOptions flow;
    flow.initialTime = settings.at("initial_time").number(0);
    flow.maximumStep = settings.at("maximum_step").number(std::min(1e-2, T / 200));
    flow.auxiliarySpacing = settings.at("auxiliary_spacing").number(0);
    std::vector<std::string> directions{"forward", "backward"};
    if (!settings.at("directions").isNull()) {
        directions.clear();
        for (const auto& d : settings.at("directions").array()) directions.push_back(d.string());
    }
    const auto& lcs = settings.at("lcs");
    const bool curves = g.axes.size() == 2 && lcs.at("enabled").boolean(true);
    c::StrainlineOptions strain;
    strain.threshold = lcs.at("threshold").number(settings.at("ridge_threshold").number(0.3));
    strain.seedStride = std::size_t(lcs.at("seed_stride").number(3));
    strain.step = lcs.at("step").number(0.25);
    strain.minimumLength = lcs.at("minimum_length").number(0);
    strain.maximumLength = lcs.at("maximum_length").number(0);
    strain.failureLength = lcs.at("failure_length").number(0);
    strain.separation = lcs.at("separation").number(0);

    c::GridField field;
    field.grid = g;
    field.status.assign(g.size(), c::GridStatus::valid);
    c::NpzRecord extra;
    std::vector<double> points, offsets{0}, direction, length, stretching;
    J::Array summary;
    for (const auto& name : directions) {
        if (name != "forward" && name != "backward") throw std::invalid_argument("ftle directions are forward and backward");
        flow.integrationTime = name == "forward" ? T : -T;
        const auto cg = c::cauchyGreenField(*system, g, flow);
        const auto sigma = cg.ftle();
        for (std::size_t p = 0; p < g.size(); ++p) if (cg.status[p] != c::GridStatus::valid) field.status[p] = cg.status[p];
        addOutput(field, "ftle_" + name, sigma);
        extra.arrays["cauchy_green_eigenvalues_" + name] = c::NpyArray::fromDoubles({g.size(), cg.axes}, cg.eigenvalues);
        double top = 0;
        for (const auto s : sigma) if (std::isfinite(s)) top = std::max(top, s);
        J::Object entry{{"direction", name}, {"integration_time", flow.integrationTime}, {"maximum_ftle", top}};
        if (g.axes.size() == 2) {
            const auto ridges = c::ridgeNodes(g, sigma, strain.threshold);
            addOutput(field, "ridge_" + name, std::vector<double>(ridges.begin(), ridges.end()));
            entry["ridge_nodes"] = double(std::count(ridges.begin(), ridges.end(), 1));
        }
        if (curves) {
            const auto found = c::strainlines(cg, strain);
            const auto mask = c::rasterizeCurves(g, found);
            addOutput(field, name == "forward" ? "lcs_repelling" : "lcs_attracting", std::vector<double>(mask.begin(), mask.end()));
            for (const auto& curve : found) {
                points.insert(points.end(), curve.points.begin(), curve.points.end());
                offsets.push_back(double(points.size() / 2));
                direction.push_back(name == "forward" ? 1 : -1);
                length.push_back(curve.length);
                stretching.push_back(curve.meanStretching);
            }
            entry["lcs_curves"] = double(found.size());
        }
        summary.push_back(entry);
    }
    if (curves) {
        extra.arrays["lcs_points"] = c::NpyArray::fromDoubles({points.size() / 2, 2}, points);
        extra.arrays["lcs_offsets"] = c::NpyArray::fromDoubles({offsets.size()}, offsets);
        extra.arrays["lcs_direction"] = c::NpyArray::fromDoubles({direction.size()}, direction);
        extra.arrays["lcs_length"] = c::NpyArray::fromDoubles({length.size()}, length);
        extra.arrays["lcs_mean_stretching"] = c::NpyArray::fromDoubles({stretching.size()}, stretching);
    }
    c::GridFieldSeries series;
    series.frames = {field};
    series.frameTime = {flow.initialTime};
    auto meta = json::parse(metadata(source, "ftle", label, ""));
    meta.object()["ftle_settings"] = settings;
    meta.object()["ftle"] = summary;
    meta.object()["initial_time"] = flow.initialTime;
    meta.object()["references"] = J::Array{
        "Haller 2001, Physica D 149:248-277", "Shadden, Lekien and Marsden 2005, Physica D 212:271-304",
        "Haller 2011, Physica D 240:574-598", "Farazmand and Haller 2012, Chaos 22:013128",
        "Haller 2015, Annu. Rev. Fluid Mech. 47:137-162"};
    series.metadata = std::move(meta);
    auto archive = c::gridFieldArchive(series);
    for (auto& [key, array] : extra.arrays) archive.arrays[key] = std::move(array);
    return publishObject(root, source, c::object_kind::field, "ftle", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
