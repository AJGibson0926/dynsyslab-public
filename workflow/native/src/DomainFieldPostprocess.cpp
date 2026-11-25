// Postprocess method domain_field: a function sampled on a grid of its domain, or a coordinate map
// acting on it, archived as dynsyslab.grid-field/v1 for Viewer (core/DomainGrid.hpp,
// core/CoordinateMaps.hpp, core/GridFieldArchive.hpp).
//
//   {"name": "domain_field", "settings": {
//      "function": {"observable": "angle"} | {"observables": ["echo", "echo_dot"]}
//                | {"expression": "x_1*cos(x_2)", "variables": [...]}       (default: state names)
//                | {"library": "libration_phase", "reference_state": [...]}
//                | {"krom": "path.npz", "modes": [0, 1]},
//      "domain": {"lower": [...], "upper": [...], "points": 101, "axes": [0, 1], "anchor": [...]},
//                default: the source trajectory's bounds, every coordinate up to three
//      "operation": "evaluate" (default) | "koopman_action" | "push_forward" | "pull_back"
//                 | "round_trip" | "fold_set",
//      "time_step": 0.1, "steps": 8, "maximum_step": 0.01,     koopman_action: g o Phi^{k dt}
//      "map": {"observables": ["echo", "echo_dot"]} | {"expression": [...]},   the coordinate map h
//      "target": {"lower": [...], "upper": [...], "points": 101},  default: bounds of h(domain)
//      "tolerance": 1e-3, "oversampling": 3,                   push_forward / round_trip conflicts
//                                                              and source refinement
//      "animation": {"frame_seconds": 0.4, "loop": true, "autoplay": true}}}   several frames
//
// pull_back evaluates G o h on the domain, with `function` written in the map's output
// coordinates (an expression whose variables default to the map's output names).
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/CoordinateMaps.hpp>
#include <dynsyslab/core/GridFieldArchive.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <limits>
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
std::vector<std::string> strings(const J& value) {
    std::vector<std::string> out;
    if (std::holds_alternative<std::string>(value.storage())) return {value.string()};
    for (const auto& item : value.array()) out.push_back(item.string());
    return out;
}

// A real or complex function of `variables`, from a catalog, an expression, the library or a K-ROM.
// `pointFunction` is set instead when a K-ROM is a function of coordinates the state cannot give.
struct Function {
    std::optional<c::Observable> observable;
    c::PointFunction pointFunction;
    std::size_t outputs{}, dimension{};
};

Function function(const J& spec, const c::ODESystem& system, const std::vector<std::string>& variables) {
    const auto catalog = system.observableCatalog();
    const auto& parameters = system.configuration().parameters;
    if (!spec.at("observable").isNull()) return {catalog.at(spec.at("observable").string()), {}, 0, 0};
    if (!spec.at("observables").isNull()) return {catalog.coordinates(strings(spec.at("observables")), "domain_function"), {}, 0, 0};
    if (!spec.at("expression").isNull()) {
        auto names = spec.at("variables").isNull() ? variables : strings(spec.at("variables"));
        std::vector<c::AnalyticExpression> parts;
        for (const auto& text : strings(spec.at("expression"))) parts.push_back(c::AnalyticExpression::parse(text, names, parameters));
        const auto label = spec.at("label").string("g");
        return {c::Observable::fromExpressions({label, label, "", "native", {}, false}, std::move(parts)), {}, 0, 0};
    }
    if (!spec.at("library").isNull()) {
        const auto reference = spec.at("reference_state").isNull() ? std::vector<double>{} : numbers(spec.at("reference_state"));
        auto psi = c::libraryEigenfunctions(system, strings(spec.at("library")), reference);
        std::vector<c::Observable> parts(psi.begin(), psi.end());
        return {parts.size() == 1 ? parts.front() : c::Observable::stack("library", "Library eigenfunctions", parts), {}, 0, 0};
    }
    if (!spec.at("krom").isNull()) {
        std::vector<std::size_t> modes;
        for (const auto m : numbers(spec.at("modes"))) modes.push_back(std::size_t(m));
        auto model = std::make_shared<c::KoopmanModel>(
            c::loadKoopmanModel(std::filesystem::u8path(spec.at("krom").string()), system.configuration(), modes));
        if (model->requiresExternalCoordinates())
            return {std::nullopt, [model](const double* q) { return model->evaluateCoordinates(q); },
                    model->eigenfunctionCount(), model->coordinateDimension()};
        std::vector<c::Observable> parts;
        for (std::size_t j = 0; j < model->eigenfunctionCount(); ++j) parts.push_back(model->eigenfunction(j));
        return {parts.size() == 1 ? parts.front() : c::Observable::stack("psi", "K-ROM eigenfunctions", parts), {}, 0, 0};
    }
    throw std::invalid_argument("domain_field function needs observable, observables, expression, library or krom");
}

c::DomainGrid grid(const J& spec, c::FunctionDomain domain, const std::vector<double>& samples) {
    const auto d = domain.dimension;
    const auto points = std::size_t(spec.at("points").number(d <= 2 ? 101 : 41));
    auto out = spec.at("lower").isNull() ? c::DomainGrid::around(domain, samples, points)
                                         : c::DomainGrid::box(domain, numbers(spec.at("lower")), numbers(spec.at("upper")), points);
    std::vector<std::size_t> axes;
    if (!spec.at("axes").isNull()) for (const auto a : numbers(spec.at("axes"))) axes.push_back(std::size_t(a));
    else for (std::size_t j = 0; j < std::min<std::size_t>(d, 3); ++j) axes.push_back(j);
    for (std::size_t j = 0; j < d; ++j) out.anchor[j] = (out.axes[j].lower + out.axes[j].upper) / 2;  // default slice: the middle
    if (!spec.at("anchor").isNull()) out.anchor = numbers(spec.at("anchor"));
    std::vector<c::GridAxis> kept;
    for (const auto a : axes) kept.push_back(out.axes.at(a));
    out.axes = std::move(kept);
    out.validate();
    return out;
}

c::GridField evaluate(const Function& f, const c::DomainGrid& g) {
    return f.observable ? c::evaluateOnGrid(*f.observable, g) : c::evaluateOnGrid(f.pointFunction, f.outputs, g);
}

}  // namespace

std::filesystem::path runDomainField(const std::filesystem::path& sourcePath, const J& settings,
                                     const std::filesystem::path& root, const std::string& label) {
    const auto source = loadSource(sourcePath);
    const auto configuration = systemConfiguration(source);
    const std::shared_ptr<const c::ODESystem> system = c::createSystem(configuration);
    const auto& stateNames = system->descriptor().stateNames;
    const auto operation = settings.at("operation").string("evaluate");
    const auto mapSpec = settings.at("map");
    std::optional<c::Observable> map;
    if (!mapSpec.isNull()) map = function(mapSpec, *system, stateNames).observable;
    if (!map && operation != "evaluate" && operation != "koopman_action")
        throw std::invalid_argument("domain_field " + operation + " needs a coordinate map");

    // The domain: the state space, or the K-ROM's own coordinates when the state cannot give them.
    auto f = operation == "pull_back" ? Function{} : function(settings.at("function"), *system, stateNames);
    const bool external = !f.observable && operation != "pull_back";
    auto domain = external ? c::FunctionDomain::supplied(f.dimension) : c::FunctionDomain::physicalState(source.dimension);
    if (!external) domain.names = stateNames;
    if (external && operation != "evaluate") throw std::invalid_argument("domain_field " + operation + " needs a function of the state");
    const auto g = grid(settings.at("domain"), domain, external ? std::vector<double>{} : source.values);

    c::GridFieldSeries series;
    auto target = [&]() {
        std::vector<double> image;
        const auto pts = g.points();
        for (std::size_t k = 0; k < g.size(); ++k) {
            try {
                const auto v = map->evaluateReal(pts.data() + k * g.dimension());
                image.insert(image.end(), v.begin(), v.end());
            } catch (const std::exception&) {}
        }
        auto q = c::FunctionDomain::supplied(map->outputDimension());
        q.names = map->outputDimension() == 1 ? std::vector<std::string>{map->descriptor().name} : map->descriptor().outputNames;
        if (q.names.size() != q.dimension) q.names.clear();
        return grid(settings.at("target"), q, image);
    };
    const c::PushForwardOptions options{settings.at("tolerance").number(1e-3),
                                        std::size_t(settings.at("oversampling").number(3))};
    if (operation == "evaluate") {
        series.frames = {evaluate(f, g)};
        series.frameTime = {0};
    } else if (operation == "koopman_action") {
        if (configuration.driver.amplitude != 0) throw std::invalid_argument("koopman_action needs the unforced system");
        const double dt = settings.at("time_step").number(0.1);
        const auto steps = std::size_t(settings.at("steps").number(8));
        series.frames = c::koopmanAction(*f.observable, *system, g, dt, steps, settings.at("maximum_step").number(1e-2));
        for (std::size_t k = 0; k <= steps; ++k) series.frameTime.push_back(double(k) * dt);
    } else if (operation == "push_forward") {
        series.frames = {c::pushForward(*f.observable, *map, g, target(), options)};
        series.frameTime = {0};
    } else if (operation == "round_trip") {
        series.frames = {c::roundTripError(*f.observable, *map, g, target(), options)};
        series.frameTime = {0};
    } else if (operation == "fold_set") {
        series.frames = {c::mapVolumeFactor(*map, g)};
        series.frameTime = {0};
    } else if (operation == "pull_back") {
        auto names = map->descriptor().outputNames;
        if (map->outputDimension() == 1) names = {map->descriptor().name};
        const auto G = function(settings.at("function"), *system, names);
        if (!G.observable) throw std::invalid_argument("pull_back needs a function of the map's coordinates");
        series.frames = {c::evaluateOnGrid(c::pullBack(*G.observable, *map), g)};
        series.frameTime = {0};
    } else {
        throw std::invalid_argument("Unknown domain_field operation " + operation);
    }
    auto meta = json::parse(metadata(source, "domain-field", label, "\"operation\":\"" + escaped(operation) + "\""));
    meta.object()["domain_field_settings"] = settings;
    if (!settings.at("animation").isNull()) meta.object()["animation"] = settings.at("animation");  // frame_seconds, loop, autoplay
    series.metadata = std::move(meta);
    return publish(root, source, "domain-field", label, c::gridFieldArchive(series));
}

}  // namespace dynsyslab::workflow::postprocess_detail
