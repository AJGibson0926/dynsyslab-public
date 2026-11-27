// Post-processing method ssm: a spectral submanifold and its reduced dynamics
// (core/SpectralSubmanifold.hpp), published as a dynamical_model object (the reduced vector field in
// real chart coordinates, the dictionary-and-coefficients form every polynomial model uses), with
// its backbone curve as a one-axis field object Viewer draws as line plots.
//
//   {"name": "ssm", "settings": {
//      "mode": "parameterization" (default: from the vector field) | "data" (from the sources),
//      parameterization:
//        "vector_field": {"expressions": ["v", "-x - x^3"], "variables": ["x", "v"]},
//                        default: the source system's autonomous drift
//        "equilibrium": [...]        default: Newton's method from the first source's last sample
//        "master_modes": [0], "order": 5, "style": "normal_form" | "graph",
//        "resonance_tolerance": 1e-6, "imaginary_resonance": true
//      data (every source is one trajectory, its components the observables):
//        "dimension": 2, "order": 3, "dynamics_order": 3, "normal_form_order": 0, "ridge": 0,
//        "skip_time": 0, "equilibrium": [...], "derivative_method": "central"
//      "backbone": {"radius": 0 (from the data), "points": 101, "coordinates": ["position"]}}}
//
// The archive holds dynamics_model_json and coefficients (reduced dynamics), ssm_parameterization
// (features x state over polynomialExponentTable(d, order, true): x as a polynomial in the chart),
// eigenvalues, the spectral quotient and invariance or fit errors; the data mode adds the tangent
// space, the singular values and each source's model prediction (prediction_<k>_values).
#include "GeometryPostprocess.hpp"
#include "ObjectPublishing.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/GridFieldArchive.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/SpectralSubmanifold.hpp>

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

J complexArray(const std::vector<std::complex<double>>& values) {
    J::Array out;
    for (const auto& z : values) out.push_back(J::Array{z.real(), z.imag()});
    return out;
}

// The reduced vector field as a polynomial dynamical model (dictionary and coefficients).
c::SindyDynamicsModel polynomialModel(std::size_t d, std::size_t degree, std::vector<double> coefficients) {
    const c::PolynomialDictionary dictionary{d, c::polynomialExponentTable(d, degree, true), std::vector<double>(d, 0),
                                             std::vector<double>(d, 1)};
    return {dictionary.asObservable(), std::move(coefficients)};
}

std::vector<std::size_t> coordinateIndices(const J& given, const std::vector<std::string>& names, std::size_t n) {
    std::vector<std::size_t> out;
    if (given.isNull()) { for (std::size_t i = 0; i < n; ++i) out.push_back(i); return out; }
    for (const auto& item : given.array()) {
        if (std::holds_alternative<double>(item.storage())) { out.push_back(std::size_t(item.number())); continue; }
        const auto at = std::find(names.begin(), names.end(), item.string());
        if (at == names.end()) throw std::invalid_argument("ssm backbone: unknown coordinate " + item.string());
        out.push_back(std::size_t(at - names.begin()));
    }
    for (const auto i : out) if (i >= n) throw std::invalid_argument("ssm backbone: coordinate index out of range");
    return out;
}

// Backbone of a two-dimensional normal-form SSM as a field over rho (empty path when not applicable).
std::filesystem::path publishBackbone(const std::filesystem::path& root, const SourceData& source, const c::SpectralSubmanifold& ssm,
                                      const J& spec, double radius, const c::StateOutputs& outputs,
                                      const std::vector<std::string>& names, J& summary) {
    if (ssm.dimension != 2 || std::abs(ssm.chart[1].imag()) != 1) return {};
    const auto points = std::max<std::size_t>(2, std::size_t(spec.at("points").number(101)));
    radius = spec.at("radius").number(radius);
    if (!(radius > 0)) return {};
    c::DomainGrid grid;
    grid.domain = c::FunctionDomain::supplied(1);
    grid.domain.names = {"rho"};
    grid.anchor = {0};
    grid.axes = {c::GridAxis{0, 0, radius, points}};
    grid.validate();
    c::BackboneCurve curve;
    try {
        curve = c::backboneCurve(ssm, grid.axisValues(0), outputs, names.size());
    } catch (const std::invalid_argument&) {
        return {};  // graph-style reduced dynamics have no backbone in this form
    }
    c::GridField field;
    field.grid = grid;
    field.outputs = 2 + names.size();
    field.outputNames = {"frequency", "damping_ratio"};
    for (const auto& name : names) field.outputNames.push_back("amplitude_" + name);
    field.status.assign(points, c::GridStatus::valid);
    for (std::size_t k = 0; k < points; ++k) {
        field.values.push_back(curve.frequency[k]);
        field.values.push_back(curve.damping[k]);
        for (std::size_t r = 0; r < names.size(); ++r) field.values.push_back(curve.amplitude[k * names.size() + r]);
        if (!std::isfinite(curve.damping[k])) field.status[k] = c::GridStatus::undefined;
    }
    c::GridFieldSeries series;
    series.frames = {field};
    series.frameTime = {0};
    auto meta = json::parse(metadata(source, "ssm-backbone", "SSM backbone curve", ""));
    meta.object()["backbone"] = J::Object{{"radius", radius}, {"points", double(points)}, {"coordinates", [&] {
        J::Array a; for (const auto& n : names) a.push_back(n); return a; }()}};
    series.metadata = std::move(meta);
    summary = J::Object{{"radius", radius}, {"frequency_at_zero", curve.frequency.front()}, {"frequency_at_radius", curve.frequency.back()}};
    return publishObject(root, source, c::object_kind::field, "ssm-backbone", c::gridFieldArchive(series));
}

void ssmArrays(c::NpzRecord& archive, J::Object& meta, const c::SpectralSubmanifold& ssm, const std::string& prefix) {
    const auto features = c::polynomialExponentTable(ssm.dimension, ssm.order, true).size() / ssm.dimension;
    archive.arrays[prefix + "parameterization"] = c::NpyArray::fromDoubles({features, ssm.stateDimension}, ssm.realParameterization());
    std::vector<double> re, im;
    for (const auto& z : ssm.eigenvalues) { re.push_back(z.real()); im.push_back(z.imag()); }
    archive.arrays[prefix + "eigenvalues_real"] = c::NpyArray::fromDoubles({re.size()}, re);
    archive.arrays[prefix + "eigenvalues_imag"] = c::NpyArray::fromDoubles({im.size()}, im);
    J::Array masters, warnings;
    std::vector<std::complex<double>> masterValues;
    for (const auto m : ssm.masterModes) { masters.push_back(double(m)); masterValues.push_back(ssm.eigenvalues[m]); }
    for (const auto& w : ssm.warnings) warnings.push_back(w);
    meta[prefix + "master_modes"] = masters;
    meta[prefix + "master_eigenvalues"] = complexArray(masterValues);
    meta[prefix + "spectral_quotient"] = std::isfinite(ssm.spectralQuotient) ? J(ssm.spectralQuotient) : J();
    meta[prefix + "warnings"] = warnings;
    meta[prefix + "order"] = double(ssm.order);
    meta[prefix + "dimension"] = double(ssm.dimension);
}

}  // namespace

std::filesystem::path runSsm(const std::vector<std::filesystem::path>& sources, const J& settings, const std::filesystem::path& root,
                             const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("ssm needs a source");
    const auto mode = settings.at("mode").string("parameterization");
    const auto first = loadSource(sources.front());
    const std::shared_ptr<const c::ODESystem> system = c::createSystem(systemConfiguration(first));
    c::NpzRecord archive;
    auto meta = json::parse(metadata(first, "ssm", label, "\"schema\":\"dynsyslab.spectral-submanifold/v1\"")).object();
    meta["settings"] = settings;
    meta["mode"] = mode;
    meta["model_type"] = "autonomous-vector-field";
    meta["references"] = J::Array{"Haller and Ponsioen 2016, Nonlinear Dyn. 86:1493-1534",
        "Jain and Haller 2022, Nonlinear Dyn. 107:1417-1450",
        "Cenedese, Axas, Bauerlein, Avila and Haller 2022, Nat. Commun. 13:872"};
    std::filesystem::path backbone;
    J backboneSummary;
    const auto& spec = settings.at("backbone");
    std::vector<std::string> stateNames = system->descriptor().stateNames;
    if (mode == "parameterization") {
        std::vector<c::AnalyticExpression> field;
        const auto& given = settings.at("vector_field");
        if (!given.isNull()) {
            std::vector<std::string> variables = stateNames;
            if (!given.at("variables").isNull()) { variables.clear(); for (const auto& v : given.at("variables").array()) variables.push_back(v.string()); }
            for (const auto& e : given.at("expressions").array())
                field.push_back(c::AnalyticExpression::parse(e.string(), variables, system->configuration().parameters));
            stateNames = variables;
        } else {
            field = system->driftExpressions();
            if (field.empty()) throw std::invalid_argument("ssm: the source system has no analytic drift; give settings.vector_field");
        }
        const auto n = field.size();
        if (stateNames.size() != n) { stateNames.clear(); for (std::size_t i = 0; i < n; ++i) stateNames.push_back("x_" + std::to_string(i + 1)); }
        std::vector<double> guess(n, 0);
        if (first.samples() && first.dimension() == n)
            guess.assign(first.values().end() - std::ptrdiff_t(n), first.values().end());
        const auto x0 = settings.at("equilibrium").isNull() ? c::equilibriumNear(field, guess) : numbers(settings.at("equilibrium"));
        c::SpectralSubmanifoldOptions options;
        options.order = std::size_t(settings.at("order").number(5));
        if (!settings.at("master_modes").isNull()) for (const auto m : numbers(settings.at("master_modes"))) options.masterModes.push_back(std::size_t(m));
        const auto style = settings.at("style").string("normal_form");
        if (style != "normal_form" && style != "graph") throw std::invalid_argument("ssm style is normal_form or graph");
        options.style = style == "graph" ? c::SsmStyle::graph : c::SsmStyle::normalForm;
        options.resonanceTolerance = settings.at("resonance_tolerance").number(1e-6);
        options.imaginaryResonance = settings.at("imaginary_resonance").boolean(true);
        const auto ssm = c::spectralSubmanifold(field, x0, options);
        const auto model = polynomialModel(ssm.dimension, ssm.order, ssm.realReducedDynamics());
        archive.arrays["dynamics_model_json"] = c::NpyArray::fromString(json::dump(model.definition()));
        archive.arrays["coefficients"] = c::NpyArray::fromDoubles({model.coefficients().size() / ssm.dimension, ssm.dimension}, model.coefficients());
        ssmArrays(archive, meta, ssm, "ssm_");
        archive.arrays["equilibrium"] = c::NpyArray::fromDoubles({n}, x0);
        // Radius: the linear part reaches the source's largest deviation from the equilibrium.
        double deviation = 0, linear = 0;
        if (first.dimension() == n)
            for (std::size_t k = 0; k < first.samples(); ++k)
                for (std::size_t i = 0; i < n; ++i) deviation = std::max(deviation, std::abs(first.values()[k * n + i] - x0[i]));
        for (std::size_t i = 0; i < n; ++i) linear = std::max(linear, 2 * std::abs(ssm.W[n + i]));
        const double radius = linear > 0 && deviation > 0 ? deviation / linear : 0.1;
        meta["ssm_invariance_error"] = ssm.invarianceError(field, radius);
        meta["ssm_invariance_radius"] = radius;
        const auto indices = coordinateIndices(spec.at("coordinates"), stateNames, n);
        std::vector<std::string> names;
        for (const auto i : indices) names.push_back(stateNames[i]);
        backbone = publishBackbone(root, first, ssm, spec, radius, [indices](const double* x) {
            std::vector<double> out; for (const auto i : indices) out.push_back(x[i]); return out; }, names, backboneSummary);
        J::Array coordinates;
        for (std::size_t j = 0; j < ssm.dimension; ++j) coordinates.push_back("q_" + std::to_string(j + 1));
        meta["coordinate_names"] = coordinates;
    } else if (mode == "data") {
        std::vector<c::TimeSeries> trajectories;
        for (const auto& path : sources) trajectories.push_back(loadSource(path).record.series);
        c::SsmFitOptions options;
        options.dimension = std::size_t(settings.at("dimension").number(2));
        options.order = std::size_t(settings.at("order").number(3));
        options.dynamicsOrder = std::size_t(settings.at("dynamics_order").number(3));
        options.normalFormOrder = std::size_t(settings.at("normal_form_order").number(0));
        options.ridge = settings.at("ridge").number(0);
        options.skipTime = settings.at("skip_time").number(0);
        if (!settings.at("equilibrium").isNull()) options.equilibrium = numbers(settings.at("equilibrium"));
        options.derivative.method = settings.at("derivative_method").string("central");
        const auto fit = c::fitSpectralSubmanifold(trajectories, options);
        const auto n = fit.observableDimension, d = fit.dimension;
        const auto model = polynomialModel(d, fit.dynamicsOrder, fit.dynamics);
        archive.arrays["dynamics_model_json"] = c::NpyArray::fromString(json::dump(model.definition()));
        archive.arrays["coefficients"] = c::NpyArray::fromDoubles({fit.dynamics.size() / d, d}, fit.dynamics);
        archive.arrays["ssm_parameterization"] = c::NpyArray::fromDoubles({fit.manifold.size() / n, n}, fit.manifold);
        archive.arrays["ssm_tangent"] = c::NpyArray::fromDoubles({n, d}, fit.tangent);
        archive.arrays["singular_values"] = c::NpyArray::fromDoubles({fit.singularValues.size()}, fit.singularValues);
        archive.arrays["equilibrium"] = c::NpyArray::fromDoubles({n}, fit.equilibrium);
        ssmArrays(archive, meta, fit.normalForm, "normal_form_");
        meta["manifold_order"] = double(fit.order);
        meta["dynamics_order"] = double(fit.dynamicsOrder);
        meta["manifold_relative_error"] = fit.manifoldError;
        meta["dynamics_relative_error"] = fit.dynamicsError;
        // Each source predicted from its first sample past the transient (normalized trajectory error).
        J::Array predictions;
        double radius = 0;
        for (std::size_t k = 0; k < trajectories.size(); ++k) {
            const auto& t = trajectories[k];
            const auto test = t.cropTime(t.time(0) + options.skipTime, t.time(t.sampleCount() - 1));
            const auto predicted = fit.predict(test);
            double error = 0, size = 0;
            for (std::size_t s = 0; s < test.sampleCount(); ++s)
                for (std::size_t i = 0; i < n; ++i) {
                    error += std::pow(predicted.value(s, i) - test.value(s, i), 2);
                    size += std::pow(test.value(s, i) - fit.equilibrium[i], 2);
                    const auto eta = fit.reduce(test.valuesData() + s * n);
                    for (const auto e : eta) radius = std::max(radius, std::abs(e));
                }
            archive.arrays["prediction_" + std::to_string(k) + "_values"] = c::NpyArray::fromDoubles({test.sampleCount(), n}, predicted.valueVector());
            archive.arrays["prediction_" + std::to_string(k) + "_time"] = c::NpyArray::fromDoubles({test.sampleCount()}, predicted.timeVector());
            predictions.push_back(J::Object{{"source", sources[k].generic_u8string()},
                                            {"relative_error", size > 0 ? std::sqrt(error / size) : 0.0}});
        }
        meta["predictions"] = predictions;
        std::vector<std::string> names = stateNames.size() == n ? stateNames : std::vector<std::string>{};
        if (names.empty()) for (std::size_t i = 0; i < n; ++i) names.push_back("y_" + std::to_string(i + 1));
        const auto indices = coordinateIndices(spec.at("coordinates"), names, n);
        std::vector<std::string> chosen;
        for (const auto i : indices) chosen.push_back(names[i]);
        // Normal-form radius reaching the data: eta ~ 2 Re(p v), |v| = 1.
        const auto& nf = fit.normalForm;
        double linear = 0;
        if (nf.dimension) for (std::size_t i = 0; i < d; ++i) linear = std::max(linear, 2 * std::abs(nf.W[d + i]));
        backbone = publishBackbone(root, first, nf, spec, linear > 0 ? radius / linear : 0, [&fit, indices](const double* eta) {
            const auto y = fit.lift(eta);
            std::vector<double> out; for (const auto i : indices) out.push_back(y[i]); return out; }, chosen, backboneSummary);
        J::Array coordinates;
        for (std::size_t j = 0; j < d; ++j) coordinates.push_back("eta_" + std::to_string(j + 1));
        meta["coordinate_names"] = coordinates;
    } else {
        throw std::invalid_argument("ssm mode is parameterization or data");
    }
    if (!backbone.empty()) { meta["backbone_archive"] = backbone.generic_u8string(); meta["backbone"] = backboneSummary; }
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    const auto path = publishObject(root, first, c::object_kind::dynamicalModel, "ssm", std::move(archive));
    if (!backbone.empty()) registerAdditionalObject(backbone);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
