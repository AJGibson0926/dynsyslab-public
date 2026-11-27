// history_control_regression: measurement-only constant input map for a
// causal-history (xi-coordinate) KROM.
//
// Each forced source is a simulator archive with an archived control history on
// the transform clock. Only the named measurement is used: it is pushed through
// the KROM's frozen history transform to xi_k, lifted to psi_k = psi(xi_k), and
// the held-input relation psi_{k+1,j} = mu_j psi_{k,j} + D_j u_k is fitted per
// mode by (ridge) least squares. The continuous map stored for KLQR is the exact
// zero-order-hold inverse B_j = D_j lambda_j / (mu_j - 1) (D_j / dt as lambda->0).
#include "MeasurementMap.hpp"
#include "PostprocessInternal.hpp"
#include "KromEvidence.hpp"

#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runHistoryControlRegression(
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label) {
    namespace c = dynsyslab::core;
    using Complex = std::complex<double>;
    const auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    if (kromPath.empty()) throw std::invalid_argument("history_control_regression requires a KROM");
    if (sources.empty()) throw std::invalid_argument("history_control_regression needs forced sources");
    auto model = c::loadKoopmanModel(kromPath);
    if (!model.coordinates().computableFromMeasurements())
        throw std::invalid_argument("history_control_regression requires a causal-history KROM");
    if (model.inputDimension() != 1)
        throw std::invalid_argument("history_control_regression supports one scalar input");
    const auto measurementMap = loadMeasurementMap(model, kromPath, settings);
    const auto& transformPath = measurementMap.path;
    const auto& transform = measurementMap.transform;
    const double dt = transform.sampleInterval();
    const auto& rates = model.spectrum().continuousRates();
    const std::size_t modes = model.liftedDimension();
    const double ridge = settings.at("ridge").number(0.0);
    std::vector<Complex> numerator(modes), mu(modes);
    for (std::size_t j = 0; j < modes; ++j) mu[j] = std::exp(rates[j] * dt);
    double denominator = ridge, residualEnergy = 0.0, targetEnergy = 0.0;
    std::size_t pairs = 0;
    struct Pair { std::vector<Complex> psi, next; double u; };
    std::vector<Pair> rows;
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        if (source.inputs().empty() || source.inputDimension() != 1)
            throw std::invalid_argument("history_control_regression sources need a scalar control history");
        if (std::abs(uniformDeltaTime(source) - dt) > 1e-8 * dt)
            throw std::invalid_argument("forced source must be sampled on the history transform clock");
        const auto system = c::createSystem(systemConfiguration(source));
        const auto name = settings.at("measurement_observable").string(transform.measurement());
        const auto measurement = system->observableCatalog().at(name);
        c::CausalHistoryProjector projector(transform);
        std::optional<std::vector<Complex>> previous; double previousInput = 0.0;
        const std::size_t inputRows = source.inputs().size() / source.inputDimension();
        for (std::size_t k = 0; k < source.samples(); ++k) {
            const auto y = measurement.evaluateReal(source.values().data() + k * source.dimension()).front();
            const auto q = projector.push(source.time()[k], y);
            if (!q) continue;
            auto psi = model.evaluateCoordinates(q->data());
            if (previous) rows.push_back({*previous, psi, previousInput});
            previous = std::move(psi);
            previousInput = k < inputRows ? source.inputs()[k] : 0.0;
        }
    }
    for (const auto& row : rows) {
        denominator += row.u * row.u;
        for (std::size_t j = 0; j < modes; ++j) numerator[j] += row.u * (row.next[j] - mu[j] * row.psi[j]);
        ++pairs;
    }
    if (pairs < 3 || !(denominator > 0)) throw std::invalid_argument("history_control_regression has no excited pairs");
    c::ConstantControlMap map;
    map.modes = modes; map.inputs = 1; map.matrix.resize(modes);
    map.provenance = "history multistep held-input regression (measurement only), ZOH-inverted";
    map.discreteMatrix.resize(modes); map.sampleTime = dt;
    for (std::size_t j = 0; j < modes; ++j) {
        const Complex D = numerator[j] / denominator;
        map.discreteMatrix[j] = D;
        const Complex z = rates[j] * dt;
        map.matrix[j] = std::abs(z) < 1e-12 ? D / dt : D * rates[j] / (mu[j] - 1.0);
    }
    for (const auto& row : rows)
        for (std::size_t j = 0; j < modes; ++j) {
            const Complex target = row.next[j] - mu[j] * row.psi[j];
            residualEnergy += std::norm(target - map.discreteMatrix[j] * row.u);
            targetEnergy += std::norm(target);
        }
    const auto representation = settings.at("representation").string("constant");
    model.setControlMap(representation, map);
    json::Value::Array sourceList; for (const auto& path : sources) sourceList.emplace_back(path.generic_u8string());
    json::Value::Object meta{
        {"schema", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("history-control-regression")},
        {"label", json::Value(label)}, {"source_krom", json::Value(kromPath.generic_u8string())},
        {"source_archive", json::Value(kromPath.generic_u8string())},
        {"forced_sources", json::Value(std::move(sourceList))},
        {"history_transform", json::Value(transformPath)}, {"representation", json::Value(representation)},
        {"regression_pairs", json::Value(double(pairs))}, {"ridge", json::Value(ridge)},
        {"relative_residual", json::Value(targetEnergy > 0 ? std::sqrt(residualEnergy / targetEnergy) : 0.0)},
        {"control_map_time_domain", json::Value("continuous-time")}, {"delta_t", json::Value(dt)},
        {"independent_archive", json::Value(true)}};
    const auto directory = kindFolder(root, dynsyslab::core::object_kind::krom);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string()
        + "__history-control-regression__" + timestamp()) + ".npz");
    // A child of the krom: its modes and the input map B fitted here (the parent keeps its training samples).
    c::saveKoopmanModel(path, model, withDataObject(json::dump(json::Value(std::move(meta))), dynsyslab::core::object_kind::krom, {kromPath}));
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
