// input_map: the KROM control matrix B_psi as one switchable post-processing step.
//
// The eigenfunctions psi come from autonomous (unforced) training and are never
// refitted here. This step only decides how the input enters psi' = Lambda psi + B_psi u:
//
//   representation "exact"          analytic grad(psi) . g(x) (state/observable KROMs)
//   representation "history_exact"  analytic grad_xi(psi) . newest-sample weight . L_g(y)
//                                   (causal-history KROMs; controller reads the simulator state)
//   representation "constant"       B_psi = const, regressed from FORCED records
//   representation "bilinear"       B_psi(psi) = b0 + sum_i b_i psi_i, regressed from FORCED records
//
// Regression ("fit"):
//   "held_input"           any KROM: psi_k from the state, measured coordinates or frozen history
//                          transform, q_k = psi_{k+1} - mu psi_k = [a(psi_k)] + D(psi_k) u_k per mode.
//                          The eigenvalues stay frozen. ZOH inversion, first order in u: an offset
//                          divides by phi(lambda_j, 0), an operator on psi_i by phi(lambda_j, lambda_i),
//                          phi(a, b) = int_0^dt exp(a (dt - s)) exp(b s) ds (Van Loan, IEEE TAC 23
//                          (1978) 395-404). "conjugate_complete" adds conj(psi_i) regressors.
//                          "drift": "affine" fits a nuisance a(psi) = a0 + sum a_i psi_i jointly
//                          (it absorbs zero-input closure error; it is not used by the controller).
//   "controlled_residual"  state/observable KROMs: existing native controlled-residual regression.
//   "multistep_output"     state/observable KROMs, constant only: existing multistep output regression.
//   "exact_map"            least-squares fit of the constant/bilinear family to the exact map
//                          sampled on the sources' states (state/observable KROMs).
//   measured-coordinate KROMs (coordinate_transform measured_derivatives): "controlled_residual" on
//                          the forced records' measured coordinates. psi from the coordinates,
//                          psi' from the "derivative" differentiator (default central), and the input
//                          passed through the same smoother so both sides see the same filtering.
//
// The output is a KROM archive with the map stored under "map_name" (default: the
// representation) so a controller selects it with control_representation alone.
#include "MeasurementMap.hpp"
#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"
#include "DerivativeSettings.hpp"
#include "HeldInputRecords.hpp"

#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/ControlRegression.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using Complex = std::complex<double>;

// Solve (X^H X + ridge I) beta = X^H y for small complex systems (row-major X: rows x n).
std::vector<Complex> leastSquares(const std::vector<Complex>& X, const std::vector<Complex>& y,
                                  std::size_t rows, std::size_t n, double ridge) {
    std::vector<Complex> A(n * n), b(n);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t i = 0; i < n; ++i) {
            const auto xi = std::conj(X[r * n + i]);
            b[i] += xi * y[r];
            for (std::size_t j = 0; j < n; ++j) A[i * n + j] += xi * X[r * n + j];
        }
    double scale = 0; for (std::size_t i = 0; i < n; ++i) scale = std::max(scale, std::abs(A[i * n + i]));
    for (std::size_t i = 0; i < n; ++i) A[i * n + i] += ridge * (scale > 0 ? scale : 1.0);
    for (std::size_t k = 0; k < n; ++k) {  // Gaussian elimination with partial pivoting
        std::size_t p = k;
        for (std::size_t i = k + 1; i < n; ++i) if (std::abs(A[i * n + k]) > std::abs(A[p * n + k])) p = i;
        if (!(std::abs(A[p * n + k]) > 1e-300)) throw std::invalid_argument("input_map regression is singular (no excitation?)");
        if (p != k) { for (std::size_t j = 0; j < n; ++j) std::swap(A[k * n + j], A[p * n + j]); std::swap(b[k], b[p]); }
        for (std::size_t i = k + 1; i < n; ++i) {
            const auto f = A[i * n + k] / A[k * n + k];
            for (std::size_t j = k; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }
    std::vector<Complex> x(n);
    for (std::size_t k = n; k-- > 0;) {
        Complex s = b[k];
        for (std::size_t j = k + 1; j < n; ++j) s -= A[k * n + j] * x[j];
        x[k] = s / A[k * n + k];
    }
    return x;
}

// phi(a, b) = int_0^h exp(a (h - s)) exp(b s) ds, the held-input weight of a term that
// evolves at rate b inside a mode that evolves at rate a.
Complex heldKernel(Complex a, Complex b, double h) {
    const Complex d = a - b;
    if (std::abs(d) * h < 1e-6) return h * std::exp(b * h) * (1.0 + d * h / 2.0 + d * d * h * h / 6.0);
    return (std::exp(a * h) - std::exp(b * h)) / d;
}

Complex inverseKernel(Complex a, Complex b, double h) {
    const Complex phi = heldKernel(a, b, h);
    if (!(std::abs(phi) > 1e-12 * h))
        throw std::invalid_argument("input_map held-input conversion is aliased at this sample time");
    return 1.0 / phi;
}

json::Value::Object heldInputRegression(c::KoopmanModel& model, const HeldInputRecords& records,
                                        const json::Value& settings, const std::string& representation,
                                        const std::string& mapName) {
    const auto& rows = records.rows;
    if (rows.size() < 8) throw std::invalid_argument("input_map has too few forced regression pairs");
    const double dt = records.sampleTime, ridge = settings.at("ridge").number(0.0);
    const auto drift = settings.at("drift").string("none");
    if (drift != "none" && drift != "affine") throw std::invalid_argument("input_map drift must be none or affine");
    const auto& rates = model.spectrum().continuousRates();
    const std::size_t m = model.liftedDimension();
    const bool bilinear = representation == "bilinear";
    const bool conjugate = bilinear && settings.at("conjugate_complete").boolean(false);
    // conj(psi_i) duplicates psi_i for a real mode, so only complex modes get a column.
    std::vector<std::size_t> conjugateModes;
    if (conjugate) for (std::size_t i = 0; i < m; ++i) if (std::abs(rates[i].imag()) > 1e-12) conjugateModes.push_back(i);
    const std::size_t nc = conjugateModes.size();
    const std::size_t nb = bilinear ? 1 + m + nc : 1, nd = drift == "affine" ? m + 1 : 0, n = nb + nd;
    std::vector<Complex> X(rows.size() * n);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        auto* x = X.data() + r * n;
        x[0] = rows[r].u;
        if (bilinear) for (std::size_t i = 0; i < m; ++i) x[1 + i] = rows[r].u * rows[r].psi[i];
        for (std::size_t c = 0; c < nc; ++c) x[1 + m + c] = rows[r].u * std::conj(rows[r].psi[conjugateModes[c]]);
        if (nd) { x[nb] = 1.0; for (std::size_t i = 0; i < m; ++i) x[nb + 1 + i] = rows[r].psi[i]; }
    }
    c::ConstantControlMap constant; c::BilinearControlMap bil;
    constant.modes = bil.modes = m; constant.inputs = bil.inputs = 1;
    constant.matrix.resize(m); constant.discreteMatrix.resize(m); constant.sampleTime = dt;
    bil.offsets.resize(m); bil.operators.assign(m * m, 0.0);
    if (conjugate) bil.conjugateOperators.assign(m * m, 0.0);
    double residual = 0, energy = 0, inputEnergy = 0;
    for (const auto& row : rows) inputEnergy += row.u * row.u;
    if (!(inputEnergy > 0)) throw std::invalid_argument("input_map forced records carry no input");
    json::Value::Array perMode;
    for (std::size_t j = 0; j < m; ++j) {
        const Complex mu = std::exp(rates[j] * dt);
        std::vector<Complex> q(rows.size());
        for (std::size_t r = 0; r < rows.size(); ++r) q[r] = rows[r].next[j] - mu * rows[r].psi[j];
        const auto beta = leastSquares(X, q, rows.size(), n, ridge);
        double res = 0, en = 0;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            Complex fit = 0; for (std::size_t i = 0; i < n; ++i) fit += X[r * n + i] * beta[i];
            res += std::norm(q[r] - fit); en += std::norm(q[r]);
        }
        residual += res; energy += en;
        const Complex toContinuous = inverseKernel(rates[j], 0.0, dt);
        if (bilinear) {
            bil.offsets[j] = beta[0] * toContinuous;
            for (std::size_t i = 0; i < m; ++i)
                bil.operators[j * m + i] = beta[1 + i] * inverseKernel(rates[j], rates[i], dt);
            for (std::size_t c = 0; c < nc; ++c) {
                const auto i = conjugateModes[c];
                bil.conjugateOperators[j * m + i] = beta[1 + m + c] * inverseKernel(rates[j], std::conj(rates[i]), dt);
            }
        } else { constant.discreteMatrix[j] = beta[0]; constant.matrix[j] = beta[0] * toContinuous; }
        perMode.emplace_back(en > 0 ? std::sqrt(res / en) : 0.0);
    }
    const std::string provenance = std::string("input_map held_input ") + representation
        + (conjugate ? " conjugate-complete" : "") + (nd ? " with affine drift nuisance" : "")
        + " (forced records, frozen eigenvalues), ZOH-inverted to first order in u";
    if (bilinear) { bil.provenance = provenance; model.setControlMap(mapName, bil); }
    else { constant.provenance = provenance; model.setControlMap(mapName, constant); }
    return {{"regression_pairs", json::Value(double(rows.size()))}, {"ridge", json::Value(ridge)},
            {"drift", json::Value(drift)}, {"relative_residual", json::Value(energy > 0 ? std::sqrt(residual / energy) : 0.0)},
            {"relative_residual_by_mode", json::Value(std::move(perMode))}, {"delta_t", json::Value(dt)},
            {"conjugate_complete", json::Value(conjugate)},
            {"pairs_dropped_at_input_change", json::Value(double(records.droppedAtInputChange))},
            {"input_rms", json::Value(std::sqrt(inputEnergy / double(rows.size())))}};
}

json::Value::Object coordinateRegression(c::KoopmanModel& model, const std::vector<std::filesystem::path>& sources,
                                         const json::Value& settings, const std::string& representation,
                                         const std::string& fit, const std::string& mapName) {
    if (representation != "constant" || fit != "controlled_residual")
        throw std::invalid_argument("measured-coordinate KROMs support constant with controlled_residual, or constant/bilinear with held_input");
    const auto derivative = parseDerivativeSpec(settings.at("derivative"));
    const double ridge = settings.at("ridge").number(0.0);
    const auto rates = model.continuousEigenvalues();
    const std::size_t m = rates.size(), q = model.coordinateDimension();
    std::vector<Complex> numerator(m, 0.0); double uu = 0; std::vector<double> rr(m, 0.0); std::size_t rows = 0;
    std::vector<std::vector<Complex>> residuals(m); std::vector<double> inputsUsed;
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        if (source.inputs.empty() || datasetDynamics(source) != "forced" || source.inputDimension != 1)
            throw std::invalid_argument("measured-coordinate input_map needs forced coordinate records with one input: " + path.generic_u8string());
        if (source.dimension != q) throw std::invalid_argument("forced record coordinates differ from the KROM coordinates");
        const std::size_t n = source.samples; const double dt = source.time[1] - source.time[0];
        std::vector<std::vector<double>> re(m, std::vector<double>(n)), im(m, std::vector<double>(n));
        for (std::size_t i = 0; i < n; ++i) {
            const auto psi = model.evaluateCoordinates(&source.values[i * q]);
            for (std::size_t j = 0; j < m; ++j) { re[j][i] = psi[j].real(); im[j][i] = psi[j].imag(); }
        }
        const auto u = c::differentiate(source.inputs, dt, derivative).values;  // same smoothing as psi
        const std::size_t skip = std::max<std::size_t>(2, n / 50);
        for (std::size_t j = 0; j < m; ++j) {
            const auto dre = c::differentiate(re[j], dt, derivative), dim = c::differentiate(im[j], dt, derivative);
            for (std::size_t i = skip; i + skip < n; ++i) {
                const Complex psi(dre.values[i], dim.values[i]), dpsi(dre.derivative[i], dim.derivative[i]);
                const Complex r = dpsi - rates[j] * psi;
                numerator[j] += r * u[i]; rr[j] += std::norm(r);
                residuals[j].push_back(r); if (j == 0) inputsUsed.push_back(u[i]);
            }
        }
        for (std::size_t i = skip; i + skip < n; ++i) uu += u[i] * u[i];
        rows += n - 2 * skip;
    }
    c::ConstantControlMap map{{}, model.liftedDimension(), 1, "measured-coordinate controlled residual"};
    map.matrix.resize(m);
    double residual = 0, energy = 0; json::Value::Array perMode;
    for (std::size_t j = 0; j < m; ++j) {
        map.matrix[j] = numerator[j] / (uu + ridge);
        double e = 0; for (std::size_t k = 0; k < residuals[j].size(); ++k) e += std::norm(residuals[j][k] - map.matrix[j] * inputsUsed[k]);
        residual += e; energy += rr[j]; perMode.emplace_back(rr[j] > 0 ? std::sqrt(e / rr[j]) : 0.0);
    }
    map.provenance = "input_map controlled_residual on measured coordinates (derivative " + derivative.method + "), forced records";
    model.setControlMap(mapName, map);
    return {{"relative_residual", json::Value(energy > 0 ? std::sqrt(residual / energy) : 0.0)}, {"ridge", json::Value(ridge)},
            {"regression_rows", json::Value(double(rows))}, {"derivative_method", json::Value(derivative.method)},
            {"relative_residual_by_mode", json::Value(std::move(perMode))}};
}

json::Value::Object stateRegression(c::KoopmanModel& model, const std::vector<std::filesystem::path>& sources,
                                    const json::Value& settings, const std::string& representation,
                                    const std::string& fit, const std::string& mapName) {
    std::vector<SourceData> loaded; std::vector<c::ControlledData> data; std::vector<double> states;
    for (const auto& path : sources) {
        auto source = loadSource(path);
        if (fit != "exact_map" && (source.inputs.empty() || datasetDynamics(source) != "forced"))
            throw std::invalid_argument("input_map regression needs forced records with a control history: " + path.generic_u8string());
        c::requireCompatibleSystems(model.systemConfiguration(), systemConfiguration(source));
        states.insert(states.end(), source.values.begin(), source.values.end());
        loaded.push_back(std::move(source));
    }
    const auto& window = settings.at("state_window");
    const double fitStart = settings.at("fit_start").number(0.0);
    const bool wholeRecords = window.isNull() && std::all_of(loaded.begin(), loaded.end(),
        [&](const SourceData& s) { return s.time.front() >= fitStart; });
    if (wholeRecords) {
        for (const auto& s : loaded) data.push_back({s.time, s.values, s.inputs, s.samples, s.dimension, s.inputDimension});
    } else {
        // Restrict the regression to samples at or after fit_start and inside a state box ([[lo, hi], ...]
        // per state component): each maximal run of inside samples is its own record, so no pair spans a gap.
        states.clear();
        for (const auto& s : loaded) {
            if (!window.isNull() && window.array().size() != s.dimension)
                throw std::invalid_argument("input_map state_window needs one [low, high] pair per state component");
            const bool intervals = s.inputs.size() == (s.samples - 1) * s.inputDimension;
            auto inside = [&](std::size_t k) {
                if (s.time[k] < fitStart) return false;
                if (window.isNull()) return true;
                for (std::size_t j = 0; j < s.dimension; ++j) {
                    const auto& b = window.array().at(j).array(); const double v = s.values[k * s.dimension + j];
                    if (v < b.at(0).number() || v > b.at(1).number()) return false;
                }
                return true;
            };
            for (std::size_t k = 0; k < s.samples;) {
                if (!inside(k)) { ++k; continue; }
                std::size_t e = k; while (e + 1 < s.samples && inside(e + 1)) ++e;
                if (e > k) {
                    c::ControlledData d; d.samples = e - k + 1; d.stateDimension = s.dimension; d.inputDimension = s.inputDimension;
                    d.time.assign(s.time.begin() + k, s.time.begin() + e + 1);
                    d.states.assign(s.values.begin() + k * s.dimension, s.values.begin() + (e + 1) * s.dimension);
                    if (!s.inputs.empty()) d.inputs.assign(s.inputs.begin() + k * s.inputDimension,
                                    s.inputs.begin() + (intervals ? e : e + 1) * s.inputDimension);
                    states.insert(states.end(), d.states.begin(), d.states.end());
                    data.push_back(std::move(d));
                }
                k = e + 1;
            }
        }
        if (data.empty()) throw std::invalid_argument("input_map fit_start/state_window keeps no forced samples");
    }
    const double ridge = settings.at("ridge").number(0.0);
    const std::size_t dimension = loaded.front().dimension, count = states.size() / dimension;
    double relative = 0;
    if (representation == "constant") {
        auto fitted = [&] {
            if (fit == "exact_map") return c::fitConstantExactControlMap(model, states, count);
            if (fit == "controlled_residual") return c::fitConstantControlResidual(model, data, ridge);
            if (fit != "multistep_output") throw std::invalid_argument("constant input_map fit must be multistep_output, controlled_residual or exact_map");
            const auto reconstruction = settings.at("reconstruction").string("state");
            auto system = c::createSystem(model.systemConfiguration());
            return c::fitConstantControlMultistepOutput(model, data, system->observableCatalog().at(reconstruction),
                                                        model.reconstruction(reconstruction), ridge);
        }();
        relative = fitted.diagnostics.relativeResidual; model.setControlMap(mapName, fitted.representation);
    } else {
        const bool conjugate = settings.at("conjugate_complete").boolean(false);
        auto fitted = fit == "exact_map" ? c::fitBilinearExactControlMap(model, states, count, ridge, conjugate)
            : fit == "controlled_residual" ? c::fitBilinearControlResidual(model, data, ridge, conjugate)
            : throw std::invalid_argument("bilinear input_map fit must be controlled_residual or exact_map");
        relative = fitted.diagnostics.relativeResidual; model.setControlMap(mapName, fitted.representation);
    }
    return {{"relative_residual", json::Value(relative)}, {"ridge", json::Value(ridge)},
            {"sampled_states", json::Value(double(count))}};
}
}  // namespace

std::filesystem::path runInputMap(const std::vector<std::filesystem::path>& sources, const json::Value& given,
                                  const std::filesystem::path& root, const std::string& label) {
    // A pipeline passes the whole B_psi parameter as "map"; its fields sit beside "krom".
    json::Value settings = given;
    if (const auto* map = std::get_if<json::Value::Object>(&given.at("map").storage()))
        for (const auto& [key, value] : *map) settings.object()[key] = value;
    const auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    if (kromPath.empty()) throw std::invalid_argument("input_map requires a KROM");
    auto model = c::loadKoopmanModel(kromPath);
    if (model.inputDimension() != 1) throw std::invalid_argument("input_map supports one scalar input");
    const bool history = model.coordinates().computableFromMeasurements();
    const auto representation = settings.at("representation").string(history ? "history_exact" : "exact");
    const bool analytic = representation == "exact" || representation == "history_exact";
    if (!analytic && representation != "constant" && representation != "bilinear")
        throw std::invalid_argument("input_map representation must be exact, history_exact, constant or bilinear");
    if (representation == "history_exact" && !history)
        throw std::invalid_argument("history_exact needs a causal-history KROM; use exact");
    if (representation == "exact" && history)
        throw std::invalid_argument("a causal-history KROM has no state-free exact map; use history_exact");
    const auto fit = settings.at("fit").string(analytic ? "analytic" : history ? "held_input" : "controlled_residual");
    const auto mapName = settings.at("map_name").string(representation);
    json::Value::Object diagnostics;
    std::string transformPath;
    if (!analytic) {
        if (sources.empty()) throw std::invalid_argument("a regressed input_map needs forced training records");
        // History regression historically used every pair; strict_hold opts in to the same check as other KROMs.
        const bool strictHold = settings.at("strict_hold").boolean(!history);
        if (history) {
            if (fit != "held_input") throw std::invalid_argument("causal-history input_map regression uses fit held_input");
            auto map = loadMeasurementMap(model, kromPath, settings);
            transformPath = map.path;
            const auto transform = std::move(map.transform);
            diagnostics = heldInputRegression(model, heldInputRecords(model, &transform, sources, settings, strictHold),
                                              settings, representation, mapName);
        } else if (fit == "held_input") {
            diagnostics = heldInputRegression(model, heldInputRecords(model, nullptr, sources, settings, strictHold),
                                              settings, representation, mapName);
        } else if (model.requiresExternalCoordinates()) {
            diagnostics = coordinateRegression(model, sources, settings, representation, fit, mapName);
        } else diagnostics = stateRegression(model, sources, settings, representation, fit, mapName);
    }
    json::Value::Array sourceList; for (const auto& path : sources) sourceList.emplace_back(path.generic_u8string());
    json::Value::Object inputMap{{"representation", json::Value(representation)}, {"fit", json::Value(fit)},
        {"map_name", json::Value(mapName)}, {"control_representation", json::Value(analytic ? representation : mapName)},
        {"forced_sources", json::Value(analytic ? json::Value::Array{} : std::move(sourceList))},
        {"diagnostics", json::Value(std::move(diagnostics))}};
    json::Value::Object meta{
        {"format", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("input-map")},
        {"label", json::Value(label)}, {"source_krom", json::Value(kromPath.generic_u8string())},
        {"source_archive", json::Value(kromPath.generic_u8string())}, {"input_map", json::Value(std::move(inputMap))},
        {"control_map_time_domain", json::Value("continuous-time")}, {"independent_archive", json::Value(true)}};
    if (!transformPath.empty()) meta["history_transform"] = json::Value(transformPath);
    const auto directory = objectFolder(root, "koopman-input-map");
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__input-map-" + representation
        + "__" + timestamp()) + ".npz");
    std::map<std::string, c::NpyArray> extra;  // keep the KROM's training samples for Viewer
    const auto parent = c::readNpz(kromPath);
    for (const char* name : {"training_coordinates", "training_time", "training_states", "source_mode_indices"})
        if (parent.contains(name)) extra[name] = parent.at(name);
    c::saveKoopmanModel(path, model, json::dump(json::Value(std::move(meta))), std::move(extra));
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
