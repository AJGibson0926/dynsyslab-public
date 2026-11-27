// mode_selection, part 3: exhaustive subset search, Pareto front, archive and optional K-ROM.
#include "ModeSelectionInternal.hpp"
#include "ObjectPublishing.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace ms = mode_selection;
namespace c = dynsyslab::core;
using ms::Complex;
namespace {

// Objectives, all minimized.
constexpr const char* kMetricNames[] = {"echo_error", "predicted_state_error", "input_map_residual",
    "state_reach", "observability_balance", "pbh_margin", "hsv_echo_min", "hsv_echo_max",
    "hsv_state_min", "hsv_state_max", "state_error_r", "state_error_v", "order"};
constexpr std::size_t kMetrics = sizeof(kMetricNames) / sizeof(kMetricNames[0]);

std::vector<double> objectives(const ms::Score& s) {
    return {s.readouts.echoError, s.predictedStateError, s.inputMapResidual, -s.stateReach};
}

void subsets(std::size_t pool, std::size_t maxSize, std::vector<std::size_t>& current, std::size_t start,
             std::vector<std::vector<std::size_t>>& out) {
    if (!current.empty()) out.push_back(current);
    if (current.size() == maxSize) return;
    for (std::size_t i = start; i < pool; ++i) { current.push_back(i); subsets(pool, maxSize, current, i + 1, out); current.pop_back(); }
}

json::Value unitJson(const ms::Unit& u, std::size_t id) {
    json::Value::Array modes; for (auto m : u.modes) modes.emplace_back(double(m));
    return json::Value::Object{{"unit", double(id)}, {"modes", std::move(modes)},
        {"lambda", json::Value::Array{u.lambda.real(), u.lambda.imag()}},
        {"input_map", json::Value::Array{u.inputMap.real(), u.inputMap.imag()}},
        {"input_map_residual", u.inputMapResidual}, {"invariance_residual", u.invariance},
        {"echo_share", u.echoShare}, {"scale", u.scale}};
}

// K-ROM with the chosen units, readouts "<measurement>" and "state", and constant map "constant".
std::filesystem::path writeKrom(const std::filesystem::path& kromPath, const std::vector<ms::Unit>& units,
    const std::vector<std::size_t>& chosen, const ms::Score& score, const ms::Options& options,
    const std::filesystem::path& root, const std::string& label, json::Value::Object provenance) {
    std::vector<std::size_t> modes;
    for (auto u : chosen) for (auto m : units[u].modes) modes.push_back(m);
    const std::size_t n = modes.size(), d = score.readouts.stateStd.size();
    const auto directory = kindFolder(root, dynsyslab::core::object_kind::krom);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__mode-selection-krom__" + timestamp()) + ".npz");
    std::map<std::string, c::NpyArray> extra;
    // 1. The subset, with each pair's lower mode made exactly conj(upper). H-DMD normalizes each
    //    mode separately (lower = c conj(upper), arbitrary c); conj(upper) is equally valid for
    //    conj(lambda) and is what the conjugate readouts and B_psi below assume. These are new
    //    eigenfunctions, so the K-ROM is a whole krom, not a child (interpolations store values,
    //    expansions coefficients), rebuilt before any map is attached, because maps record the
    //    identity of the eigenfunctions they belong to.
    auto subset = c::loadKoopmanModel(kromPath, modes);
    subset.clearParent();
    auto record = c::makeKoopmanArchive(subset);
    for (const char* name : {"eigenfunction_values", "expansion_coefficients"}) {
        if (!record.contains(name)) continue;
        auto& array = record.arrays[name];
        auto values = array.complexDoubles();
        const bool modesInColumns = std::string(name) != "expansion_coefficients";
        const std::size_t rows = array.shape.at(0), cols = array.shape.size() > 1 ? array.shape[1] : 1;
        if ((modesInColumns ? cols : rows) != n || values.empty()) continue;
        std::size_t at = 0;
        for (auto u : chosen) {
            if (units[u].pair) {
                if (modesInColumns) for (std::size_t r = 0; r < rows; ++r) values[r * cols + at + 1] = std::conj(values[r * cols + at]);
                else for (std::size_t k = 0; k < cols; ++k) values[(at + 1) * cols + k] = std::conj(values[at * cols + k]);
            }
            at += units[u].pair ? 2 : 1;
        }
        array = c::NpyArray::fromComplex(array.shape, values);
    }
    // 2. Readouts and the constant map on the conjugate-consistent model. A rewritten lower mode takes
    //    its upper partner's residual: for real data a residual of (conj psi, conj lambda) equals that of
    //    (psi, lambda), under any definition.
    auto model = c::readKoopmanModel(record);
    for (std::size_t at = 0, u = 0; u < chosen.size(); at += units[chosen[u]].pair ? 2 : 1, ++u)
        if (units[chosen[u]].pair) model.setResidual(at + 1, model.spectrum().residual(at));
    c::ReconstructionResult echo{std::vector<Complex>(n), 1, n}, state{std::vector<Complex>(d * n), d, n};
    echo.offset = {score.readouts.echoOffset}; state.offset = score.readouts.stateOffset;
    echo.relativeResidual = score.readouts.echoError;
    c::ConstantControlMap map{std::vector<Complex>(n), n, 1, "mode_selection constant B_psi (forced records, per-mode controlled residual)"};
    std::size_t column = 0, slot = 0;
    for (auto u : chosen) {
        const auto& unit = units[u];
        // Real coefficients (a, b) on (Re psi, Im psi)/s equal nu psi + conj(nu) conj(psi), nu = (a - i b) / (2 s).
        auto nu = [&](const std::vector<double>& row, std::size_t offset) {
            return unit.pair ? Complex(row[offset + column], -row[offset + column + 1]) / (2.0 * unit.scale)
                             : Complex(row[offset + column] / unit.scale, 0.0);
        };
        const std::size_t k = score.readouts.echo.size();
        echo.matrix[slot] = nu(score.readouts.echo, 0);
        for (std::size_t i = 0; i < d; ++i) state.matrix[i * n + slot] = nu(score.readouts.state, i * k);
        map.matrix[slot] = unit.pair ? unit.inputMap : Complex(unit.inputMap.real(), 0.0);
        if (unit.pair) {
            echo.matrix[slot + 1] = std::conj(echo.matrix[slot]);
            for (std::size_t i = 0; i < d; ++i) state.matrix[i * n + slot + 1] = std::conj(state.matrix[i * n + slot]);
            map.matrix[slot + 1] = std::conj(map.matrix[slot]);
        }
        slot += unit.pair ? 2 : 1; column += unit.pair ? 2 : 1;
    }
    model.setReconstruction(options.measurement, echo);
    model.setReconstruction("state", state);
    model.setControlMap("constant", map);
    provenance["schema"] = "dynsyslab.derived/v1"; provenance["method"] = "mode-selection-krom";
    provenance["label"] = label; provenance["source_krom"] = kromPath.generic_u8string();
    provenance["source_archive"] = kromPath.generic_u8string(); provenance["independent_archive"] = true;
    provenance["conjugate_pairs"] = "lower mode of each pair rewritten as conj(upper)";
    c::saveKoopmanModel(path, model, withDataObject(json::dump(json::Value(std::move(provenance))), dynsyslab::core::object_kind::krom, {kromPath}), std::move(extra));
    return path;
}
}  // namespace

std::filesystem::path runModeSelection(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                       const std::filesystem::path& root, const std::string& label) {
    auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    if (kromPath.empty()) throw std::invalid_argument("mode_selection requires a candidate pool K-ROM (krom)");
    // Optional smooth pool: replace scattered (H-DMD) eigenfunctions by a polynomial least-squares fit of
    // the same values, so every score, readout, B_psi and the written K-ROM use the smooth surface.
    json::Value::Object polynomialFit;
    if (!settings.at("eigenfunction_polynomial_degree").isNull())
        kromPath = fitPolynomialEigenfunctions(kromPath,
            static_cast<std::size_t>(settings.at("eigenfunction_polynomial_degree").number(0)), root, polynomialFit);
    const auto model = c::loadKoopmanModel(kromPath);
    ms::Options options;
    options.measurement = settings.at("measurement").string(options.measurement);
    options.sensorNoise = settings.at("sensor_noise").number(options.sensorNoise);
    options.ridge = settings.at("ridge").number(options.ridge);
    options.stabilityMargin = settings.at("stability_margin").number(options.stabilityMargin);
    const auto records = ms::loadRecords(model, sources, settings, options);
    if (std::none_of(records.begin(), records.end(), [](const auto& r) { return r.forced; }) ||
        std::all_of(records.begin(), records.end(), [](const auto& r) { return r.forced; }))
        throw std::invalid_argument("mode_selection needs unforced records (readouts) and forced records (B_psi)");
    auto units = ms::buildUnits(model, records, settings, options);
    for (std::size_t i = 0; i < units.size(); ++i) {
        const auto e = ms::fitReadouts(units, {i}, records, options.ridge).echoError;
        units[i].echoShare = 1.0 - e * e;
    }
    // Pool: the units that explain the most echo on their own.
    std::vector<std::size_t> order(units.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return units[a].echoShare > units[b].echoShare; });
    const auto poolSize = std::min<std::size_t>(units.size(), std::size_t(settings.at("pool_size").number(10)));
    const auto maxUnits = std::min<std::size_t>(poolSize, std::size_t(settings.at("max_units").number(4)));
    std::vector<std::size_t> pool(order.begin(), order.begin() + std::ptrdiff_t(poolSize));
    std::vector<std::vector<std::size_t>> candidates; std::vector<std::size_t> current;
    subsets(poolSize, maxUnits, current, 0, candidates);
    for (auto& subset : candidates) for (auto& i : subset) i = pool[i];
    std::vector<ms::Score> scores;
    for (const auto& subset : candidates) scores.push_back(ms::scoreSubset(units, subset, records, options));
    // Pareto front over the four objectives; rank-sum composite to pick one.
    const std::size_t S = candidates.size();
    std::vector<bool> pareto(S, false);
    std::vector<std::vector<double>> obj(S);
    for (std::size_t i = 0; i < S; ++i) if (scores[i].valid) obj[i] = objectives(scores[i]);
    for (std::size_t i = 0; i < S; ++i) {
        if (!scores[i].valid) continue;
        bool dominated = false;
        for (std::size_t j = 0; j < S && !dominated; ++j) {
            if (i == j || !scores[j].valid) continue;
            bool better = false, worse = false;
            for (std::size_t o = 0; o < obj[i].size(); ++o) { better |= obj[j][o] < obj[i][o]; worse |= obj[j][o] > obj[i][o]; }
            dominated = better && !worse;
        }
        pareto[i] = !dominated;
    }
    std::vector<double> rankSum(S, INFINITY);
    for (std::size_t o = 0; o < 4; ++o) {
        std::vector<std::size_t> idx;
        for (std::size_t i = 0; i < S; ++i) if (scores[i].valid) idx.push_back(i);
        std::stable_sort(idx.begin(), idx.end(), [&](auto a, auto b) { return obj[a][o] < obj[b][o]; });
        for (std::size_t r = 0; r < idx.size(); ++r) rankSum[idx[r]] = (o ? rankSum[idx[r]] : 0.0) + double(r);
    }
    std::size_t best = S;
    for (std::size_t i = 0; i < S; ++i) if (pareto[i] && (best == S || rankSum[i] < rankSum[best])) best = i;
    if (best == S) throw std::invalid_argument("mode_selection: no subset could be scored");

    c::NpzRecord out;
    const std::size_t U = units.size();
    std::vector<double> uLambda, uB, uRes, uInv, uShare, uScale, uModes;
    for (const auto& u : units) {
        uLambda.insert(uLambda.end(), {u.lambda.real(), u.lambda.imag()});
        uB.insert(uB.end(), {u.inputMap.real(), u.inputMap.imag()});
        uRes.push_back(u.inputMapResidual); uInv.push_back(u.invariance); uShare.push_back(u.echoShare); uScale.push_back(u.scale);
        uModes.insert(uModes.end(), {double(u.modes[0]), u.pair ? double(u.modes[1]) : -1.0});
    }
    out.arrays["unit_eigenvalues"] = c::NpyArray::fromDoubles({U, 2}, uLambda);
    out.arrays["unit_input_map"] = c::NpyArray::fromDoubles({U, 2}, uB);
    out.arrays["unit_input_map_residual"] = c::NpyArray::fromDoubles({U}, uRes);
    out.arrays["unit_invariance_residual"] = c::NpyArray::fromDoubles({U}, uInv);
    out.arrays["unit_echo_share"] = c::NpyArray::fromDoubles({U}, uShare);
    out.arrays["unit_scale"] = c::NpyArray::fromDoubles({U}, uScale);
    out.arrays["unit_modes"] = c::NpyArray::fromDoubles({U, 2}, uModes);
    std::vector<double> members(S * maxUnits, -1.0), metrics(S * kMetrics, NAN), flags(S * 3, 0.0);
    for (std::size_t i = 0; i < S; ++i) {
        for (std::size_t j = 0; j < candidates[i].size(); ++j) members[i * maxUnits + j] = double(candidates[i][j]);
        const auto& s = scores[i];
        flags[i * 3] = s.valid; flags[i * 3 + 1] = pareto[i]; flags[i * 3 + 2] = rankSum[i];
        if (!s.valid) continue;
        const double row[kMetrics] = {s.readouts.echoError, s.predictedStateError, s.inputMapResidual, s.stateReach,
            s.observabilityBalance, s.pbhMargin, s.hsvEchoMin, s.hsvEchoMax, s.hsvStateMin, s.hsvStateMax,
            s.predictedStateErrors.size() > 0 ? s.predictedStateErrors[0] : NAN,
            s.predictedStateErrors.size() > 1 ? s.predictedStateErrors[1] : NAN, 0.0};
        std::copy(row, row + kMetrics, metrics.begin() + std::ptrdiff_t(i * kMetrics));
        std::size_t order = 0; for (auto u : candidates[i]) order += units[u].pair ? 2 : 1;
        metrics[i * kMetrics + kMetrics - 1] = double(order);
    }
    out.arrays["subset_units"] = c::NpyArray::fromDoubles({S, maxUnits}, members);
    out.arrays["subset_metrics"] = c::NpyArray::fromDoubles({S, kMetrics}, metrics);
    out.arrays["subset_flags"] = c::NpyArray::fromDoubles({S, 3}, flags);  // valid, pareto, rank sum

    json::Value::Array metricNames; for (auto* name : kMetricNames) metricNames.emplace_back(name);
    json::Value::Array unitList; for (std::size_t i = 0; i < U; ++i) unitList.push_back(unitJson(units[i], i));
    json::Value::Array chosen; for (auto u : candidates[best]) chosen.emplace_back(double(u));
    json::Value::Array provenance;
    for (const auto& r : records)
        provenance.push_back(json::Value::Object{{"record", r.path}, {"forced", r.forced}, {"samples", double(r.samples)},
            {"dt", r.dt}, {"e_dot", r.derivativeEvidence}});
    json::Value::Object bestJson{{"units", chosen}, {"echo_error", scores[best].readouts.echoError},
        {"predicted_state_error", scores[best].predictedStateError}, {"input_map_residual", scores[best].inputMapResidual},
        {"state_reach", scores[best].stateReach}, {"rank_sum", rankSum[best]}};
    json::Value::Object meta{{"method", "mode-selection"},
        {"schema", "dynsyslab.mode-selection/v1"}, {"label", label}, {"source_krom", kromPath.generic_u8string()},
        {"source_archive", kromPath.generic_u8string()}, {"independent_archive", true},
        {"measurement", options.measurement}, {"sensor_noise", options.sensorNoise},
        {"metric_names", std::move(metricNames)}, {"objectives", json::Value::Array{"echo_error", "predicted_state_error",
            "input_map_residual", "state_reach (maximized)"}},
        {"units", std::move(unitList)}, {"pool_units", double(poolSize)}, {"max_units", double(maxUnits)},
        {"subsets", double(S)}, {"pareto_count", double(std::count(pareto.begin(), pareto.end(), true))},
        {"best", bestJson}, {"records", std::move(provenance)}, {"settings", settings},
        {"eigenfunction_polynomial_fit", polynomialFit.empty() ? json::Value() : json::Value(polynomialFit)},
        {"e_dot_entry_points", json::Value::Array{
            "coordinate_transform: e' of each (e, e') coordinate record (numerical derivative of the recorded echo)",
            "psi(e, e') evaluated on records for C_e, C_x, invariance and B_psi",
            "invariance and B_psi use psi(t+h) - e^{lambda h} psi(t) over invariance_horizon (no numerical psi')",
            "online: none (the echo-only KLQE measures y = e)"}}};
    if (!settings.at("write_krom").isNull() && settings.at("write_krom").boolean(true)) {
        // "krom_units": write this unit subset instead of the composite pick.
        auto subset = candidates[best]; auto score = scores[best]; json::Value::Object pick = bestJson;
        if (!settings.at("krom_units").isNull()) {
            subset.clear();
            for (const auto& v : settings.at("krom_units").array()) {
                const double u = v.number(-1);
                if (!(u >= 0) || u >= double(units.size()) || u != std::floor(u))
                    throw std::invalid_argument("mode_selection krom_units must index units");
                subset.push_back(static_cast<std::size_t>(u));
            }
            score = ms::scoreSubset(units, subset, records, options);
            if (!score.valid) throw std::invalid_argument("mode_selection krom_units subset failed: " + score.failure);
            json::Value::Array ids; for (auto u : subset) ids.emplace_back(double(u));
            pick = {{"units", std::move(ids)}, {"echo_error", score.readouts.echoError},
                    {"predicted_state_error", score.predictedStateError}, {"input_map_residual", score.inputMapResidual},
                    {"state_reach", score.stateReach}, {"chosen_by", "krom_units"}};
        }
        json::Value::Array chosenUnits;
        for (auto u : subset) chosenUnits.push_back(unitJson(units[u], u));
        json::Value::Object krom{{"mode_selection", pick}, {"units", std::move(chosenUnits)}};
        meta["selected_krom"] = writeKrom(kromPath, units, subset, score, options, root, label,
                                          std::move(krom)).generic_u8string();
    }
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const auto directory = kindFolder(root, dynsyslab::core::object_kind::report);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__mode-selection__" + timestamp()) + ".npz");
    stampDataObject(out, dynsyslab::core::object_kind::report, "mode-selection", {});
    c::writeNpzAtomic(path, out);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
