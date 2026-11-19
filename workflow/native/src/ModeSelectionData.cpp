// mode_selection, part 1: records and per-unit quantities.
//
// Where a numerical echo derivative e' enters (for later noise studies):
//  * external-coordinate pools psi(e, e'): e' is the record's coordinate column, produced by
//    coordinate_transform (measured derivative); Record::derivativeEvidence names the method;
//  * the invariance and B_psi regressions use psi(t+h) - e^{lambda h} psi(t) over a horizon h
//    (integral form), so psi itself carries e' but psi is never differentiated numerically.
// Nothing here is used online: the echo-only KLQE measures y = e alone.
#include "ModeSelectionInternal.hpp"
#include "KromDatasetSupport.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail::mode_selection {
namespace c = dynsyslab::core;

std::vector<Record> loadRecords(const c::KoopmanModel& model, const std::vector<std::filesystem::path>& sources,
                                const json::Value& settings, const Options& options) {
    // Coordinate pools (psi(e, e') polynomials) are evaluated on a coordinate record's MEASURED
    // columns (coordinate_transform archives with source_states). On a plain state record, a pool in
    // ANALYTIC coordinates (e.g. recentred (r, r') or (e, e') with the analytic e' of the true state)
    // is evaluated through that analytic map, offline only; the record's e_dot evidence says so.
    // Measured / causal-history coordinates cannot be computed from a state and need coordinate records.
    // State pools (H-DMD / state SINDy) and scattered eigenfunctions are evaluated on states.
    bool polynomial = true;
    try { (void)model.evaluateCoordinates(std::vector<double>(model.coordinateDimension(), 0.0).data()); }
    catch (const std::logic_error&) { polynomial = false; }
    const auto kind = model.coordinates().kind;
    const bool coordinateKind = polynomial && kind != c::EigenfunctionCoordinateKind::physicalState;
    const bool analyticKind = kind == c::EigenfunctionCoordinateKind::analyticObservable;
    const double maxSamples = settings.at("max_samples_per_record").number(20000);
    // Readout records: "unforced" (default), "forced", "all", or a list of source indices. The
    // readouts are a static map psi -> (e, x); fitting them where the controller operates removes
    // the local bias of a readout fitted over a wider range.
    const auto& readoutSpec = settings.at("readout_records");
    std::vector<bool> readoutIndex(sources.size(), false);
    const std::string readoutMode = readoutSpec.isNull() ? "unforced"
        : std::holds_alternative<json::Value::Array>(readoutSpec.storage()) ? "list" : readoutSpec.string();
    if (readoutMode == "list")
        for (const auto& v : readoutSpec.array()) {
            const double k = v.number(-1);
            if (!(k >= 0) || k >= double(sources.size()) || k != std::floor(k))
                throw std::invalid_argument("mode_selection readout_records must index sources");
            readoutIndex[static_cast<std::size_t>(k)] = true;
        }
    else if (readoutMode != "unforced" && readoutMode != "forced" && readoutMode != "all")
        throw std::invalid_argument("mode_selection readout_records must be unforced, forced, all or a list");
    std::vector<Record> records;
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        const auto raw = c::readNpz(path);
        Record r; r.path = path.generic_u8string();
        r.forced = datasetDynamics(source) == "forced" || !source.inputs.empty();
        r.readout = readoutMode == "list" ? readoutIndex[records.size()]
            : readoutMode == "all" || (readoutMode == "forced") == r.forced;
        if (r.forced && source.inputDimension != 1)
            throw std::invalid_argument("mode_selection: forced records need one scalar input: " + r.path);
        // "time_window" selects unforced rows (readouts, invariance); "forced_time_window" forced rows (B_psi).
        const auto window = settings.at(r.forced ? "forced_time_window" : "time_window").array();
        const double t0 = window.size() == 2 ? window[0].number() : -INFINITY;
        const double t1 = window.size() == 2 ? window[1].number() : INFINITY;
        std::vector<double> states; std::size_t d = 0, echoColumn = 0;
        std::optional<c::Observable> echo;
        const bool coordinatePool = coordinateKind && (raw.contains("source_states") || !analyticKind);
        if (coordinatePool) {
            if (!raw.contains("source_states") || source.dimension != model.coordinateDimension())
                throw std::invalid_argument("mode_selection: a coordinate pool needs coordinate records (coordinate_transform "
                    "archives with source_states) matching its coordinates: " + r.path);
            states = raw.at("source_states").doubles(); d = raw.at("source_states").shape.at(1);
            const auto meta = json::parse(source.metadataJson);
            const auto& names = meta.at("coordinate_names").array();
            std::size_t column = names.size();
            for (std::size_t i = 0; i < names.size(); ++i) if (names[i].string() == options.measurement) column = i;
            if (column == names.size())
                throw std::invalid_argument("mode_selection: measurement " + options.measurement + " is not a coordinate of " + r.path);
            echoColumn = column;
            r.derivativeEvidence = meta.at("coordinate_transform_kind").string("") + " / "
                + meta.at("derivative_evidence").string("none");
        } else {
            states = source.values; d = source.dimension;
            echo = c::createSystem(systemConfiguration(source))->observableCatalog().at(options.measurement);
            r.derivativeEvidence = coordinateKind
                ? "analytic coordinates of the recorded true state (offline; no numerical derivative)"
                : "none (state pool)";
        }
        std::size_t kept = 0;
        for (std::size_t k = 0; k < source.samples; ++k) if (source.time[k] >= t0 && source.time[k] <= t1) ++kept;
        const std::size_t stride = std::max<std::size_t>(1, static_cast<std::size_t>(
            settings.at("sample_stride").number(std::ceil(double(kept) / maxSamples))));
        const std::size_t inputRows = source.inputDimension ? source.inputs.size() / source.inputDimension : 0;
        r.stateDimension = d; r.modes = model.liftedDimension();
        std::size_t seen = 0;
        for (std::size_t k = 0; k < source.samples; ++k) {
            if (source.time[k] < t0 || source.time[k] > t1) continue;
            if (seen++ % stride) continue;
            const double* state = states.data() + k * d;
            const auto psi = coordinatePool ? model.evaluateCoordinates(source.values.data() + k * source.dimension)
                                            : model.evaluate(state);
            r.time.push_back(source.time[k]);
            r.echo.push_back(coordinatePool ? source.values[k * source.dimension + echoColumn] : echo->evaluateReal(state).front());
            r.states.insert(r.states.end(), state, state + d);
            r.input.push_back(k < inputRows ? source.inputs[k] : 0.0);
            r.psi.insert(r.psi.end(), psi.begin(), psi.end());
        }
        r.samples = r.time.size();
        if (r.samples < 16) throw std::invalid_argument("mode_selection: too few samples in " + r.path);
        r.dt = r.time[1] - r.time[0];
        records.push_back(std::move(r));
    }
    return records;
}

std::vector<Unit> buildUnits(const c::KoopmanModel& model, const std::vector<Record>& records,
                             const json::Value& settings, const Options& options) {
    const auto rates = model.continuousEigenvalues();
    const std::size_t m = rates.size();
    const double tolerance = settings.at("imaginary_tolerance").number(1e-9);
    // Derivative-free (integral) form over a horizon h: no numerical psi' anywhere.
    //   unforced: psi(t+h) - e^{lambda h} psi(t)                       (invariance)
    //   forced:   psi(t+h) - e^{lambda h} psi(t) = B int_0^h e^{lambda (h-s)} u(t+s) ds
    const double horizon = settings.at("invariance_horizon").number(0.25);
    if (!(horizon > 0)) throw std::invalid_argument("mode_selection invariance_horizon must be positive");
    std::vector<Unit> units;
    std::vector<bool> used(m, false);
    for (std::size_t j = 0; j < m; ++j) {
        if (used[j] || rates[j].real() > -options.stabilityMargin) continue;  // growing or marginal
        Unit u; u.lambda = rates[j];
        if (std::abs(rates[j].imag()) <= tolerance) { u.modes = {j}; u.lambda = rates[j].real(); }
        else {
            std::size_t best = j; double gap = INFINITY;
            for (std::size_t k = 0; k < m; ++k)
                if (k != j && !used[k] && std::abs(rates[k] - std::conj(rates[j])) < gap) { gap = std::abs(rates[k] - std::conj(rates[j])); best = k; }
            if (best == j || gap > 1e-7 * std::max(1.0, std::abs(rates[j]))) continue;  // no stored conjugate
            u.pair = true;
            u.modes = rates[j].imag() > 0 ? std::vector<std::size_t>{j, best} : std::vector<std::size_t>{best, j};
            u.lambda = rates[u.modes[0]];
            used[best] = true;
        }
        used[j] = true;
        const std::size_t mode = u.modes[0];
        double psiEnergy = 0, samples = 0, miss = 0, change = 0, windows = 0, windowTime = 0;
        Complex numerator = 0; double gain = 0, forcedEnergy = 0;
        std::vector<std::pair<Complex, Complex>> forcedRows;
        for (const auto& r : records) {
            auto psi = [&](std::size_t i) { return r.psi[i * r.modes + mode]; };
            const std::size_t H = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(horizon / r.dt)));
            const double h = double(H) * r.dt;
            const Complex flow = std::exp(u.lambda * h);
            const std::size_t step = std::max<std::size_t>(1, H / 4);
            if (!r.forced) for (std::size_t i = 0; i < r.samples; ++i) { psiEnergy += std::norm(psi(i)); samples += 1; }
            for (std::size_t i = 0; i + H < r.samples; i += step) {
                const Complex delta = psi(i + H) - flow * psi(i);
                if (r.forced) {
                    Complex g = 0;
                    for (std::size_t k = 0; k < H; ++k)
                        g += std::exp(u.lambda * (h - (double(k) + 0.5) * r.dt)) * r.input[i + k] * r.dt;
                    numerator += delta * std::conj(g); gain += std::norm(g);
                    forcedEnergy += std::norm(delta); forcedRows.emplace_back(delta, g);
                } else {
                    miss += std::norm(delta); change += std::norm(psi(i + H) - psi(i));
                    windows += 1; windowTime = h;
                }
            }
        }
        if (!(samples > 0) || !(windows > 0)) throw std::invalid_argument("mode_selection needs unforced records longer than the horizon");
        u.scale = std::sqrt(psiEnergy / samples);
        if (!(u.scale > 0)) continue;
        u.invariance = change > 0 ? std::sqrt(miss / change) : 1.0;
        // White-noise intensity that accumulates the horizon miss: var(miss) / h, in scaled units.
        u.processIntensity = miss / windows / windowTime / (u.scale * u.scale);
        if (gain > 0) {
            u.inputMap = numerator / gain;
            double left = 0; for (const auto& [delta, g] : forcedRows) left += std::norm(delta - u.inputMap * g);
            u.inputMapResidual = forcedEnergy > 0 ? std::sqrt(left / forcedEnergy) : 1.0;
        }
        units.push_back(std::move(u));
    }
    if (units.empty()) throw std::invalid_argument("mode_selection: the pool has no stable units");
    return units;
}

}  // namespace dynsyslab::workflow::postprocess_detail::mode_selection
