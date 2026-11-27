// K-ROM tool methods (pooled post-processing stages; docs/architecture/KROM_TOOLS.md):
//   control_error_compare  rank control-map candidates by their error field against the exact
//                          map; no winner is selected (ControlRegression evaluateControlMapError2D).
//   krom_rollout           autonomous K-ROM prediction psi(t) = exp(Lambda (t - t0)) psi(x0) at any
//                          increasing times, decoded by a readout (time_series object).
//   pod_galerkin           POD-Galerkin ROM from trajectory snapshots (core/PodGalerkin.hpp).
//   krom_refit             eigenfunction refit at fixed multipliers, gauge alignment and term
//                          pruning, re-validated (core/EigenfunctionRefit.hpp).
//   reconstruction_transfer  C of one K-ROM moved to another after per-mode scale alignment.
#include "KromTools.hpp"
#include "PostprocessInternal.hpp"
#include "SampledDomain.hpp"

#include <dynsyslab/core/ControlRegression.hpp>
#include <dynsyslab/core/EigenfunctionRefit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/PodGalerkin.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;
using Complex = std::complex<double>;

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

std::filesystem::path publishReport(const std::filesystem::path& root, const std::string& label, const char* method,
                                    c::NpzRecord record, json::Value::Object meta, const char* kind = nullptr) {
    meta["schema"] = "dynsyslab.derived/v1"; meta["method"] = method; meta["label"] = label;
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const char* objectKind = kind ? kind : c::object_kind::report;
    const auto directory = kindFolder(root, objectKind);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + method + "__" + timestamp()) + ".npz");
    stampDataObject(record, objectKind, method, {});
    c::writeNpzAtomicNew(path, record);
    return path;
}

std::filesystem::path publishModel(const std::filesystem::path& root, const std::string& label, const char* method,
                                   const std::filesystem::path& parent, const c::KoopmanModel& model,
                                   json::Value::Object meta) {
    meta["schema"] = "dynsyslab.derived/v1"; meta["method"] = method; meta["label"] = label;
    meta["source_krom"] = parent.generic_u8string();
    const auto directory = kindFolder(root, c::object_kind::krom);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(parent.stem().u8string() + "__" + method + "__" + timestamp()) + ".npz");
    c::saveKoopmanModel(path, model, withDataObject(json::dump(json::Value(std::move(meta))), c::object_kind::krom, {parent}));
    return path;
}

struct Run { std::vector<double> time, values; std::size_t samples{}, dimension{}; };
Run loadRun(const std::filesystem::path& path) {
    const auto archive = c::readNpz(path);
    Run run;
    run.time = archive.at("time").doubles();
    run.values = archive.at("values").doubles();
    run.samples = run.time.size();
    run.dimension = run.samples ? run.values.size() / run.samples : 0;
    return run;
}

std::filesystem::path controlErrorCompare(const json::Value& settings, const std::filesystem::path& root,
                                          const std::string& label) {
    const auto& candidates = settings.at("candidates").array();
    if (candidates.empty()) throw std::invalid_argument("control_error_compare needs candidates [{krom, map, label}]");
    const auto& grid = settings.at("grid");
    const auto a0 = numbers(grid.at("axis0")), a1 = numbers(grid.at("axis1"));
    if (a0.size() != 2 || a1.size() != 2) throw std::invalid_argument("control_error_compare grid needs axis0 and axis1 [low, high]");
    const auto columns = static_cast<std::size_t>(grid.at("columns").number(41));
    const auto rows = static_cast<std::size_t>(grid.at("rows").number(41));
    const auto fixed = numbers(settings.at("fixed_state"));
    const double floor = settings.at("relative_floor").number(0.0);
    const std::vector<std::string> metrics{"relative_mean", "relative_median", "relative_p95", "relative_maximum",
                                           "absolute_mean", "absolute_p95", "absolute_maximum"};
    std::vector<double> table(candidates.size() * metrics.size());
    json::Value::Array names;
    for (std::size_t k = 0; k < candidates.size(); ++k) {
        const auto& item = candidates[k];
        const auto model = c::loadKoopmanModel(c::wide(item.at("krom").string()));
        const auto map = item.at("map").string();
        const auto e = c::evaluateControlMapError2D(model, map, a0[0], a0[1], a1[0], a1[1], columns, rows, fixed, floor);
        const double values[] = {e.relativeMean, e.relativeMedian, e.relativeP95, e.relativeMaximum, e.absoluteMean, e.absoluteP95,
                                 e.absoluteMaximum};
        std::copy(std::begin(values), std::end(values), table.begin() + k * metrics.size());
        names.emplace_back(item.at("label").string(item.at("krom").string() + ":" + map));
    }
    const auto metric = settings.at("metric").string("relative_p95");
    const auto column = std::find(metrics.begin(), metrics.end(), metric);
    if (column == metrics.end()) throw std::invalid_argument("control_error_compare metric is unknown: " + metric);
    std::vector<std::int64_t> order(candidates.size());
    std::iota(order.begin(), order.end(), 0);
    const auto m = static_cast<std::size_t>(column - metrics.begin());
    std::stable_sort(order.begin(), order.end(), [&](std::int64_t x, std::int64_t y) {
        return table[std::size_t(x) * metrics.size() + m] < table[std::size_t(y) * metrics.size() + m]; });
    c::NpzRecord record;
    record.arrays["candidate_metrics"] = c::NpyArray::fromDoubles({candidates.size(), metrics.size()}, table);
    record.arrays["ranking"] = c::NpyArray::fromInt64s({candidates.size()}, order);
    json::Value::Array metricNames;
    for (const auto& name : metrics) metricNames.emplace_back(name);
    return publishReport(root, label, "control-error-compare", std::move(record),
        {{"candidates", json::Value(std::move(names))}, {"metric_columns", json::Value(std::move(metricNames))},
         {"ranking_metric", metric}, {"winner_selected", false}, {"settings", settings}});
}

std::filesystem::path kromRollout(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                  const std::filesystem::path& root, const std::string& label) {
    const auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    const auto model = c::loadKoopmanModel(kromPath);
    if (!model.spectrum().hasContinuousRates()) throw std::invalid_argument("krom_rollout needs continuous eigenvalues");
    auto x0 = numbers(settings.at("initial_state"));
    std::vector<double> times = numbers(settings.at("times"));
    if ((x0.empty() || times.empty()) && sources.empty())
        throw std::invalid_argument("krom_rollout needs initial_state and times, or a source trajectory");
    std::vector<double> truth;
    if (!sources.empty()) {
        const auto run = loadRun(sources.front());
        if (x0.empty()) x0.assign(run.values.begin(), run.values.begin() + run.dimension);
        if (times.empty()) times = run.time;
        truth = run.values;
    }
    for (std::size_t k = 1; k < times.size(); ++k)
        if (!(times[k] > times[k - 1])) throw std::invalid_argument("krom_rollout times must increase");
    const auto rates = model.spectrum().continuousRates();
    const auto psi0 = model.evaluate(x0.data());
    const std::size_t modes = model.liftedDimension();
    std::vector<Complex> lifted(times.size() * modes);
    for (std::size_t k = 0; k < times.size(); ++k)
        for (std::size_t j = 0; j < modes; ++j) lifted[k * modes + j] = std::exp(rates[j] * (times[k] - times[0])) * psi0[j];
    c::NpzRecord record;
    record.arrays["time"] = c::NpyArray::fromDoubles({times.size()}, times);
    record.arrays["lifted_values"] = c::NpyArray::fromComplex({times.size(), modes}, lifted);
    const auto readout = settings.at("reconstruction").string("state");
    json::Value::Object meta{{"source_krom", kromPath.generic_u8string()}, {"initial_state", settings.at("initial_state")},
                             {"rollout", "psi(t) = exp(lambda (t - t0)) psi(x0)"}};
    if (model.hasReconstruction(readout)) {
        const auto& C = model.reconstruction(readout);
        const auto decoded = c::reconstructObservable(C, lifted, times.size());
        std::vector<double> values(decoded.size());
        for (std::size_t k = 0; k < decoded.size(); ++k) values[k] = decoded[k].real();
        record.arrays["values"] = c::NpyArray::fromDoubles({times.size(), C.outputs}, values);
        meta["reconstruction"] = readout;
        if (truth.size() == values.size()) {
            double squares = 0, norm = 0;
            for (std::size_t k = 0; k < values.size(); ++k) { squares += (values[k] - truth[k]) * (values[k] - truth[k]); norm += truth[k] * truth[k]; }
            meta["relative_l2_error"] = norm > 0 ? std::sqrt(squares / norm) : 0.0;
        }
    }
    return publishReport(root, label, "krom-rollout", std::move(record), std::move(meta), c::object_kind::timeSeries);
}

std::filesystem::path podGalerkin(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                  const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("pod_galerkin needs trajectory snapshots");
    std::vector<double> snapshots; std::size_t n = 0;
    for (const auto& path : sources) {
        const auto run = loadRun(path);
        if (n && run.dimension != n) throw std::invalid_argument("pod_galerkin snapshots differ in dimension");
        n = run.dimension;
        snapshots.insert(snapshots.end(), run.values.begin(), run.values.end());
    }
    const auto source = loadSource(sources.front());
    auto system = c::createSystem(systemConfiguration(source));
    const auto model = c::fitPodGalerkin(snapshots, snapshots.size() / n, n,
        static_cast<std::size_t>(settings.at("rank").number(0)), settings.at("energy").number(0.999),
        settings.at("reference").string("mean"), *system);
    c::NpzRecord record;
    record.arrays["pod_modes"] = c::NpyArray::fromDoubles({n, model.rank}, model.modes);
    record.arrays["reference_state"] = c::NpyArray::fromDoubles({n}, model.reference);
    record.arrays["singular_values"] = c::NpyArray::fromDoubles({model.singularValues.size()}, model.singularValues);
    record.arrays["reduced_A"] = c::NpyArray::fromDoubles({model.rank, model.rank}, model.A);
    if (!model.B.empty()) record.arrays["reduced_B"] = c::NpyArray::fromDoubles({model.rank, model.inputs}, model.B);
    return publishReport(root, label, "pod-galerkin", std::move(record),
        {{"rank", double(model.rank)}, {"energy_fraction", model.energyFraction},
         {"reference", settings.at("reference").string("mean")}, {"system", systemConfiguration(source).identifier},
         {"model", "z' = Phi^T f(x_r + Phi z); A_r = Phi^T J(x_r) Phi, B_r = Phi^T B(x_r)"}},
        c::object_kind::dynamicalModel);
}

// psi at every sample of every run (q from the K-ROM's coordinate observable).
struct Lifted { std::vector<std::vector<double>> theta; std::vector<std::vector<Complex>> psi; std::vector<std::vector<double>> time; };

// Terms of one mode with |c| >= floor * max|c|, largest first: index, kind, orders and the value.
json::Value rankedTerms(const c::ExpansionDictionary& dictionary, const std::vector<Complex>& cj, double floor) {
    double largest = 0.0;
    for (const auto& v : cj) largest = std::max(largest, std::abs(v));
    std::vector<std::size_t> order;
    for (std::size_t t = 0; t < cj.size(); ++t) if (std::abs(cj[t]) >= floor * largest && largest > 0.0) order.push_back(t);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return std::abs(cj[a]) > std::abs(cj[b]); });
    json::Value::Array rows;
    for (const auto t : order) {
        const auto& term = dictionary.terms()[t];
        json::Value::Array orders;
        for (const auto o : term.orders) orders.emplace_back(double(o));
        rows.emplace_back(json::Value::Object{{"term", double(t)}, {"kind", term.kind}, {"orders", json::Value(std::move(orders))},
            {"real", cj[t].real()}, {"imag", cj[t].imag()}, {"magnitude", std::abs(cj[t])}});
    }
    return json::Value(std::move(rows));
}

std::filesystem::path kromRefit(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                const std::filesystem::path& root, const std::string& label) {
    const auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    auto model = c::loadKoopmanModel(kromPath);
    const auto& expansion = model.expansion();
    const auto& dictionary = expansion.dictionary();
    const std::size_t terms = dictionary.size(), modes = model.liftedDimension();
    const auto rates = model.spectrum().continuousRates();
    if (sources.empty()) throw std::invalid_argument("krom_refit needs sampled runs");
    const auto coordinates = model.coordinateObservable();
    std::vector<double> now, next; std::size_t pairs = 0;
    std::vector<std::vector<double>> runQ; std::vector<std::vector<double>> runTime, runStates;
    double dt = 0;
    for (const auto& path : sources) {
        const auto run = loadRun(path);
        std::vector<double> qs;
        for (std::size_t k = 0; k < run.samples; ++k) {
            const auto q = coordinates.evaluateReal(run.values.data() + k * run.dimension);
            qs.insert(qs.end(), q.begin(), q.end());
            if (k) {
                const double step = run.time[k] - run.time[k - 1];
                if (dt == 0) dt = step;
                if (std::abs(step - dt) > 1e-9 * dt) throw std::invalid_argument("krom_refit needs one uniform sample step");
            }
        }
        const std::size_t qd = model.coordinateDimension();
        for (std::size_t k = 0; k + 1 < run.samples; ++k) {
            const auto a = dictionary.evaluate(qs.data() + k * qd), b = dictionary.evaluate(qs.data() + (k + 1) * qd);
            now.insert(now.end(), a.begin(), a.end()); next.insert(next.end(), b.begin(), b.end()); ++pairs;
        }
        runQ.push_back(std::move(qs)); runTime.push_back(run.time); runStates.push_back(run.values);
    }
    auto coefficients = expansion.coefficients();  // modes x terms
    const bool refit = settings.at("refit").boolean(true);
    const double prune = settings.at("prune_threshold").number(0.0);
    const bool gauge = settings.at("gauge").boolean(true);
    // "orient": fix each eigenfunction phase on the training runs; "coefficient_summary": the terms
    // with |c| >= floor * max|c|, largest first, per mode.
    const bool orient = settings.at("orient").boolean(false);
    const double rankFloor = settings.at("coefficient_summary").number(0.0);
    json::Value::Array ranked;
    // "gauge_reference": another K-ROM (e.g. a hand-derived eigenfunction) whose scale and phase the
    // eigenfunctions take, mode by nearest eigenvalue, on the sampled runs.
    std::optional<c::KoopmanModel> reference;
    if (!settings.at("gauge_reference").string().empty())
        reference.emplace(c::loadKoopmanModel(std::filesystem::u8path(settings.at("gauge_reference").string())));
    const auto evaluateMode = [&](const std::vector<Complex>& cj) {
        std::vector<Complex> values;
        const std::size_t qd = model.coordinateDimension();
        for (const auto& qs : runQ)
            for (std::size_t k = 0; k < qs.size() / qd; ++k) {
                const auto theta = dictionary.evaluate(qs.data() + k * qd);
                Complex v = 0.0;
                for (std::size_t t = 0; t < terms; ++t) v += cj[t] * theta[t];
                values.push_back(v);
            }
        return values;
    };
    const auto anchored = [&](const std::vector<Complex>& cj, Complex rate) {
        double squares = 0, norm = 0;
        const std::size_t qd = model.coordinateDimension();
        for (std::size_t r = 0; r < runQ.size(); ++r) {
            const auto& qs = runQ[r];
            Complex first = 0.0;
            for (std::size_t k = 0; k < qs.size() / qd; ++k) {
                const auto theta = dictionary.evaluate(qs.data() + k * qd);
                Complex v = 0.0;
                for (std::size_t t = 0; t < terms; ++t) v += cj[t] * theta[t];
                if (!k) first = v;
                squares += std::norm(v - std::exp(rate * (runTime[r][k] - runTime[r][0])) * first);
                norm += std::norm(v);
            }
        }
        return norm > 0 ? std::sqrt(squares / norm) : 0.0;
    };
    std::vector<double> before(modes), after(modes), residuals(modes, 0.0), removed(modes, 0.0);
    std::vector<Complex> scales(modes, 1.0);
    for (std::size_t j = 0; j < modes; ++j) {
        std::vector<Complex> old(coefficients.begin() + j * terms, coefficients.begin() + (j + 1) * terms);
        before[j] = anchored(old, rates[j]);
        auto updated = old;
        if (refit) updated = c::refitAtMultiplier(now, next, pairs, terms, std::exp(rates[j] * dt), &residuals[j]);
        if (prune > 0.0) removed[j] = double(c::pruneCoefficients(updated, prune));
        if (reference) {  // gauge_reference: align to the reference K-ROM's mode of nearest eigenvalue
            const auto rr = reference->spectrum().continuousRates();
            std::size_t match = 0;
            for (std::size_t k = 1; k < rr.size(); ++k) if (std::abs(rr[k] - rates[j]) < std::abs(rr[match] - rates[j])) match = k;
            std::vector<Complex> target;
            for (std::size_t r = 0; r < runStates.size(); ++r) {
                const auto width = runStates[r].size() / runTime[r].size();
                for (std::size_t k = 0; k < runTime[r].size(); ++k) target.push_back(reference->evaluate(runStates[r].data() + k * width)[match]);
            }
            scales[j] = c::gaugeFactor(evaluateMode(updated), target);
            for (auto& value : updated) value *= scales[j];
        } else if (gauge || refit) {
            scales[j] = c::gaugeFactor(evaluateMode(updated), evaluateMode(old));
            for (auto& value : updated) value *= scales[j];
        }
        if (orient) {  // phase: the largest-magnitude value on the training runs becomes real positive
            const auto values = evaluateMode(updated);
            Complex pivot = 0.0;
            for (const auto& v : values) if (std::abs(v) > std::abs(pivot)) pivot = v;
            if (std::abs(pivot) > std::numeric_limits<double>::min()) {
                const auto turn = std::conj(pivot) / std::abs(pivot);
                for (auto& value : updated) value *= turn;
                scales[j] *= turn;
            }
        }
        after[j] = anchored(updated, rates[j]);
        std::copy(updated.begin(), updated.end(), coefficients.begin() + j * terms);
        if (rankFloor > 0.0) ranked.push_back(rankedTerms(dictionary, updated, rankFloor));
    }
    model.replaceRepresentation(std::make_shared<const c::FunctionExpansion>(
        expansion.dictionaryHandle(), coefficients, expansion.coordinates(), expansion.outputNames()));
    json::Value::Array b, a, s, p;
    for (std::size_t j = 0; j < modes; ++j) {
        b.emplace_back(before[j]); a.emplace_back(after[j]);
        s.emplace_back(json::Value::Array{scales[j].real(), scales[j].imag()}); p.emplace_back(removed[j]);
    }
    return publishModel(root, label, "krom-refit", kromPath, model,
        {{"refit_at_fixed_multiplier", refit}, {"gauge_aligned", gauge || refit || reference.has_value()},
         {"gauge_reference", settings.at("gauge_reference")}, {"prune_threshold", prune},
         {"anchored_error_before", json::Value(std::move(b))}, {"anchored_error_after", json::Value(std::move(a))},
         {"gauge_factors", json::Value(std::move(s))}, {"pruned_terms", json::Value(std::move(p))},
         {"sample_pairs", double(pairs)}, {"oriented", orient}, {"coefficient_summary", json::Value(std::move(ranked))}});
}

std::filesystem::path reconstructionTransfer(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                             const std::filesystem::path& root, const std::string& label) {
    const auto targetPath = std::filesystem::u8path(settings.at("krom").string());
    const auto sourcePath = std::filesystem::u8path(settings.at("source_krom").string());
    auto target = c::loadKoopmanModel(targetPath);
    const auto source = c::loadKoopmanModel(sourcePath);
    const auto name = settings.at("reconstruction").string("state");
    if (!source.hasReconstruction(name)) throw std::invalid_argument("source K-ROM has no readout " + name);
    if (sources.empty()) throw std::invalid_argument("reconstruction_transfer needs sample states (a trajectory)");
    const auto run = loadRun(sources.front());
    const auto& C = source.reconstruction(name);
    const std::size_t tm = target.liftedDimension(), sm = source.liftedDimension();
    const auto tr = target.spectrum().continuousRates(), sr = source.spectrum().continuousRates();
    std::vector<std::vector<Complex>> tv(tm), sv(sm);
    for (std::size_t k = 0; k < run.samples; ++k) {
        const auto a = target.evaluate(run.values.data() + k * run.dimension);
        const auto b = source.evaluate(run.values.data() + k * run.dimension);
        for (std::size_t j = 0; j < tm; ++j) tv[j].push_back(a[j]);
        for (std::size_t j = 0; j < sm; ++j) sv[j].push_back(b[j]);
    }
    c::ReconstructionResult moved;
    moved.outputs = C.outputs; moved.modes = tm; moved.offset = C.offset;
    moved.matrix.assign(C.outputs * tm, 0.0);
    json::Value::Array report;
    for (std::size_t j = 0; j < tm; ++j) {
        std::size_t match = 0;
        for (std::size_t k = 1; k < sm; ++k) if (std::abs(sr[k] - tr[j]) < std::abs(sr[match] - tr[j])) match = k;
        const auto alpha = c::gaugeFactor(tv[j], sv[match]);  // alpha psi_target ~= psi_source
        double residual = 0, norm = 0;
        for (std::size_t k = 0; k < run.samples; ++k) { residual += std::norm(alpha * tv[j][k] - sv[match][k]); norm += std::norm(sv[match][k]); }
        for (std::size_t o = 0; o < C.outputs; ++o) moved.matrix[o * tm + j] = C.matrix[o * sm + match] * alpha;
        report.emplace_back(json::Value::Object{{"mode", double(j)}, {"source_mode", double(match)},
            {"scale", json::Value(json::Value::Array{alpha.real(), alpha.imag()})},
            {"eigenvalue_distance", std::abs(sr[match] - tr[j])},
            {"alignment_relative_residual", norm > 0 ? std::sqrt(residual / norm) : 0.0}});
    }
    target.setReconstruction(name, moved);
    return publishModel(root, label, "reconstruction-transfer", targetPath, target,
        {{"reconstruction", name}, {"transferred_from", sourcePath.generic_u8string()},
         {"scale_alignment", json::Value(std::move(report))}});
}

}  // namespace

bool isKromToolMethod(const std::string& name) {
    return name == "control_error_compare" || name == "krom_rollout" || name == "pod_galerkin"
        || name == "krom_refit" || name == "reconstruction_transfer";
}

std::filesystem::path runKromTool(const std::string& name, const std::vector<std::filesystem::path>& sources,
                                  const json::Value& settings, const std::filesystem::path& root, const std::string& label) {
    if (name == "control_error_compare") return controlErrorCompare(settings, root, label);
    if (name == "krom_rollout") return kromRollout(sources, settings, root, label);
    if (name == "pod_galerkin") return podGalerkin(sources, settings, root, label);
    if (name == "krom_refit") return kromRefit(sources, settings, root, label);
    if (name == "reconstruction_transfer") return reconstructionTransfer(sources, settings, root, label);
    throw std::invalid_argument("unknown K-ROM tool method: " + name);
}

}  // namespace dynsyslab::workflow::postprocess_detail
