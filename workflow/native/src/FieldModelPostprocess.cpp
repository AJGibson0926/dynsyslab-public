// Post-processing method field_model: a data-driven vector field of any field-model kind
// (core/FieldModels.hpp, docs/architecture/FIELD_MODELS.md), pooled over the job's sources and
// scored on held-out records by rate error and free rollouts.
//
//   {"name": "field_model", "settings": {
//      "kind": "polynomial" | "compact_radial" | "local_tangent" | "paired_phase"
//            | "second_order" | "pair_network" | "interior",
//      "rates": "finite_difference" | "stored",   stored: the sources' coordinate_derivatives
//      "library": {"family": "polynomial", "degree": 2, "centers": 24, "harmonics": 5,
//                  "rational_widths": [...], "ridge": 1e-6, "threshold": 0.003, "time_weighted": false},
//      "options": {...},                          the kind's options (FIELD_MODELS.md)
//      "validation_sources": [archives],          held-out records, never fitted
//      "rollout": {"horizon": 2.0, "starts": 4, "maximum_step": 0},
//      "calibrate": {"duration": 1.0, ...},       polynomial only: fixed-support shooting refit
//      "controls": {"maximum_step": h, "restart_time": t, "compare_from": t,
//                   "visits": [[t1, t2], ...]}}}  replay controls on the held-out records
//
// Rates are fourth-order central differences of each record (two samples dropped at each end)
// unless stored. Publishes a dynamical_model object: field_model_json (the definition), metadata
// sections "training" and "validation" (rate relative L2, rollout relative RMS and completion) and,
// when calibrated, "calibration"; with controls, "controls" (halved-step change, restart against
// the original rollout on common times, visit separations).
#include "FieldModelPostprocess.hpp"
#include "KromDatasetSupport.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/FieldModels.hpp>
#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

c::FieldRecord fieldRecord(const std::filesystem::path& path, const std::string& rates) {
    const auto source = loadSource(path);
    if (datasetDynamics(source) == "forced") throw std::invalid_argument("field_model fits autonomous fields; " + path.u8string() + " is forced");
    const auto n = source.dimension(), N = source.samples();
    if (N < 8) throw std::invalid_argument("field_model: too few samples in " + path.u8string());
    const auto series = c::TimeSeries::view(source.time(), source.values().data(), N, n);
    const auto q = c::cropSamples(series, 2, N - 4);
    c::FieldRecord r{q.timeVector(), q.valueVector(), {}, n, n};
    if (rates == "stored") {
        const auto a = c::readNpz(path);
        if (!a.contains("coordinate_derivatives")) throw std::invalid_argument("field_model rates stored: no coordinate_derivatives in " + path.u8string());
        const auto& v = a.at("coordinate_derivatives");
        if (v.shape != std::vector<std::size_t>{N, n}) throw std::invalid_argument("field_model: stored derivatives must align with the samples");
        r.label = c::cropSamples(c::TimeSeries(source.time(), v.doubles(), n), 2, N - 4).valueVector();
    } else if (rates == "finite_difference") {
        r.label = c::differentiateTimeSeries(series).valueVector();
    } else throw std::invalid_argument("field_model rates are finite_difference or stored");
    if (r.label.size() != r.state.size()) throw std::invalid_argument("field_model: rates do not align with the states");
    return r;
}

c::MapFitOptions mapOptions(const J& s) {
    c::MapFitOptions o;
    o.family = s.at("family").string(o.family);
    o.degree = std::size_t(s.at("degree").number(double(o.degree)));
    o.centers = std::size_t(s.at("centers").number(double(o.centers)));
    o.harmonics = std::size_t(s.at("harmonics").number(double(o.harmonics)));
    if (!s.at("rational_widths").isNull()) { o.rationalWidths.clear(); for (const auto& w : s.at("rational_widths").array()) o.rationalWidths.push_back(w.number()); }
    o.ridge = s.at("ridge").number(o.ridge);
    o.threshold = s.at("threshold").number(o.threshold);
    o.timeWeighted = s.at("time_weighted").boolean(false);
    return o;
}

std::vector<double> numbers(const J& v) {
    std::vector<double> out;
    if (v.isNull()) return out;
    if (std::holds_alternative<double>(v.storage())) return {v.number()};
    for (const auto& x : v.array()) out.push_back(x.number());
    return out;
}

J fit(const std::string& kind, const std::vector<c::FieldRecord>& records, const J& library, const J& o) {
    const auto map = mapOptions(library);
    if (kind == "polynomial") return c::fitMap(records, map);
    if (kind == "compact_radial") {
        c::CompactRadialOptions co;
        co.stride = std::size_t(o.at("stride").number(double(co.stride))); co.neighbors = std::size_t(o.at("neighbors").number(double(co.neighbors)));
        co.asinhScale = o.at("asinh_scale").number(0); co.affineTail = o.at("affine_tail").boolean(false);
        co.ridge = o.at("ridge").number(co.ridge); co.threshold = o.at("threshold").number(co.threshold);
        co.timeWeighted = o.at("time_weighted").boolean(false); co.speedScale = o.at("speed_scale").number(0);
        return c::fitCompactRadialField(records, co);
    }
    if (kind == "local_tangent") {
        c::LocalTangentOptions lo;
        lo.stride = std::size_t(o.at("stride").number(double(lo.stride))); lo.samples = std::size_t(o.at("samples").number(double(lo.samples)));
        lo.neighbors = std::size_t(o.at("neighbors").number(double(lo.neighbors))); lo.decay = o.at("decay").number(lo.decay);
        lo.ridge = o.at("ridge").number(lo.ridge); lo.threshold = o.at("threshold").number(lo.threshold);
        return c::fitLocalTangentField(records, lo);
    }
    if (kind == "paired_phase") {
        c::PairedPhaseOptions po;
        for (const auto h : numbers(o.at("harmonics"))) po.harmonics.push_back(int(h));
        po.phase = o.at("phase").boolean(true);
        for (const auto w : numbers(o.at("weak_windows"))) po.weakWindows.push_back(std::size_t(w));
        po.instantWeight = o.at("instant_weight").number(1); po.ridge = o.at("ridge").number(po.ridge); po.threshold = o.at("threshold").number(po.threshold);
        return c::fitPairedPhaseField(records, po);
    }
    if (kind == "second_order") return c::fitSecondOrderField(records, map);
    if (kind == "pair_network") return c::fitPairNetwork(records, numbers(o.at("lags")), numbers(o.at("weights")), map, o.at("tolerance").number(1e-12));
    if (kind == "interior") {
        c::InteriorOptions io;
        io.gamma = numbers(o.at("gamma"));
        io.outer = o.at("outer").isNull() ? map : mapOptions(o.at("outer"));
        io.drive = o.at("drive").isNull() ? map : mapOptions(o.at("drive"));
        io.waveformDrive = o.at("waveform_drive").boolean(false);
        io.feedback = !o.at("feedback").isNull();
        if (io.feedback) io.feedbackOptions = mapOptions(o.at("feedback"));
        return c::fitInteriorField(records, io);
    }
    throw std::invalid_argument("field_model kind " + kind + " is not known");
}

// Rate relative L2 and free-rollout relative RMS from evenly spaced starts over `horizon`.
J score(const c::Observable& field, const std::vector<c::FieldRecord>& records, const std::vector<std::string>& names, const J& rollout) {
    const double horizon = rollout.at("horizon").number(0), step = rollout.at("maximum_step").number(0);
    const auto starts = std::size_t(rollout.at("starts").number(4));
    double error = 0, energy = 0, rollError = 0, rollEnergy = 0;
    std::size_t completed = 0, attempted = 0;
    J::Array rows;
    for (std::size_t r = 0; r < records.size(); ++r) {
        const auto& rec = records[r];
        const auto n = rec.dimension;
        double e = 0, w = 0;
        for (std::size_t k = 0; k < rec.samples(); ++k) {
            const auto v = field.evaluateReal(rec.state.data() + k * n);
            for (std::size_t j = 0; j < n; ++j) { e += std::pow(v[j] - rec.label[k * n + j], 2); w += std::pow(rec.label[k * n + j], 2); }
        }
        error += e; energy += w;
        J::Object row{{"source", names[r]}, {"rate_relative_l2", w > 0 ? std::sqrt(e / w) : 0.0}};
        if (horizon > 0 && starts) {
            double re = 0, rw = 0;
            std::size_t done = 0, tried = 0;
            const double last = rec.time.back() - horizon;
            for (std::size_t s = 0; s < starts && rec.time.front() <= last; ++s) {
                const double t0 = rec.time.front() + (starts == 1 ? 0.0 : (last - rec.time.front()) * double(s) / double(starts - 1));
                std::size_t i = 0;
                while (i + 1 < rec.samples() && rec.time[i] < t0) ++i;
                std::size_t stop = i;
                while (stop < rec.samples() && rec.time[stop] <= rec.time[i] + horizon) ++stop;
                if (stop < i + 2) continue;
                ++tried;
                const std::vector<double> clock(rec.time.begin() + std::ptrdiff_t(i), rec.time.begin() + std::ptrdiff_t(stop));
                const std::vector<double> x0(rec.state.begin() + std::ptrdiff_t(i * n), rec.state.begin() + std::ptrdiff_t((i + 1) * n));
                const auto path = c::rolloutField(field, x0, clock, step);
                if (path.complete) ++done;
                for (std::size_t k = 0; k < clock.size(); ++k)
                    for (std::size_t j = 0; j < n; ++j) {
                        const double truth = rec.state[(i + k) * n + j];
                        // A failed rollout counts its missing samples as the whole signal.
                        const double d = k * n + j < path.state.size() ? path.state[k * n + j] - truth : truth;
                        re += d * d; rw += truth * truth;
                    }
            }
            rollError += re; rollEnergy += rw; completed += done; attempted += tried;
            row["rollout_relative_rms"] = rw > 0 ? std::sqrt(re / rw) : 0.0;
            row["rollouts_complete"] = double(done); row["rollouts"] = double(tried);
        }
        rows.push_back(J(std::move(row)));
    }
    J::Object out{{"records", J(std::move(rows))}, {"rate_relative_l2", energy > 0 ? std::sqrt(error / energy) : 0.0}};
    if (horizon > 0) {
        out["rollout_relative_rms"] = rollEnergy > 0 ? std::sqrt(rollError / rollEnergy) : 0.0;
        out["rollout_completion"] = attempted ? double(completed) / double(attempted) : 0.0;
    }
    return J(std::move(out));
}

// Per-component relative RMS of `prediction` against rows [first, first + count) of `truth`.
J relativeRms(const std::vector<double>& truth, std::size_t first, const std::vector<double>& prediction, std::size_t predFirst,
              std::size_t count, std::size_t n) {
    J::Array out;
    for (std::size_t j = 0; j < n; ++j) {
        double e = 0, w = 0;
        for (std::size_t k = 0; k < count; ++k) {
            const double t = truth[(first + k) * n + j];
            e += std::pow(prediction[(predFirst + k) * n + j] - t, 2); w += t * t;
        }
        out.push_back(w > 0 ? J(std::sqrt(e / w)) : J());
    }
    return J(std::move(out));
}

std::size_t sampleAt(const std::vector<double>& time, double t) {
    const auto at = std::lower_bound(time.begin(), time.end(), t);
    if (at == time.end()) throw std::invalid_argument("field_model controls: a time outside the record");
    return std::size_t(at - time.begin());
}

// Replay controls of a fitted field on one record (the side-car audit.replay_controls and
// controls.run_controls): the rollout with the step halved (integration-error control), a
// restart from the record's state at restart_time compared with the original rollout on common
// times, and the separation of the record's states at pairs of visit times (cubic Hermite
// interpolation with the record's rates; Fritsch and Carlson 1980).
J controls(const c::Observable& field, const c::FieldRecord& rec, const J& s) {
    const auto n = rec.dimension, N = rec.samples();
    double step = s.at("maximum_step").number(0);
    if (!(step > 0)) step = (rec.time[1] - rec.time[0]) / 4;
    const std::vector<double> x0(rec.state.begin(), rec.state.begin() + std::ptrdiff_t(n));
    const auto base = c::rolloutField(field, x0, rec.time, step), fine = c::rolloutField(field, x0, rec.time, step / 2);
    J::Object out{{"maximum_step", step}, {"baseline_complete", base.complete}, {"halved_step_complete", fine.complete}};
    if (base.complete && fine.complete) {
        double e = 0, w = 0;
        for (std::size_t i = 0; i < rec.state.size(); ++i) { e += std::pow(base.state[i] - fine.state[i], 2); w += rec.state[i] * rec.state[i]; }
        out["step_change_over_record_rms"] = w > 0 ? std::sqrt(e / w) : 0.0;
    }
    if (!s.at("restart_time").isNull()) {
        const auto i = sampleAt(rec.time, s.at("restart_time").number());
        if (i + 1 >= N) throw std::invalid_argument("field_model controls: restart falls outside the usable record");
        const std::vector<double> clock(rec.time.begin() + std::ptrdiff_t(i), rec.time.end());
        const std::vector<double> xi(rec.state.begin() + std::ptrdiff_t(i * n), rec.state.begin() + std::ptrdiff_t((i + 1) * n));
        const auto restart = c::rolloutField(field, xi, clock, step);
        out["restart_at"] = rec.time[i];
        out["restart_complete"] = restart.complete;
        if (base.complete && restart.complete) {
            const auto j = sampleAt(rec.time, std::max(rec.time[i], s.at("compare_from").number(rec.time[i])));
            out["common_start"] = rec.time[j];
            out["original_on_common"] = relativeRms(rec.state, j, base.state, j, N - j, n);
            out["restart_on_common"] = relativeRms(rec.state, j, restart.state, j - i, N - j, n);
        }
    }
    J::Array visits;
    for (const auto& pair : s.at("visits").array()) {
        std::vector<double> at[2];
        for (int v = 0; v < 2; ++v) {
            const double t = pair.array().at(std::size_t(v)).number();
            if (t < rec.time.front() || t > rec.time.back()) throw std::invalid_argument("field_model controls: visits must lie on the record");
            auto k = std::min(sampleAt(rec.time, t), N - 1);
            if (k > 0 && rec.time[k] > t) --k;
            const auto k1 = std::min(k + 1, N - 1);
            const double h = rec.time[k1] - rec.time[k], u = h > 0 ? (t - rec.time[k]) / h : 0.0;
            const double h00 = 2 * u * u * u - 3 * u * u + 1, h10 = u * u * u - 2 * u * u + u, h01 = -2 * u * u * u + 3 * u * u, h11 = u * u * u - u * u;
            for (std::size_t j = 0; j < n; ++j)
                at[v].push_back(h00 * rec.state[k * n + j] + h10 * h * rec.label[k * n + j] + h01 * rec.state[k1 * n + j] + h11 * h * rec.label[k1 * n + j]);
        }
        J::Array difference, distance;
        double cumulative = 0;
        for (std::size_t j = 0; j < n; ++j) {
            const double d = at[0][j] - at[1][j];
            cumulative += d * d;
            difference.push_back(d); distance.push_back(std::sqrt(cumulative));
        }
        visits.push_back(J::Object{{"times", pair}, {"difference", J(std::move(difference))}, {"distance", J(std::move(distance))}});
    }
    if (!visits.empty()) out["visits"] = J(std::move(visits));
    return J(std::move(out));
}

}  // namespace

std::filesystem::path runFieldModel(const std::vector<std::filesystem::path>& sources, const J& settings,
                                    const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("field_model needs sources");
    const auto kind = settings.at("kind").string("polynomial"), rates = settings.at("rates").string("finite_difference");
    std::vector<c::FieldRecord> training, validation;
    std::vector<std::string> trainingNames, validationNames;
    for (const auto& s : sources) { training.push_back(fieldRecord(s, rates)); trainingNames.push_back(s.generic_u8string()); }
    for (const auto& v : settings.at("validation_sources").array()) {
        validation.push_back(fieldRecord(std::filesystem::u8path(v.string()), rates)); validationNames.push_back(v.string());
    }
    auto definition = fit(kind, training, settings.at("library"), settings.at("options"));
    J calibration;
    if (!settings.at("calibrate").isNull()) {
        if (kind != "polynomial") throw std::invalid_argument("field_model calibrate refits polynomial fields only");
        const auto& s = settings.at("calibrate");
        c::CalibrationOptions co;
        co.duration = s.at("duration").number(co.duration);
        co.startsPerRecord = std::size_t(s.at("starts_per_record").number(double(co.startsPerRecord)));
        co.maximumEvaluations = std::size_t(s.at("maximum_evaluations").number(double(co.maximumEvaluations)));
        co.coefficientPenalty = s.at("coefficient_penalty").number(co.coefficientPenalty);
        co.maximumStep = s.at("maximum_step").number(0); co.changeBound = s.at("change_bound").number(co.changeBound);
        co.globalWeights = numbers(s.at("global_weights")); co.localWeights = numbers(s.at("local_weights"));
        co.floors = numbers(s.at("floors")); co.waveformDuration = s.at("waveform_duration").number(0);
        const auto result = c::calibrateField(definition, training, co);
        definition = result.definition;
        calibration = J::Object{{"recipe", "training-shooting"}, {"accepted", result.accepted}, {"shots", double(result.shots)},
            {"fitted_terms", double(result.terms)}, {"initial_loss", result.initialLoss}, {"candidate_loss", result.candidateLoss},
            {"evaluations", double(result.evaluations)}};
    }
    const auto field = c::fieldFromDefinition(definition);
    const auto first = loadSource(sources.front());
    auto meta = json::parse(metadata(first, "field-model", label, "\"schema\":\"dynsyslab.field-model/v1\"")).object();
    meta["settings"] = settings;
    meta["field_kind"] = kind;
    meta["training"] = score(field, training, trainingNames, settings.at("rollout"));
    if (!validation.empty()) meta["validation"] = score(field, validation, validationNames, settings.at("rollout"));
    if (!calibration.isNull()) meta["calibration"] = calibration;
    if (!settings.at("controls").isNull()) {
        // On the held-out records when there are any, else on the training records.
        const auto& records = validation.empty() ? training : validation;
        const auto& names = validation.empty() ? trainingNames : validationNames;
        J::Object byRecord;
        for (std::size_t r = 0; r < records.size(); ++r) byRecord[names[r]] = controls(field, records[r], settings.at("controls"));
        meta["controls"] = J(std::move(byRecord));
    }
    J::Array inputs;
    for (const auto& name : trainingNames) inputs.push_back(name);
    meta["sources"] = inputs;
    meta["references"] = J::Array{"Brunton, Proctor and Kutz 2016, PNAS 113:3932-3937", "Wendland 1995, Adv. Comput. Math. 4:389-396",
        "Zhang and Zha 2004, SIAM J. Sci. Comput. 26:313-338", "Messenger and Bortz 2021, J. Comput. Phys. 443:110525"};
    c::NpzRecord archive;
    archive.arrays["field_model_json"] = c::NpyArray::fromString(json::dump(definition));
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, first, c::object_kind::dynamicalModel, "field-model-" + kind, std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
