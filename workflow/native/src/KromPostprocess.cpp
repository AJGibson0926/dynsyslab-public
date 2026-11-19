#include "PostprocessInternal.hpp"
#include "KromEvidence.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/KromValidationPlan.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

// Rule-based selection over the parent's stored (energy-sorted) mode order:
//   {"rule":"oscillatory_pairs","count":N}  first N nonreal modes + conjugates
//   {"rule":"leading","count":N}            first N modes, conjugates completed
// Optional "max_real_part": R skips modes whose continuous rate has real part above R
// (growing modes); absent or null keeps every mode.
// so a job can say "the dominant oscillatory pair" without inspecting indices.
std::vector<std::size_t> ruleIndices(const json::Value& rule,
                                     const dynsyslab::core::NpzRecord& parent) {
    const auto name = rule.at("rule").string();
    const double countValue = rule.at("count").number(1);
    if (!(countValue >= 1) || std::floor(countValue) != countValue)
        throw std::invalid_argument("mode_selection.count must be a positive integer");
    const auto wanted = static_cast<std::size_t>(countValue);
    const char* key = parent.contains("continuous_eigenvalues") ? "continuous_eigenvalues" : "eigenvalues";
    if (!parent.contains(key)) throw std::invalid_argument("mode_selection needs a stored spectrum");
    const auto spectrum = parent.at(key).complexDoubles();
    const double tolerance = rule.at("imaginary_tolerance").number(1e-9);
    auto partner = [&](std::size_t mode) -> std::optional<std::size_t> {
        for (std::size_t other = 0; other < spectrum.size(); ++other)
            if (other != mode && std::abs(spectrum[other] - std::conj(spectrum[mode]))
                                     <= 1e-9 + 1e-7 * std::abs(spectrum[mode])) return other;
        return std::nullopt;
    };
    const auto& maxRealSetting = rule.at("max_real_part");
    if (!maxRealSetting.isNull() && std::string(key) != "continuous_eigenvalues")
        throw std::invalid_argument("mode_selection.max_real_part needs continuous eigenvalues");
    const double maxReal = maxRealSetting.isNull() ? INFINITY : maxRealSetting.number();
    std::vector<std::size_t> result; std::size_t taken = 0;
    auto chosen = [&](std::size_t mode) {
        return std::find(result.begin(), result.end(), mode) != result.end(); };
    for (std::size_t mode = 0; mode < spectrum.size() && taken < wanted; ++mode) {
        if (chosen(mode)) continue;
        const bool oscillatory = std::abs(spectrum[mode].imag()) > tolerance;
        if (name == "oscillatory_pairs" && !oscillatory) continue;
        if (spectrum[mode].real() > maxReal) continue;
        if (name != "oscillatory_pairs" && name != "leading")
            throw std::invalid_argument("mode_selection.rule must be oscillatory_pairs or leading");
        result.push_back(mode); ++taken;
        if (oscillatory) {
            const auto other = partner(mode);
            if (!other) throw std::invalid_argument("selected oscillatory mode has no stored conjugate");
            result.push_back(*other);
        }
    }
    if (taken < wanted) throw std::invalid_argument("parent KROM has too few modes for mode_selection");
    return result;
}

// {"rule":"observable_energy","reconstructions":[names],"thresholds":[...],"caps":[...],
//  "max_real_part":0,"exclude_constant":false,"constant_tolerance":1e-6}
// For each named reconstruction C_o, rank modes by |C_o,j| * rms(psi_j) over the parent's
// training eigenfunction values and keep modes until `threshold` of that energy (capped).
// Modes whose continuous rate has real part above max_real_part (growing modes) are
// never kept; exclude_constant drops |lambda| < constant_tolerance (pair it with
// refit_reconstructions.intercept). Conjugates are completed. Mirrors the MATLAB
// suite's extract(method="observable", threshold, caps).
std::vector<std::size_t> observableEnergyIndices(const json::Value& rule,
                                                 const dynsyslab::core::NpzRecord& parent) {
    if (!parent.contains("continuous_eigenvalues") || !parent.contains("eigenfunction_values"))
        throw std::invalid_argument("observable_energy needs continuous eigenvalues and training eigenfunction values");
    const auto rates = parent.at("continuous_eigenvalues").complexDoubles();
    const auto& valuesArray = parent.at("eigenfunction_values");
    const auto values = valuesArray.complexDoubles();
    const std::size_t samples = valuesArray.shape.at(0), modes = valuesArray.shape.at(1);
    if (rates.size() != modes) throw std::invalid_argument("observable_energy: spectrum and eigenfunction values differ");
    std::vector<double> rms(modes, 0.0);
    for (std::size_t i = 0; i < samples; ++i)
        for (std::size_t j = 0; j < modes; ++j) rms[j] += std::norm(values[i * modes + j]);
    for (auto& v : rms) v = std::sqrt(v / double(samples));
    const double maxReal = rule.at("max_real_part").number(0.0);
    const bool excludeConstant = rule.at("exclude_constant").boolean(false);
    const double constantTolerance = rule.at("constant_tolerance").number(1e-6);
    const double imagTolerance = rule.at("imaginary_tolerance").number(1e-9);
    auto allowed = [&](std::size_t j) {
        return rates[j].real() <= maxReal && !(excludeConstant && std::abs(rates[j]) < constantTolerance);
    };
    const auto& names = rule.at("reconstructions").array();
    if (names.empty()) throw std::invalid_argument("observable_energy needs reconstructions");
    const auto& thresholds = rule.at("thresholds").array();
    const auto& caps = rule.at("caps").array();
    std::vector<std::size_t> result;
    auto has = [&](std::size_t j) { return std::find(result.begin(), result.end(), j) != result.end(); };
    for (std::size_t o = 0; o < names.size(); ++o) {
        const auto key = "reconstruction__" + safeArrayKey(names[o].string());
        if (!parent.contains(key)) throw std::invalid_argument("observable_energy: parent has no " + key);
        const auto& C = parent.at(key);
        if (C.shape.size() != 2 || C.shape[1] != modes) throw std::invalid_argument("observable_energy: bad " + key);
        const auto c = C.complexDoubles();
        std::vector<double> weight(modes, 0.0);
        for (std::size_t j = 0; j < modes; ++j) {
            if (!allowed(j)) continue;
            double w = 0; for (std::size_t r = 0; r < C.shape[0]; ++r) w += std::norm(c[r * modes + j]);
            weight[j] = w * rms[j] * rms[j];
        }
        std::vector<std::size_t> order(modes);
        for (std::size_t j = 0; j < modes; ++j) order[j] = j;
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return weight[a] > weight[b]; });
        double total = 0; for (double w : weight) total += w;
        if (!(total > 0)) continue;
        const double threshold = o < thresholds.size() ? thresholds[o].number() : rule.at("threshold").number(0.99);
        const double capValue = o < caps.size() ? caps[o].number() : rule.at("cap").number(double(modes));
        double cumulative = 0; std::size_t taken = 0;
        for (const auto j : order) {
            if (!(weight[j] > 0) || double(taken) >= capValue) break;
            if (!has(j)) result.push_back(j);
            ++taken; cumulative += weight[j];
            if (cumulative >= threshold * total) break;
        }
    }
    // complete conjugate pairs
    const auto chosen = result;
    for (const auto j : chosen) {
        if (std::abs(rates[j].imag()) <= imagTolerance) continue;
        std::size_t best = j; double gap = INFINITY;
        for (std::size_t k = 0; k < modes; ++k)
            if (k != j && std::abs(rates[k] - std::conj(rates[j])) < gap) { gap = std::abs(rates[k] - std::conj(rates[j])); best = k; }
        if (best == j || gap > 1e-9 + 1e-7 * std::abs(rates[j]))
            throw std::invalid_argument("observable_energy: selected oscillatory mode has no stored conjugate");
        if (!has(best)) result.push_back(best);
    }
    if (result.empty()) throw std::invalid_argument("observable_energy selected no modes");
    std::sort(result.begin(), result.end());
    return result;
}

// Refit named reconstructions on the selected modes' training eigenfunction values, optionally
// with a real intercept: h(x) ~= offset + C psi(x). Uses the parent's aligned training states.
json::Value::Array refitReconstructions(dynsyslab::core::KoopmanModel& reduced, const json::Value& refit,
                                        const dynsyslab::core::NpzRecord& parent,
                                        const std::vector<std::size_t>& indices) {
    namespace c = dynsyslab::core;
    using Complex = std::complex<double>;
    json::Value::Array report;
    if (refit.isNull()) return report;
    if (!parent.contains("eigenfunction_values") || !parent.contains("training_states"))
        throw std::invalid_argument("refit_reconstructions needs training eigenfunction values and states");
    const auto& valuesArray = parent.at("eigenfunction_values");
    const auto all = valuesArray.complexDoubles();
    const std::size_t samples = valuesArray.shape.at(0), parentModes = valuesArray.shape.at(1), k = indices.size();
    const auto states = parent.at("training_states").doubles();
    const std::size_t stateWidth = parent.at("training_states").shape.at(1);
    if (parent.at("training_states").shape.at(0) != samples)
        throw std::invalid_argument("refit_reconstructions: training states and eigenfunction values are not aligned");
    const bool intercept = refit.at("intercept").boolean(true);
    const double ridge = refit.at("ridge").number(1e-10);
    // Centered design (intercept) or raw design.
    std::vector<Complex> mean(k, 0.0), Z(samples * k);
    for (std::size_t i = 0; i < samples; ++i)
        for (std::size_t j = 0; j < k; ++j) Z[i * k + j] = all[i * parentModes + indices[j]];
    if (intercept) {
        for (std::size_t i = 0; i < samples; ++i) for (std::size_t j = 0; j < k; ++j) mean[j] += Z[i * k + j];
        for (auto& m : mean) m /= double(samples);
        for (std::size_t i = 0; i < samples; ++i) for (std::size_t j = 0; j < k; ++j) Z[i * k + j] -= mean[j];
    }
    // Normal equations G = Z^H Z (+ ridge * max diag).
    std::vector<Complex> G(k * k, 0.0);
    for (std::size_t i = 0; i < samples; ++i)
        for (std::size_t a = 0; a < k; ++a) {
            const auto za = std::conj(Z[i * k + a]);
            for (std::size_t b = 0; b < k; ++b) G[a * k + b] += za * Z[i * k + b];
        }
    double diag = 0; for (std::size_t a = 0; a < k; ++a) diag = std::max(diag, std::abs(G[a * k + a]));
    for (std::size_t a = 0; a < k; ++a) G[a * k + a] += ridge * (diag > 0 ? diag : 1.0);
    auto solve = [&](std::vector<Complex> b) {  // Gaussian elimination on a copy of G
        auto A = G;
        for (std::size_t c0 = 0; c0 < k; ++c0) {
            std::size_t p = c0;
            for (std::size_t r = c0 + 1; r < k; ++r) if (std::abs(A[r * k + c0]) > std::abs(A[p * k + c0])) p = r;
            if (!(std::abs(A[p * k + c0]) > 1e-300)) throw std::invalid_argument("refit_reconstructions is singular; add ridge");
            if (p != c0) { for (std::size_t j = 0; j < k; ++j) std::swap(A[c0 * k + j], A[p * k + j]); std::swap(b[c0], b[p]); }
            for (std::size_t r = c0 + 1; r < k; ++r) {
                const auto f = A[r * k + c0] / A[c0 * k + c0];
                for (std::size_t j = c0; j < k; ++j) A[r * k + j] -= f * A[c0 * k + j];
                b[r] -= f * b[c0];
            }
        }
        std::vector<Complex> x(k);
        for (std::size_t r = k; r-- > 0;) { Complex v = b[r]; for (std::size_t j = r + 1; j < k; ++j) v -= A[r * k + j] * x[j]; x[r] = v / A[r * k + r]; }
        return x;
    };
    const auto system = c::createSystem(reduced.systemConfiguration());
    const auto catalog = system->observableCatalog();
    for (const auto& item : refit.at("observables").array()) {
        const auto name = item.string();
        const auto observable = catalog.at(name);
        const auto out = observable.evaluateRealBatch(states, samples);
        const std::size_t q = out.size() / samples;
        std::vector<double> ymean(q, 0.0);
        if (intercept) { for (std::size_t i = 0; i < samples; ++i) for (std::size_t r = 0; r < q; ++r) ymean[r] += out[i * q + r];
                         for (auto& v : ymean) v /= double(samples); }
        c::ReconstructionResult fitted; fitted.outputs = q; fitted.modes = k; fitted.matrix.assign(q * k, 0.0);
        double residual = 0, energy = 0, worst = 0;
        for (std::size_t r = 0; r < q; ++r) {
            std::vector<Complex> b(k, 0.0);
            for (std::size_t i = 0; i < samples; ++i)
                for (std::size_t a = 0; a < k; ++a) b[a] += std::conj(Z[i * k + a]) * (out[i * q + r] - ymean[r]);
            const auto x = solve(b);
            for (std::size_t a = 0; a < k; ++a) fitted.matrix[r * k + a] = x[a];
        }
        if (intercept) {
            fitted.offset.assign(q, 0.0);
            for (std::size_t r = 0; r < q; ++r) {
                Complex shift = ymean[r];
                for (std::size_t a = 0; a < k; ++a) shift -= fitted.matrix[r * k + a] * mean[a];
                fitted.offset[r] = shift.real();
            }
        }
        for (std::size_t i = 0; i < samples; ++i)
            for (std::size_t r = 0; r < q; ++r) {
                Complex y = intercept ? Complex(fitted.offset[r]) : Complex(0.0);
                for (std::size_t a = 0; a < k; ++a) y += fitted.matrix[r * k + a] * all[i * parentModes + indices[a]];
                const double e = out[i * q + r] - y.real();
                residual += e * e; worst = std::max(worst, std::abs(e)); energy += std::pow(out[i * q + r] - ymean[r], 2);
            }
        fitted.residualRms = std::sqrt(residual / double(samples * q)); fitted.residualMaximum = worst;
        fitted.relativeResidual = energy > 0 ? std::sqrt(residual / energy) : 0.0;
        // The reduced model's mode order follows `indices` (KoopmanModel::extract keeps that order).
        reduced.setReconstruction(name, fitted);
        json::Value::Array offsets; for (double v : fitted.offset) offsets.emplace_back(v);
        report.emplace_back(json::Value::Object{{"observable", json::Value(name)}, {"intercept", json::Value(intercept)},
            {"offset", json::Value(std::move(offsets))}, {"relative_residual", json::Value(fitted.relativeResidual)},
            {"residual_rms", json::Value(fitted.residualRms)}});
    }
    return report;
}

std::vector<std::size_t> modeIndices(const json::Value& settings,
                                     const dynsyslab::core::NpzRecord& parent) {
    if (!settings.at("mode_selection").isNull()) {
        if (!settings.at("selected_modes").isNull())
            throw std::invalid_argument("give selected_modes or mode_selection, not both");
        if (settings.at("mode_selection").at("rule").string() == "observable_energy")
            return observableEnergyIndices(settings.at("mode_selection"), parent);
        return ruleIndices(settings.at("mode_selection"), parent);
    }
    std::vector<std::size_t> result;
    for (const auto& value : settings.at("selected_modes").array()) {
        const double number=value.number(-1.0);
        const auto index=static_cast<std::size_t>(number);
        if (number<0.0 || number!=static_cast<double>(index))
            throw std::invalid_argument("selected_modes must contain nonnegative integers");
        result.push_back(index);
    }
    if (result.empty())
        throw std::invalid_argument("KROM extraction requires selected_modes");
    return result;
}

json::Value extractionMetadata(const std::filesystem::path& source,
    const std::string& label,const std::vector<std::size_t>& indices,
    const dynsyslab::core::KoopmanModel& model,const json::Value& parent) {
    json::Value::Array selected;
    for (const auto index : indices) selected.emplace_back(static_cast<double>(index));
    json::Value::Array maps;
    if(model.hasControlMap("exact"))maps.emplace_back("exact");
    for (const auto& map : model.controlMaps()) maps.emplace_back(map.first);
    json::Value::Array reconstructions;
    for (const auto& item : model.reconstructions()) reconstructions.emplace_back(item.first);
    json::Value::Object backend;
    backend["name"]="native-cpp"; backend["compiled"]=true;
    backend["python_required"]=false;
    json::Value::Object metadata;
    metadata["format"]="dynsyslab.derived/v1";
    metadata["label"]=label;
    metadata["method"]="krom-extract";
    metadata["source_archive"]=source.generic_u8string();
    metadata["source_krom"]=source.generic_u8string();
    metadata["selected_modes"]=json::Value(std::move(selected));
    metadata["control_maps"]=json::Value(std::move(maps));
    metadata["reconstruction_observables"]=json::Value(std::move(reconstructions));
    metadata["independent_archive"]=true;
    const double deltaTime=parent.at("delta_t").number();
    if (deltaTime>0.0) metadata["delta_t"]=deltaTime;
    metadata["eigenvalue_domain"]=model.spectrum().domain()==dynsyslab::core::EigenvalueDomain::discreteTime
        ?"discrete-time":"continuous-time";
    metadata["execution_backend"]=json::Value(std::move(backend));
    return json::Value(std::move(metadata));
}

}  // namespace

std::filesystem::path runKromExtraction(
    const std::filesystem::path& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    if (!std::filesystem::is_regular_file(source))
        throw std::invalid_argument("KROM extraction source does not exist");
    const auto parentArchive=dynsyslab::core::readNpz(source);
    const auto indices=modeIndices(settings,parentArchive);
    const auto parentMetadata=parentArchive.metadataJson().empty()
        ? json::Value{}:json::parse(parentArchive.metadataJson());
    const auto parent=dynsyslab::core::loadKoopmanModel(source);
    auto reduced=parent.extract(indices);
    const auto refitReport=refitReconstructions(reduced,settings.at("refit_reconstructions"),parentArchive,indices);
    const auto directory=root/L"derived"/L"koopman-krom-extract";
    std::filesystem::create_directories(directory);
    const auto path=directory/dynsyslab::core::wide(safeName(
        source.stem().u8string()+"__krom-extract__"+timestamp())+".npz");
    const auto inherited=inheritedSourceModes(parentArchive,indices,reduced.liftedDimension());
    std::vector<double> sourceIndices;
    for (const auto index : inherited) sourceIndices.push_back(static_cast<double>(index));
    std::map<std::string,dynsyslab::core::NpyArray> extra;
    extra["source_mode_indices"]=dynsyslab::core::NpyArray::fromDoubles(
        {sourceIndices.size()},sourceIndices);
    inheritKromEvidence(parentArchive,indices,extra);
    auto meta=extractionMetadata(source,label,indices,reduced,parentMetadata);
    const auto origin=parentMetadata.at("extraction_algorithm").string(parentMetadata.at("method").string());
    meta.object()["extraction_algorithm"]=origin;
    if(!refitReport.empty())meta.object()["refit_reconstructions"]=json::Value(refitReport);
    if(origin=="edmdc") {
        const auto plan=planAutonomousKoopmanValidation(reduced,path,root);
        for(const auto& ref:dynsyslab::core::koopmanDatasetReferencesFromMetadata(plan))reduced.addDatasetReference(ref);
        meta.object()["autonomous_validation_plan"]=plan;
    }
    dynsyslab::core::saveKoopmanModel(path,reduced,json::dump(meta),std::move(extra));
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
