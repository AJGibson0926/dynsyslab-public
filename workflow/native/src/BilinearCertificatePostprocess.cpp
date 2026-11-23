// bilinear_certificate: a frozen-A discrete bilinear model, residual sector and checked
// state-feedback gain derived from a K-ROM's stored bilinear map and held-input records.
// The K-ROM is read only. Method: docs/architecture/CERTIFIED_BILINEAR_CONTROL_TOOLS.md.
// References: Strasser, Worthmann et al., SafEDMD, arXiv:2402.03145; Goswami & Paley,
// IEEE TAC 67 (2022) 2715-2728; Bruder, Fu & Vasudevan, IEEE RA-L 6 (2021) 4369-4376;
// Van Loan, IEEE TAC 23 (1978) 395-404; Boyd et al., LMIs in System and Control Theory,
// SIAM 1994; Bernussou, Peres & Geromel, Syst. Control Lett. 13 (1989) 65-72.
#include "BilinearCertificateBinding.hpp"
#include "HeldInputRecords.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/CertifiedDesign.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanBilinearBridge.hpp>
#include <dynsyslab/core/ResidualSector.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/SupportedResidualSector.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number());
    return out;
}

json::Value array(const std::vector<double>& values) {
    json::Value::Array out;
    for (double v : values) out.emplace_back(v);
    return json::Value(std::move(out));
}

// target: {"lifted": [...]} (split psi), {"coordinates": [...]} or {"state": [...]}.
std::vector<double> liftedTarget(const c::KoopmanModel& model, const json::Value& target) {
    const auto lifted = numbers(target.at("lifted")), coordinates = numbers(target.at("coordinates")),
               state = numbers(target.at("state"));
    if (int(!lifted.empty()) + int(!coordinates.empty()) + int(!state.empty()) != 1)
        throw std::invalid_argument("bilinear_certificate target needs exactly one of lifted, coordinates or state");
    if (!lifted.empty()) {
        if (lifted.size() != model.splitDimension()) throw std::invalid_argument("target lifted dimension differs from the K-ROM");
        return lifted;
    }
    if (!coordinates.empty()) {
        if (coordinates.size() != model.coordinateDimension()) throw std::invalid_argument("target coordinates differ from the K-ROM");
        return model.evaluateCoordinatesSplit(coordinates.data());
    }
    if (model.requiresExternalCoordinates())
        throw std::invalid_argument("this K-ROM has no state lift: give the target as coordinates or lifted values");
    if (state.size() != model.stateDimension()) throw std::invalid_argument("target state dimension differs from the K-ROM");
    return model.evaluateSplit(state.data());
}

struct RegionData {
    c::BilinearDataset data;
    std::size_t total{}, outside{};
    double largestNorm{};  // of z over the kept pairs
};

// Centred pairs inside the trusted region: ||z|| <= R and |v_i| <= ubar_i.
RegionData centredPairs(const c::KoopmanModel& model, const HeldInputRecords& records,
                        const std::vector<double>& target, double targetInput, double radius, double bound) {
    const auto n = model.splitDimension();
    RegionData out;
    out.total = records.rows.size();
    for (const auto& row : records.rows) {
        auto z = model.splitValues(row.psi), next = model.splitValues(row.next);
        double squared = 0.0;
        for (std::size_t i = 0; i < n; ++i) { z[i] -= target[i]; next[i] -= target[i]; squared += z[i] * z[i]; }
        const double norm = std::sqrt(squared), v = row.u - targetInput;
        if (norm > radius || std::abs(v) > bound) { ++out.outside; continue; }
        out.largestNorm = std::max(out.largestNorm, norm);
        out.data.z.insert(out.data.z.end(), z.begin(), z.end());
        out.data.zNext.insert(out.data.zNext.end(), next.begin(), next.end());
        out.data.u.push_back(v);
        ++out.data.samples;
    }
    return out;
}
}  // namespace

std::filesystem::path runBilinearCertificate(const std::vector<std::filesystem::path>& sources,
    const json::Value& settings, const std::filesystem::path& root, const std::string& label) {
    const auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    if (kromPath.empty()) throw std::invalid_argument("bilinear_certificate requires a K-ROM");
    if (sources.empty()) throw std::invalid_argument("bilinear_certificate needs forced held-input records");
    const auto model = c::loadKoopmanModel(kromPath);
    if (model.inputDimension() != 1) throw std::invalid_argument("bilinear_certificate supports one scalar input");
    const auto representation = settings.at("control_representation").string("bilinear");
    const auto operators = settings.at("operators").string("converted");
    if (operators != "converted" && operators != "fit")
        throw std::invalid_argument("bilinear_certificate operators must be converted or fit");
    const auto bounds = numbers(settings.at("input_bounds"));
    if (bounds.size() != 1 || !(bounds[0] > 0.0) || !std::isfinite(bounds[0]))
        throw std::invalid_argument("bilinear_certificate needs one positive input bound");
    const double radius = settings.at("state_radius").number(0.0);
    if (!(radius > 0.0) || !std::isfinite(radius)) throw std::invalid_argument("bilinear_certificate needs a positive state_radius");
    if (settings.at("target").isNull()) throw std::invalid_argument("bilinear_certificate needs a target");
    auto targetInput = numbers(settings.at("target_input"));
    if (targetInput.empty()) targetInput.assign(1, 0.0);
    if (targetInput.size() != 1 || !std::isfinite(targetInput[0])) throw std::invalid_argument("target_input must be one finite value");
    if (operators == "fit" && targetInput[0] != 0.0)
        throw std::invalid_argument("operators fit needs target_input 0, because A_eff would depend on the fitted B_i");

    std::optional<c::FrozenHistoryTransform> transform;
    if (model.coordinates().kind == c::EigenfunctionCoordinateKind::causalHistory) {
        const auto transformPath = settings.at("history_transform").string(model.coordinates().providerArtifact.empty()
            ? kromPath.u8string() : model.coordinates().providerArtifact);
        transform = c::loadFrozenHistoryTransform(std::filesystem::u8path(transformPath));
        if (transform->identity() != model.coordinates().providerIdentity)
            throw std::invalid_argument("History transform differs from K-ROM training coordinates");
    }
    const auto records = heldInputRecords(model, transform ? &*transform : nullptr, sources, settings, true);
    if (records.rows.empty()) throw std::invalid_argument("bilinear_certificate found no held-input pairs");
    const double requestedDt = settings.at("sample_time").number(records.sampleTime);
    if (std::abs(requestedDt - records.sampleTime) > 1e-9 * records.sampleTime)
        throw std::invalid_argument("bilinear_certificate sample_time differs from the forced records");

    c::KoopmanBilinearBridgeOptions bridgeOptions;
    bridgeOptions.representation = representation;
    bridgeOptions.sampleTime = records.sampleTime;
    bridgeOptions.targetLifted = liftedTarget(model, settings.at("target"));
    bridgeOptions.targetInput = targetInput;
    bridgeOptions.equilibriumTolerance = settings.at("equilibrium_tolerance").number(bridgeOptions.equilibriumTolerance);
    const auto bridge = c::bridgeKoopmanBilinear(model, bridgeOptions);
    const auto n = bridge.model.stateDimension;
    const auto region = centredPairs(model, records, bridge.targetLifted, targetInput[0], radius, bounds[0]);
    const auto defect = c::conversionDefect(bridge, bounds);

    auto used = bridge.model;
    json::Value::Object fitReport;
    std::string status;
    if (!bridge.targetIsEquilibrium) status = "target_not_equilibrium";
    else if (region.data.samples < 1 + n + n) status = "insufficient_data_in_region";
    if (status.empty() && operators == "fit") {
        c::BilinearFitOptions fitOptions;
        fitOptions.ridge = settings.at("ridge").number(0.0);
        fitOptions.allowRankDeficient = settings.at("allow_rank_deficient").boolean(false);
        try {
            const auto fit = c::refitBridgeOperators(bridge, region.data, fitOptions);
            used = fit.model;
            fitReport = {{"regressors", json::Value(double(fit.regressors))}, {"rank", json::Value(double(fit.rank))},
                {"identifiable", json::Value(fit.identifiable)},
                {"condition_number", std::isfinite(fit.conditionNumber) ? json::Value(fit.conditionNumber) : json::Value()},
                {"relative_residual", json::Value(fit.stats.relative)}};
        } catch (const std::invalid_argument& failure) {
            status = "operator_fit_failed";
            fitReport = {{"error", json::Value(failure.what())}};
        }
    }

    c::ResidualSectorBound sector;
    std::string sectorMessage;
    json::Value::Object sectorReport{{"evaluated", json::Value(false)}};
    c::CertifiedDesignResult design;
    if (status.empty()) {
        c::SupportedSectorOptions sectorOptions;
        sectorOptions.inflation = settings.at("inflation").number(1.0);
        sectorOptions.zeroInputTolerance = settings.at("zero_input_tolerance").number(sectorOptions.zeroInputTolerance);
        sectorOptions.zeroResidualTolerance = settings.at("zero_residual_tolerance").number(sectorOptions.zeroResidualTolerance);
        sectorOptions.minimumCoverage = settings.at("minimum_coverage").number(sectorOptions.minimumCoverage);
        sectorOptions.stateRadius = radius;
        sectorOptions.inputBound = bounds[0];
        if (const auto& declared = settings.at("one_sided_sector"); !declared.isNull())
            sectorOptions.oneSided = c::OneSidedSectorDeclaration{declared.at("channel").string(), declared.at("reason").string()};
        const auto estimated = c::estimateSupportedResidualSector(used, region.data, sectorOptions);
        sector = estimated.bound;
        sectorReport = {{"evaluated", json::Value(true)}, {"method", json::Value("zero_input_first")},
            {"supported", json::Value(estimated.supported)}, {"feasible", json::Value(sector.feasible)},
            {"c_x", json::Value(sector.cState)}, {"c_u", json::Value(sector.cInput)},
            {"inflation", json::Value(sector.inflation)}, {"samples", json::Value(double(sector.samples))},
            {"maximum_residual", json::Value(sector.maximumResidual)},
            {"reason", json::Value(estimated.reason)}, {"one_sided", json::Value(estimated.oneSided)},
            {"declaration_reason", json::Value(estimated.declarationReason)},
            {"zero_input_pairs", json::Value(double(estimated.zeroInputPairs))},
            {"forced_pairs", json::Value(double(estimated.forcedPairs))},
            {"largest_zero_input_state_norm", json::Value(estimated.largestZeroInputStateNorm)},
            {"largest_state_norm", json::Value(estimated.largestStateNorm)},
            {"largest_input_norm", json::Value(estimated.largestInputNorm)},
            {"guarantee", json::Value(sector.guarantee)}};
        if (!estimated.supported) status = "sector_unsupported";
        sectorMessage = estimated.reason;
    }
    if (status.empty()) {
        c::CertifiedDesignSpec spec;
        spec.inputBounds = bounds;
        spec.stateRadius = radius;
        spec.contraction = settings.at("contraction").number(spec.contraction);
        spec.solverMargin = settings.at("solver_margin").number(spec.solverMargin);
        spec.certificateMargin = settings.at("certificate_margin").number(spec.certificateMargin);
        spec.solver.tolerance = settings.at("tolerance").number(spec.solver.tolerance);
        if (!settings.at("max_iterations").isNull())
            spec.solver.maxIterations = static_cast<std::size_t>(settings.at("max_iterations").number());
        design = c::designCertifiedBilinearController(used, sector, spec);
        status = c::certifiedDesignStatusName(design.status);
    }

    std::vector<double> operatorValues, continuousOperators;
    for (const auto& matrix : used.B) operatorValues.insert(operatorValues.end(), matrix.begin(), matrix.end());
    for (const auto& matrix : bridge.continuous.operators)
        continuousOperators.insert(continuousOperators.end(), matrix.begin(), matrix.end());
    c::NpzRecord record;
    const auto store = [&](const char* name, std::vector<std::size_t> shape, const std::vector<double>& values) {
        record.arrays[name] = c::NpyArray::fromDoubles(std::move(shape), values);
    };
    store("frozen_A", {n, n}, bridge.frozenTransition);
    store("model_A", {n, n}, used.A);
    store("model_B0", {n, 1}, used.B0);
    store("model_B", {1, n, n}, operatorValues);
    store("continuous_generator", {n, n}, bridge.continuous.generator);
    store("continuous_offsets", {n, 1}, bridge.continuous.offsets);
    store("continuous_operators", {1, n, n}, continuousOperators);
    store("target_lifted", {n}, bridge.targetLifted);
    store("target_input", {1}, bridge.targetInput);
    store("constant_term", {n}, bridge.constantTerm);
    if (design.certified) {
        store("certificate_K", {1, n}, design.gain);
        store("certificate_X", {n, n}, design.X);
    }

    json::Value::Array forcedSources;
    for (const auto& path : sources) forcedSources.emplace_back(path.generic_u8string());
    json::Value::Object certificate{
        {"status", json::Value(status)}, {"certified", json::Value(design.certified)},
        {"message", json::Value(!design.message.empty() ? design.message : status == "sector_unsupported" ? sectorMessage : status == "target_not_equilibrium"
            ? "The target is not an equilibrium of the bilinear model at target_input; see constant_term." : status)},
        {"control_law", json::Value("u = u* + K (z - z*), z = split(psi)")},
        {"scope", json::Value(detail::kBilinearCertificateScope)},
        {"model_identity", json::Value(c::koopmanModelIdentity(model))},
        {"control_representation", json::Value(representation)},
        {"operators", json::Value(operators)}, {"operator_fit", json::Value(std::move(fitReport))},
        {"autonomous_A", json::Value("frozen: exp(generator dt) of the source K-ROM, never fitted")},
        {"held_target_input", json::Value(bridge.heldTargetInput)},
        {"sample_time", json::Value(records.sampleTime)}, {"input_bounds", array(bounds)},
        {"state_radius", json::Value(radius)},
        {"contraction", json::Value(settings.at("contraction").number(c::CertifiedDesignSpec{}.contraction))},
        {"target_is_equilibrium", json::Value(bridge.targetIsEquilibrium)},
        {"constant_term_norm", json::Value(bridge.constantNorm)},
        {"constant_term_tolerance", json::Value(bridge.constantTolerance)},
        {"conversion_defect", json::Value(json::Value::Object{
            {"transition", json::Value(defect.transition)}, {"input", json::Value(defect.input)},
            {"meaning", json::Value("O(u^2) held-input terms the discrete bilinear family omits, at the box vertices")}})},
        {"sector", json::Value(std::move(sectorReport))},
        {"data", json::Value(json::Value::Object{
            {"held_pairs", json::Value(double(region.total))}, {"pairs_in_region", json::Value(double(region.data.samples))},
            {"pairs_outside_region", json::Value(double(region.outside))},
            {"pairs_dropped_at_input_change", json::Value(double(records.droppedAtInputChange))},
            {"largest_z_norm_in_region", json::Value(region.largestNorm)},
            {"coverage_of_radius_checked", json::Value(false)}})}};
    if (design.certified) {
        certificate["mu"] = json::Value(design.mu);
        certificate["inscribed_radius"] = json::Value(design.inscribedRadius);
        certificate["worst_lmi_eigenvalue"] = json::Value(design.certificate.worstMinEigenvalue);
        certificate["verifier_margin"] = json::Value(design.certificate.margin);
    }
    json::Value::Object meta{
        {"format", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("koopman-bilinear-certificate")},
        {"label", json::Value(label)}, {"source_krom", json::Value(kromPath.generic_u8string())},
        {"source_archive", json::Value(kromPath.generic_u8string())},
        {"forced_sources", json::Value(std::move(forcedSources))},
        {"bilinear_certificate", json::Value(std::move(certificate))}, {"independent_archive", json::Value(true)}};
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));

    const auto directory = root / L"derived" / L"koopman-bilinear-certificate";
    std::filesystem::create_directories(directory);
    const auto stem = safeName(kromPath.stem().u8string() + "__bilinear-certificate__" + timestamp());
    for (std::size_t suffix = 0; suffix < 10000; ++suffix) {
        const auto path = directory / c::wide(stem + (suffix ? "__" + std::to_string(suffix) : "") + ".npz");
        if (std::filesystem::exists(path)) continue;
        c::writeNpzAtomicNew(path, record);
        return path;
    }
    throw std::runtime_error("unable to allocate a bilinear certificate archive path");
}

}  // namespace dynsyslab::workflow::postprocess_detail
