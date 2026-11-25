#include <dynsyslab/core/ExpansionDictionary.hpp>
#include "KoopmanPolynomialContract.hpp"
#include "KoopmanSpectrumArchive.hpp"
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <climits>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
constexpr const char *key = "polynomial_contract_json";
std::size_t size(const json::Value &v) {
    const double n = v.number(-1);
    if (!std::isfinite(n) || n < 0 || n != std::floor(n) || n > double(std::size_t(-1) / 2))
        throw std::invalid_argument("Polynomial contract requires nonnegative integer sizes");
    return static_cast<std::size_t>(n);
}
std::vector<int> branches(const NpzRecord &a) {
    std::vector<int> out;
    if (a.contains("spectrum_logarithm_branches"))
        for (double v : a.at("spectrum_logarithm_branches").doubles()) {
            if (!std::isfinite(v) || v != std::floor(v) || v < INT_MIN || v > INT_MAX)
                throw std::invalid_argument("Invalid logarithm branch");
            out.push_back(static_cast<int>(v));
        }
    return out;
}
} // namespace
bool hasPolynomialContract(const NpzRecord &a) {
    return a.contains(key);
}
KoopmanModel readPolynomialContract(const NpzRecord &a, const SystemConfiguration &system) {
    const auto meta = json::parse(a.at(key).string());
    const auto schema = meta.at("schema").string();
    if (schema != "dynsyslab.polynomial-eigenfunctions/v1" && schema != "dynsyslab.eigenfunction-expansion/v1")
        throw std::invalid_argument("Unsupported eigenfunction expansion contract");
    auto q = functionDomainFromJson(meta.at("coordinates"));
    const auto &c = a.at("coefficients");
    std::shared_ptr<const ExpansionDictionary> basis;
    if (a.contains("expansion_dictionary_json")) {  // any dictionary
        basis = std::make_shared<const ExpansionDictionary>(
            ExpansionDictionary::fromJson(json::parse(a.at("expansion_dictionary_json").string())));
    } else {  // monomial arrays
        const auto &e = a.at("polynomial_exponents");
        if (e.shape.size() != 2 || e.shape[1] != q.dimension)
            throw std::invalid_argument("Polynomial contract arrays have incompatible shapes");
        PolynomialDictionary d;
        d.coordinateDimension = q.dimension;
        for (double v : e.doubles())
            d.exponents.push_back(size(json::Value(v)));
        d.center = a.at("state_center").doubles();
        d.scale = a.at("state_scale").doubles();
        basis = std::make_shared<const ExpansionDictionary>(ExpansionDictionary::monomials(std::move(d)));
    }
    if (c.shape.size() != 2 || c.shape[1] != basis->size() || basis->dimension() != q.dimension)
        throw std::invalid_argument("Expansion contract arrays have incompatible shapes");
    const auto domain = meta.at("spectrum_domain").string();
    const double dt = meta.at("sample_interval").number();
    if (domain != "discrete-time" && domain != "continuous-time")
        throw std::invalid_argument("Unknown spectrum domain");
    const auto primary = a.at("spectrum_values").complexDoubles();
    auto spectrum = domain == "continuous-time"
                        ? Spectrum::fromContinuousRates(primary, dt)
                        : Spectrum::fromStoredDiscrete(
                              primary, dt, branches(a), a.at("spectrum_rates").complexDoubles());
    if (c.shape[0] != spectrum.size())
        throw std::invalid_argument("Polynomial coefficient count differs from spectrum");
    EigenfunctionFitDiagnostics fit;
    const auto &evidence = meta.at("fit");
    fit.algorithm = evidence.at("algorithm").string();
    fit.residualDefinition = evidence.at("residual_definition").string();
    fit.sampleCount = size(evidence.at("sample_count"));
    fit.directFitCount = size(evidence.at("direct_fit_count"));
    if (a.contains("training_residual_rms"))
        fit.residualRms = a.at("training_residual_rms").doubles();
    if (a.contains("training_residual_maximum"))
        fit.residualMaximum = a.at("training_residual_maximum").doubles();
    FunctionExpansion functions(std::move(basis), c.complexDoubles(), std::move(q));
    return system.identifier.empty()
               ? KoopmanModel(std::move(functions), std::move(spectrum), size(meta.at("input_dimension")), std::move(fit))
               : KoopmanModel(std::move(functions), std::move(spectrum), system, std::move(fit));
}
void writePolynomialContract(NpzRecord &a, const KoopmanModel &model,const std::string& metadataJson) {
    // Caller-supplied evidence cannot override the model's owned optional fields.
    for(const char* keyToClear:{"dynamical_state_dimension","coordinate_observable_json",
        "coordinate_observable_names","coordinate_provider","training_residual_rms",
        "training_residual_maximum","eigenvalue_sample_time","eigenvalue_origin",
        "eigenvalue_log_branches"})a.arrays.erase(keyToClear);
    const auto &p = model.expansion();
    auto s = model.spectrum();
    // Legacy callers supplied the sampling interval in metadata rather than
    // the in-memory fit. Import that explicit interval, never infer one from physics.
    const auto metadata=metadataJson.empty()?json::Value{}:json::parse(metadataJson);
    if(s.sampleInterval()==0&&!metadata.at("delta_t").isNull()) {
        const double dt=metadata.at("delta_t").number(-1);
        if(!std::isfinite(dt)||dt<=0)throw std::invalid_argument("Invalid legacy sampling interval");
        s=Spectrum::fromContinuousRates(s.continuousRates(),dt);
    }
    const auto &fit = model.fitDiagnostics();
    const auto modes = p.outputCount();
    a.arrays["coefficients"] = NpyArray::fromComplex({modes, p.dictionary().size()}, p.coefficients());
    if (p.monomial()) {  // established polynomial arrays, read by older consumers
        const auto &d = p.monomialTable();
        std::vector<double> e;
        for (auto n : d.exponents)
            e.push_back(double(n));
        a.arrays["polynomial_exponents"] =
            NpyArray::fromDoubles({d.featureCount(), d.coordinateDimension}, e);
        a.arrays["state_center"] = NpyArray::fromDoubles({d.coordinateDimension}, d.center);
        a.arrays["state_scale"] = NpyArray::fromDoubles({d.coordinateDimension}, d.scale);
        a.arrays.erase("expansion_dictionary_json");
    } else {
        a.arrays["expansion_dictionary_json"] = NpyArray::fromString(json::dump(p.dictionary().describe()));
        for (const char *legacy : {"polynomial_exponents", "state_center", "state_scale"}) a.arrays.erase(legacy);
    }
    a.arrays["spectrum_values"] = NpyArray::fromComplex({modes}, s.values());
    std::vector<std::complex<double>> rates;
    for (std::size_t j = 0; j < modes; ++j)
        rates.push_back(s.continuousRate(j).value_or(std::complex<double>{}));
    a.arrays["spectrum_rates"] = NpyArray::fromComplex({modes}, rates);
    std::vector<double> b;
    for (auto j : s.logarithmBranches())
        b.push_back(j);
    a.arrays["spectrum_logarithm_branches"] = NpyArray::fromDoubles({b.size()}, b);
    // Retain established array keys for older consumers when mathematically valid.
    if (s.hasContinuousRates()) {
        a.arrays["continuous_eigenvalues"] = NpyArray::fromComplex({modes}, s.continuousRates());
        writePolynomialSpectrumMetadata(a, model.eigenfunctions());
    } else
        a.arrays.erase("continuous_eigenvalues");
    if(s.sampleInterval()>0) {
        a.arrays["eigenvalue_sample_time"]=NpyArray::fromDoubles({1},{s.sampleInterval()});
        a.arrays["discrete_eigenvalues"]=NpyArray::fromComplex({modes},s.discreteMultipliers());
    }
    a.arrays["eigenvalue_origin"]=NpyArray::fromString(s.domain()==EigenvalueDomain::discreteTime?"discrete-time":"continuous-time");
    const auto &q = p.coordinates();
    const auto physical = q.stateMap ? std::optional<std::size_t>(q.stateMap->stateDimension()) : q.physicalStateDimension;
    if (physical)
        a.arrays["dynamical_state_dimension"] = NpyArray::fromDoubles({1}, {double(*physical)});
    if (const auto definition = q.stateMapDefinition(); !definition.empty())
        a.arrays["coordinate_observable_json"] = NpyArray::fromString(definition);
    if (!q.names.empty()) {
        std::string names;
        for (const auto &n : q.names) {
            if (!names.empty())
                names += '\n';
            names += n;
        }
        a.arrays["coordinate_observable_names"] = NpyArray::fromString(names);
    }
    if (q.requiresExternalValues())
        a.arrays["coordinate_provider"] = NpyArray::fromString(
            q.computableFromMeasurements() ? "causal-history" : "measured");
    if (!fit.residualRms.empty())
        a.arrays["training_residual_rms"] =
            NpyArray::fromDoubles({fit.residualRms.size()}, fit.residualRms);
    if (!fit.residualMaximum.empty())
        a.arrays["training_residual_maximum"] =
            NpyArray::fromDoubles({fit.residualMaximum.size()}, fit.residualMaximum);
    a.arrays[key] = NpyArray::fromString(json::dump(json::Value::Object{
        {"schema", p.monomial() ? "dynsyslab.polynomial-eigenfunctions/v1" : "dynsyslab.eigenfunction-expansion/v1"},

        {"coordinates", functionDomainJson(q)},
        {"spectrum_domain",
         s.domain() == EigenvalueDomain::discreteTime ? "discrete-time" : "continuous-time"},
        {"sample_interval", s.sampleInterval()},
        {"input_dimension", double(model.inputDimension())},
        {"physical_system_attached", model.hasPhysicalSystem()},
        {"fit", json::Value::Object{{"algorithm", fit.algorithm},
                                    {"residual_definition", fit.residualDefinition},
                                    {"sample_count", double(fit.sampleCount)},
                                    {"direct_fit_count", double(fit.directFitCount)}}}}));
}
void writeNamedInputMaps(NpzRecord &a, const KoopmanModel &model) {
    json::Value::Array maps;
    // The legacy writer can import dt from explicit metadata. Bind to the actual
    // serialized spectrum, so that read-back does not change the map's identity.
    const auto identity=hasPolynomialContract(a)
        ?koopmanModelIdentity(readPolynomialContract(a,model.hasPhysicalSystem()?model.systemConfiguration():SystemConfiguration{}))
        :koopmanModelIdentity(model);
    for (const auto &entry : model.controlMaps())
        maps.push_back(json::Value::Object{
            {"name", entry.first},
            {"model_identity", identity},
            {"kind", std::holds_alternative<ConstantControlMap>(entry.second)
                         ? "continuous-constant"
                         : "continuous-bilinear"},
            {"provenance", std::visit([](const auto &m) { return m.provenance; }, entry.second)}});
    std::size_t index = 0;
    for (const auto &entry : model.discreteInputMaps()) {
        const auto array = "discrete_input_map_" + std::to_string(index++);
        const auto &m = entry.second;
        a.arrays[array] = NpyArray::fromComplex({m.eigenfunctionCount, m.inputDimension}, m.matrix);
        maps.push_back(json::Value::Object{{"name", entry.first},
                                           {"model_identity", identity},
                                           {"kind", "discrete"},
                                           {"array", array},
                                           {"sample_interval", m.sampleInterval},
                                           {"provenance", m.provenance}});
    }
    a.arrays["input_map_manifest_json"] = NpyArray::fromString(json::dump(maps));
}
void readNamedInputMaps(const NpzRecord &a, KoopmanModel &model) {
    if (!a.contains("input_map_manifest_json"))
        return;
    const auto manifest = json::parse(a.at("input_map_manifest_json").string());
    for (const auto &entry : manifest.array()) {
        const auto name = entry.at("name").string(), kind = entry.at("kind").string(),
                   prefix = "control_" + name;
        const auto provenance = entry.at("provenance").string();
        if (kind == "discrete")
            model.setDiscreteInputMap(name, {a.at(entry.at("array").string()).complexDoubles(),
                                             model.eigenfunctionCount(), model.inputDimension(),
                                             entry.at("sample_interval").number(), provenance,
                                             entry.at("model_identity").string()});
        else if (kind == "continuous-constant") {
            ConstantControlMap m{a.at(prefix + "_matrix").complexDoubles(),
                                 model.eigenfunctionCount(), model.inputDimension(), provenance};
            m.modelIdentity=entry.at("model_identity").string();
            if (a.contains(prefix + "_discrete_matrix"))
                m.discreteMatrix = a.at(prefix + "_discrete_matrix").complexDoubles();
            if (a.contains(prefix + "_sample_time"))
                m.sampleTime = a.at(prefix + "_sample_time").doubles().at(0);
            model.setControlMap(name, std::move(m));
        } else if (kind == "continuous-bilinear") {
            BilinearControlMap m;
            m.modes = model.eigenfunctionCount();
            m.inputs = model.inputDimension();
            m.provenance = provenance;
            m.modelIdentity=entry.at("model_identity").string();
            m.operators = a.at(prefix + "_operators").complexDoubles();
            m.offsets = a.at(prefix + "_offsets").complexDoubles();
            if (a.contains(prefix + "_conjugate_operators"))
                m.conjugateOperators = a.at(prefix + "_conjugate_operators").complexDoubles();
            model.setControlMap(name, std::move(m));
        } else
            throw std::invalid_argument("Unsupported input map contract");
    }
}
} // namespace dynsyslab::core::detail
