#include "ArchiveCatalog.hpp"
#include "HierarchyCatalog.hpp"
#include "ModeWeightPayload.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cctype>
#include <map>
#include <set>

namespace dsl {
namespace {
namespace json=dynsyslab::core::json;

constexpr const char* kSourcePrefix="archive-source:";

const std::vector<dynsyslab::core::SystemDescriptor>& knownSystems() {
    // Extensions are loaded before the catalog is built. Reusing this process-
    // local snapshot avoids reconstructing and sorting the registry per archive.
    static const auto value = dynsyslab::core::availableSystems();
    return value;
}

void appendUnique(std::vector<std::string>& values,const std::string& value) {
    if (!value.empty() && std::find(values.begin(),values.end(),value)==values.end())
        values.push_back(value);
}

void prependUnique(std::vector<std::string>& values,const std::string& value) {
    values.erase(std::remove(values.begin(),values.end(),value),values.end());
    values.insert(values.begin(),value);
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::string systemLabel(const json::Value& metadata) {
    const auto explicitLabel=metadata.at("system_label").string();
    if (!explicitLabel.empty()) return explicitLabel;
    const auto identifier=metadata.at("system").string(
        metadata.at("system_configuration").at("identifier").string());
    for (const auto& descriptor : knownSystems())
        if (descriptor.identifier==identifier) return descriptor.label;
    return identifier;
}

std::string joined(const json::Value& values) {
    std::string result;
    for (const auto& value : values.array()) {
        if (!result.empty()) result+=", ";
        const auto text=value.string();
        result+=text.empty()?std::to_string(static_cast<std::size_t>(value.number())):text;
    }
    return result;
}

void relation(CatalogItem& item,const std::string& path,
              std::string kind,const std::string& label={}) {
    if (path.empty()) return;
    const auto source=std::filesystem::u8path(path);
    const std::string target=std::string(kSourcePrefix)+normalizedSource(source);
    for (const auto& existing : item.relations)
        if (existing.target==target && existing.kind==kind) return;
    item.relations.push_back({target,std::move(kind),label.empty()
        ? source.stem().u8string():label});
}

void catalogRelation(CatalogItem& item,const std::string& target,
                     std::string kind,const std::string& label) {
    if (target.empty() || target==item.id) return;
    for (const auto& existing:item.relations)
        if (existing.target==target && existing.kind==kind) return;
    item.relations.push_back({target,std::move(kind),label});
}

std::string reverseKind(const std::string& kind) {
    if (kind=="extracted from") return "has extracted KROM";
    if (kind=="uses KROM") return "used by derived model";
    if (kind=="expands training domain of") return "has domain-expansion data";
    if (kind=="controlled by KROM") return "used by control attempt";
    if (kind=="evaluates KROM") return "has validation result";
    if (kind=="initialized from spectrum") return "initializes optimized spectrum";
    if (kind=="fit from trajectory") return "used by spectral fit";
    return "produced derived result";
}

bool decomposition(const std::string& method) {
    static const std::set<std::string> methods{"svd","pca","pod","bagged-pca",
        "bagged_pca","robust-pca","robust_pca","hankel-pod","hankel_pod","hankel-svd","hankel_svd",
        "havok"};
    return methods.count(method)!=0;
}

bool kromMethod(const std::string& method) {
    static const std::set<std::string> methods{"koopman-with-sindy","sindy","koopman-control-regression",
        "control_regression","krom-extract","krom_extract","koopman-krom-extract"};
    return methods.count(method)!=0;
}

bool spectrumMethod(const std::string& method) {
    static const std::set<std::string> methods{"dmd","dmdc","edmdc","hankel-dmd","hdmd","edmd",
        "trajectory-eigenvalue-optimization","eigenvalue_optimization"};
    return methods.count(method)!=0;
}

std::size_t odeDimension(const std::string& identifier) {
    for (const auto& system : knownSystems())
        if (system.identifier == identifier && system.family == "ode")
            return system.stateDimension;
    return 0;
}

}  // namespace

void enrichArchiveCatalogItem(CatalogItem& item,
                              const dynsyslab::core::NpzIndex& archive) {
    try {
        const auto& metadataText=archive.metadataJson;
        if (metadataText.empty()) return;
        const auto metadata=json::parse(metadataText);
        const auto label=metadata.at("label").string();
        if (!label.empty()) item.label=label;
        const auto system=systemLabel(metadata);
        if (!system.empty()) item.system=system;
        const auto method=metadata.at("method").string(metadata.at("analysis").string());
        const auto format=metadata.at("format").string();
        if (!method.empty()) {
            item.properties["method"]=method;
            appendUnique(item.tags,method);
        }
        if (!format.empty()) item.properties["archive schema"]=format;
        const std::pair<const char*,const char*> properties[]{
            {"system_preset","system preset"},{"equations_version","equations version"},
            {"equations_fingerprint","equations fingerprint"},
            {"integration_form","coordinates"},{"krom_mode_count","KROM modes"},
            {"eigenfunction_representation","eigenfunction representation"},
            {"control_map_time_domain","control-map time domain"},
            {"control_coordinate_system","control coordinates"}};
        for (const auto& property : properties) {
            const auto& value=metadata.at(property.first);
            const auto text=value.string();
            if (!text.empty()) item.properties[property.second]=text;
            else if (std::holds_alternative<double>(value.storage()))
                item.properties[property.second]=std::to_string(
                    static_cast<std::size_t>(value.number()));
        }
        const auto& legacyPreset=metadata.at("system_parameters").at("preset");
        const auto presetKey=legacyPreset.at("key").string();
        const auto presetDescription=legacyPreset.at("description").string();
        if (item.properties.count("system preset")==0 && !presetKey.empty())
            item.properties["system preset"]=presetKey;
        if (!presetDescription.empty()) item.properties["physics preset"]=presetDescription;
        else if (!presetKey.empty()) item.properties["physics preset"]=presetKey;
        const auto gas=legacyPreset.at("gas").string();
        if (!gas.empty()) item.properties["gas"]=gas;
        const auto temperature=legacyPreset.at("temperature_c").isNull()
            ? legacyPreset.at("temperature") : legacyPreset.at("temperature_c");
        if (std::holds_alternative<double>(temperature.storage()))
            item.properties["temperature"]=
                std::to_string(static_cast<int>(temperature.number()))+" C";
        const auto& controllerMetadata=metadata.at("controller");
        const auto controller=controllerMetadata.at("controller").string(
            metadata.at("controller_kind").string());
        const auto representation=controllerMetadata.at("control_representation").string();
        const auto target=controllerMetadata.at("target").string();
        if (!controller.empty()) item.properties["controller"]=controller;
        if (!representation.empty()) item.properties["control representation"]=representation;
        if (!target.empty()) item.properties["control target"]=target;
        const auto targetSpace=controllerMetadata.at("target_space").string();
        if (!targetSpace.empty()) item.properties["control target space"]=targetSpace;
        const auto dataRole=metadata.at("dataset_role").string(
            metadata.at("data_role").string());
        if (!dataRole.empty()) item.properties["data role"]=dataRole;
        const auto forcingType=metadata.at("forcing_type").string(
            metadata.at("forcing").at("forcing_type").string(
                metadata.at("training").at("kind").string()));
        if (!forcingType.empty()) item.properties["forcing type"]=forcingType;
        const auto initialPolicy=metadata.at("initial_state_policy").string();
        if (!initialPolicy.empty()) item.properties["IC policy"]=initialPolicy;
        const auto controlMaps=joined(metadata.at("control_maps"));
        if (!controlMaps.empty()) item.properties["control maps"]=controlMaps;
        const auto selected=joined(metadata.at("selected_modes"));
        if (!selected.empty()) item.properties["source mode indices"]=selected;
        if (metadata.at("independent_archive").boolean()) {
            item.properties["parent required at runtime"]="no";
            appendUnique(item.tags,"independent");
        }

        const bool analyticKrom=archive.contains("coefficients")
            && archive.contains("polynomial_exponents");
        const bool legacyKrom=archive.contains("coefficients") && archive.contains("eigenvalues")
            && !metadata.at("analysis_metadata").at("krom").isNull();
        const bool scatteredKrom=archive.contains("training_states")
            && archive.contains("eigenfunction_values")
            && archive.contains("continuous_eigenvalues")
            && metadata.at("eigenfunction_representation").string()
                =="scattered-interpolation";
        const bool dashboardReady=archive.contains("training_states")
            && ((analyticKrom && archive.contains("continuous_eigenvalues")
                 && archive.contains("state_center") && archive.contains("state_scale"))
                || legacyKrom || scatteredKrom)
            && (archive.contains("discrete_eigenvalues") || archive.contains("eigenvalues")
                || metadata.at("delta_t").number()>0.0);
        if(method=="linear-systems"
            || metadata.at("schema").string()=="dynsyslab.linear-systems/v1") {
            if(item.parentId.empty())item.category="Linear systems";
            item.kind="linear_systems";
            item.subtitle="Classical LTI analysis: gramians, balanced truncation, norms, transfer functions, loop margins";
            item.views.clear();
            for(const char* view:{"hankel_spectrum","gramian_spectrum","frequency_response","pole_zero"})
                appendUnique(item.views,view);
            if(archive.contains("bode_magnitude_db"))appendUnique(item.views,"bode");
            if(std::any_of(archive.arrays.begin(),archive.arrays.end(),[](const std::string& a){return a.rfind("loop__",0)==0;}))
                for(const char* view:{"nyquist","loop_bode","sensitivity"})appendUnique(item.views,view);
            item.properties["time domain"]=metadata.at("time_domain").string();
            item.properties["reduced order"]=std::to_string(static_cast<std::size_t>(
                metadata.at("reduction").at("order").number()));
        } else if(method=="sindy-dynamics") {
            if(item.parentId.empty())item.category="Post-processing";
            item.kind="sindy_dynamics";item.subtitle="Identified autonomous vector field";
            item.views.clear();appendUnique(item.views,"sindy_dynamics");appendUnique(item.views,"sindy_dynamics_transient");
        } else if (method=="koopman-validation" || method=="koopman_validation") {
            if (item.parentId.empty()) item.category="Validation";
            item.kind="koopman_validation";
            item.subtitle="Free-trajectory Koopman semigroup validation";
            item.views.clear();appendUnique(item.views,"koopman_validation");
            const auto role=metadata.at("dataset_role").string();
            appendUnique(item.tags,role);
            if(!role.empty())item.properties["dataset role"]=role;
        } else if(archive.contains("spectral_krom_contract_json")) {
            if(item.parentId.empty())item.category="Koopman spectral models";
            item.kind="koopman_spectral_model";
            item.subtitle="Frozen dictionary operator and spectral diagnostics";
            item.views.clear();appendUnique(item.views,"koopman_spectrum");
            item.properties["control-ready"]="no inferred physical provider or input maps";
        } else if (kromMethod(method) || analyticKrom || legacyKrom
                   || scatteredKrom) {
            if (item.parentId.empty()) item.category="KROMs";
            item.kind="krom"; item.subtitle="Native Koopman reduced-order model";
            item.views.erase(std::remove(item.views.begin(), item.views.end(),
                "data_overview"), item.views.end());
            if (dashboardReady) prependUnique(item.views,"krom_dashboard");
            appendUnique(item.views,"koopman_spectrum");
            if (hasModeWeights(archive)) appendUnique(item.views,"mode_weights");
        } else if (spectrumMethod(method) || archive.contains("continuous_eigenvalues")
                   || archive.contains("eigenvalues")) {
            if (item.parentId.empty()) item.category="Koopman spectra";
            item.kind="koopman_spectrum";
            item.subtitle="Koopman spectral result (not a control-ready KROM)";
            item.views.clear(); appendUnique(item.views,"koopman_spectrum");
            if (hasModeWeights(archive)) appendUnique(item.views,"mode_weights");
            if(method=="dmdc"||method=="edmdc") {
                item.subtitle=method=="edmdc"?"Polynomial EDMDc input/output model (not a control-ready KROM)"
                    :"DMDc input/output model (not a control-ready KROM)";
                item.properties["snapshot array"]=metadata.at("snapshot_array").string();
                item.properties["eigenfunction coordinates"]=method=="edmdc"
                    ?"Polynomials of supplied snapshots (saved normalization)":"Linear functions of supplied snapshots";
            }
            if(method=="trajectory-eigenvalue-optimization"
                ||method=="eigenvalue_optimization"){
                item.subtitle="Optimized Koopman spectrum (eigenvalues only)";
                prependUnique(item.views,"eigenvalue_optimization");
                item.properties["eigenfunctions constructed"]="no";
            }
        } else if (method=="observability" || method=="klqe"
                   || archive.contains("observability_singular_values")) {
            if (item.parentId.empty()) item.category="Estimation";
            item.kind=method=="klqe"?"estimation":"observability";
            item.subtitle=method=="klqe"?"Native Koopman state estimate":
                "Native observability analysis";
            appendUnique(item.views,"observability");
            if (method=="klqe") { appendUnique(item.views,"time"); appendUnique(item.views,"phase"); }
        } else if (decomposition(method) || archive.contains("singular_values")
                   || archive.contains("low_rank")) {
            if (item.parentId.empty()) item.category="Post-processing";
            item.kind="decomposition"; item.subtitle="Native matrix decomposition";
            appendUnique(item.views,"data_overview");
        }
        if (archive.contains("control_error_axis0")) {
            appendUnique(item.views,"control_map_error");
            item.properties["B_psi error maps"]="available";
        }
        if(method=="krom-control-search"||archive.contains("candidate_selection_score")){
            if(item.parentId.empty())item.category="Control campaigns";
            item.kind="krom_control_campaign";item.subtitle="Bounded native KROM control search";
            item.views.clear();appendUnique(item.views,"krom_control_campaign");
            const auto attempts=metadata.at("attempt_count").number();
            item.properties["attempts"]=std::to_string(static_cast<std::size_t>(attempts));
            item.properties["stop reason"]=metadata.at("stop_reason").string();
        }
        const bool hasControlHistory = archive.contains("control_values") ||
            archive.contains("controls") || archive.contains("controlled_controls") ||
            archive.contains("controlled_controls_sample_aligned");
        if (hasControlHistory) {
            appendUnique(item.tags,"forced"); appendUnique(item.views,"control");
            appendUnique(item.views,"control_fourier");
            const auto controllerKind=metadata.at("controller_kind").string();
            appendUnique(item.tags,controllerKind);
            if (!controllerKind.empty()) item.properties["controller"]=controllerKind;
        }
        if (lower(controller).find("lqr")!=std::string::npos) {
            if (item.parentId.empty()) item.category="Control studies";
            item.kind="control_attempt";
            item.subtitle="Koopman LQR control attempt";
            appendUnique(item.tags,"KLQR");
        }
        const std::size_t observableDimension = odeDimension(metadata.at("system").string());
        const bool hasTrajectory = archive.contains("values") || archive.contains("states")
            || archive.contains("controlled_states") || archive.contains("free_training_values");
        if (observableDimension && hasTrajectory) {
            appendUnique(item.views, "time");
            if (observableDimension >= 2) appendUnique(item.views, "phase");
            if (observableDimension >= 3) appendUnique(item.views, "three_d");
            appendUnique(item.views, "fourier");
            appendUnique(item.views, "time_frequency");
            appendUnique(item.views, "observable_time");
            if (observableDimension == 2) appendUnique(item.views, "observable_heatmap");
        }
        if (archive.contains("measurement_values") && archive.contains("measurement_estimate")
            && archive.contains("estimator_time"))
            appendUnique(item.views, "estimator");
        if (archive.contains("lattice_initial") && archive.contains("lattice_final"))
            appendUnique(item.views, "lattice");
        relation(item,metadata.at("source_archive").string(),
            method=="krom-extract"?"extracted from":"derived from");
        for(const auto& sourceArchive:metadata.at("source_archives").array())
            if(sourceArchive.string()!=metadata.at("source_archive").string())
                relation(item,sourceArchive.string(),"fit from trajectory");
        relation(item,metadata.at("source_spectrum").string(),"initialized from spectrum");
        const auto sourceKrom=metadata.at("source_krom").string();
        if (sourceKrom!=metadata.at("source_archive").string())
            relation(item,sourceKrom,
                method=="koopman-validation"?"evaluates KROM":"uses KROM");
        relation(item,metadata.at("controller_configuration").at("settings")
            .at("krom").string(),"controlled by KROM");
        relation(item,controllerMetadata.at("rom").at("source_archive").string(),
            "controlled by KROM");
        relation(item,metadata.at("training_preflight").at("source_krom").string(),
            "expands training domain of");
        for(const auto& result:metadata.at("result_archives").array())
            relation(item,result.string(),"produced control result");
    } catch (...) {
        // A malformed or unsupported archive remains visible through the
        // conservative filesystem classification and can still be inspected.
    }
}

void enrichArchiveCatalogItem(CatalogItem& item,
                              const std::filesystem::path& source) {
    try { enrichArchiveCatalogItem(item, dynsyslab::core::readNpzIndex(source)); }
    catch (...) {
        // The filesystem item remains visible even if metadata cannot be read.
    }
}

void resolveArchiveCatalogRelations(Catalog& catalog,
                                    const std::filesystem::path& projectRoot) {
    std::map<std::string,std::string> bySource;
    std::map<std::string,std::size_t> byId;
    for (std::size_t index=0;index<catalog.items.size();++index) {
        const auto& item=catalog.items[index];byId[item.id]=index;
        auto source=std::filesystem::u8path(item.source);
        if (source.is_relative()) source=projectRoot/source;
        bySource.emplace(normalizedSource(source),item.id);
    }
    for (auto& item : catalog.items) {
        auto output=item.relations.begin();
        for (auto relation=item.relations.begin();relation!=item.relations.end();++relation) {
            if (relation->target.rfind(kSourcePrefix,0)==0) {
                const auto found=bySource.find(relation->target.substr(
                    std::char_traits<char>::length(kSourcePrefix)));
                if (found==bySource.end()) continue;
                relation->target=found->second;
            }
            *output++=*relation;
        }
        item.relations.erase(output,item.relations.end());
    }
    struct Reverse {std::size_t target;std::string source,kind,label;};
    std::vector<Reverse> reverse;
    for (auto& item:catalog.items) for (const auto& link:item.relations) {
        const auto target=byId.find(link.target);if(target==byId.end())continue;
        reverse.push_back({target->second,item.id,reverseKind(link.kind),item.label});
        if (item.parentId.empty() && item.kind=="krom" && link.kind=="extracted from"
            && catalog.items[target->second].kind=="krom") item.parentId=link.target;
    }
    for (const auto& link:reverse)
        catalogRelation(catalog.items[link.target],link.source,link.kind,link.label);
}

}  // namespace dsl
