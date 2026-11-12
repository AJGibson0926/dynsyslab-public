#include "PostprocessInternal.hpp"

#include <dynsyslab/core/DmdControl.hpp>
#include <dynsyslab/core/EdmdControl.hpp>
#include <dynsyslab/workflow/KromValidationPlan.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace core=dynsyslab::core;

std::size_t rankSetting(const json::Value& settings,const char* name) {
    const auto& value=settings.at(name);
    const double n=value.number(-1);
    if(value.isNull()) return 0;
    if(!std::isfinite(n)||n<0||n!=std::floor(n)||n>1000000000)
        throw std::invalid_argument(std::string("DMDc ")+name+" must be a nonnegative integer");
    return static_cast<std::size_t>(n);
}

// Retain the selected array's stored units. In particular, never interpret
// two latent/measurement columns as r,v or call a system model to transform them.
SourceData measuredSource(const std::filesystem::path& path,
    const core::NpzRecord& archive,const std::string& snapshotKey) {
    SourceData source; source.path=path; source.metadataJson=archive.metadataJson();
    const auto& time=archive.at("time");
    const auto& data=archive.at(snapshotKey);
    if(time.shape.size()!=1||time.shape[0]<3||data.shape.empty()
        ||data.shape.size()>2||data.shape[0]!=time.shape[0])
        throw std::invalid_argument("DMDc requires time[N] and selected snapshots[N,d] (or [N])");
    source.samples=time.shape[0]; source.dimension=data.shape.size()==1?1:data.shape[1];
    source.time=time.doubles(); source.values=data.doubles();
    if(!std::all_of(source.time.begin(),source.time.end(),[](double t){return std::isfinite(t);}))
        throw std::invalid_argument("DMDc time must be finite");
    const bool interval=archive.contains("control_interval_values");
    if(!interval&&!archive.contains("control_values"))
        throw std::invalid_argument("DMDc requires recorded control_interval_values or control_values");
    const auto& input=archive.at(interval?"control_interval_values":"control_values");
    if(input.shape.empty()||input.shape.size()>2
        ||input.shape[0]!=source.samples-(interval?1:0))
        throw std::invalid_argument("DMDc inputs must have N-1 interval rows or N sample rows");
    source.inputDimension=input.shape.size()==1?1:input.shape[1];
    source.inputs=input.doubles();
    source.inputs.resize((source.samples-1)*source.inputDimension);
    return source;
}
} // namespace

std::filesystem::path runDmdControl(const std::filesystem::path& path,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label,bool extended) {
    const auto input=core::readNpz(path);
    const auto key=settings.at("snapshot_array").string(input.contains("values")?"values":"states");
    const auto source=measuredSource(path,input,key);
    const auto pairs=source.samples-1,n=source.dimension,m=source.inputDimension;
    core::DmdControlOptions options;
    options.rank=rankSetting(settings,"rank");
    options.augmentedRank=rankSetting(settings,"augmented_rank");
    options.singularValueTolerance=settings.at("svd_tolerance").number(0);
    const double dt=uniformDeltaTime(source);
    // Tighter relative test than the legacy helper for very small physical dt.
    for(std::size_t i=1;i<source.samples;++i)
        if(source.time[i]<=source.time[i-1]||std::abs((source.time[i]-source.time[i-1])-dt)
            >1e-7*dt+16*std::numeric_limits<double>::epsilon()*std::abs(source.time[i]))
            throw std::invalid_argument("DMDc requires uniform positive sample intervals");
    const std::vector<double> x(source.values.begin(),source.values.begin()+pairs*n);
    const std::vector<double> y(source.values.begin()+n,source.values.end());
    std::optional<core::EdmdControlResult> polynomial;
    core::DmdControlResult linear;
    if(extended) {
        if(settings.at("dictionary").string("polynomial")!="polynomial")
            throw std::invalid_argument("EDMDc currently supports the SINDy polynomial dictionary only");
        core::EdmdControlOptions liftedOptions;
        liftedOptions.regression=options;
        if(!settings.at("polynomial_degree").isNull())
            liftedOptions.polynomialDegree=rankSetting(settings,"polynomial_degree");
        liftedOptions.includeConstant=settings.at("include_constant").boolean(true);
        liftedOptions.normalize=settings.at("normalize").boolean(true);
        polynomial=core::extendedDynamicModeDecompositionWithControl(
            x,y,source.inputs,pairs,n,m,dt,liftedOptions);
    } else linear=core::dynamicModeDecompositionWithControl(x,y,source.inputs,pairs,n,m,dt,options);
    const auto& fit=polynomial?polynomial->lifted:linear;
    const auto r=fit.rank,features=fit.dimension;
    const std::string method=extended?"edmdc":"dmdc";
    core::NpzRecord archive;
    auto real=[&](const char* name,std::vector<std::size_t> shape,const std::vector<double>& v) {
        archive.arrays[name]=core::NpyArray::fromDoubles(std::move(shape),v);
    };
    auto complex=[&](const char* name,std::vector<std::size_t> shape,
                     const std::vector<std::complex<double>>& v) {
        archive.arrays[name]=core::NpyArray::fromComplex(std::move(shape),v);
    };
    real("operator",{r,r},fit.operatorMatrix);
    real("control_matrix",{r,m},fit.controlMatrix);
    real("projection_basis",{features,r},fit.basis);
    real("augmented_singular_values",{fit.augmentedSingularValues.size()},fit.augmentedSingularValues);
    real("output_singular_values",{fit.outputSingularValues.size()},fit.outputSingularValues);
    complex("eigenvalues",{r},fit.discreteEigenvalues);
    complex("modes",{features,r},fit.modes);
    complex("exact_modes",{features,r},fit.exactModes);
    if(polynomial) {
        const auto& d=polynomial->dictionary;
        archive.arrays["polynomial_exponents"]=core::NpyArray::fromInt64s(
            {features,n},{d.exponents.begin(),d.exponents.end()});
        real("state_center",{n},d.center);
        real("state_scale",{n},d.scale);
    }
    // A zero multiplier has no finite logarithmic continuous-time rate.
    std::vector<std::complex<double>> continuous;
    for(const auto mu:fit.discreteEigenvalues) {
        if(std::abs(mu)==0) { continuous.clear(); break; }
        continuous.push_back(std::log(mu)/dt);
    }
    if(continuous.size()==r) complex("continuous_eigenvalues",{r},continuous);
    if(fit.hasEigenfunctions) {
        complex("eigenfunction_coefficients",{r,features},fit.eigenfunctionCoefficients);
        complex("modal_input_matrix",{r,m},fit.modalInputMatrix);
        if(settings.at("save_eigenfunction_values").boolean(true)) {
            complex("eigenfunction_values",{source.samples,r},
                polynomial?core::evaluateEdmdControlEigenfunctions(*polynomial,source.values,source.samples)
                    :core::evaluateDmdControlEigenfunctions(fit,source.values,source.samples));
            real("time",{source.samples},source.time);
        }
    }
    json::Value::Array warnings;
    if(!fit.fullRowRank) warnings.emplace_back("[X;U] is row-rank deficient: A and B are not separately identifiable in the full supplied space.");
    if(fit.augmentedRank<fit.augmentedNumericalRank) warnings.emplace_back("Augmented SVD was truncated: discarded directions are not identified.");
    if(!fit.hasEigenfunctions) warnings.emplace_back("Ill-conditioned/defective eigenbasis: modal eigenfunctions are omitted; A and B remain usable.");
    if(extended) warnings.emplace_back("Polynomial lifting does not generally preserve additive constant input coupling; validate held-out forced rollouts. No continuous B_psi is inferred.");
    json::Value::Object details{
        {"model_schema",extended?"dynsyslab.edmdc/v1":"dynsyslab.dmdc/v1"},{"eigenvalue_domain","discrete-time"},
        {"delta_t",dt},{"rank",static_cast<double>(r)},
        {"pairs",static_cast<double>(pairs)},
        {"augmented_rank",static_cast<double>(fit.augmentedRank)},
        {"augmented_numerical_rank",static_cast<double>(fit.augmentedNumericalRank)},
        {"snapshot_dimension",static_cast<double>(n)},{"input_dimension",static_cast<double>(m)},
        {"feature_count",static_cast<double>(features)},
        {"snapshot_array",key},{"input_array",input.contains("control_interval_values")?"control_interval_values":"control_values"},
        {"input_alignment","row k drives snapshot k to k+1; final sample input unused"},
        {"coordinate_policy",extended?"selected source array; normalization fitted on X only and frozen; no state/EOM access":"selected source array in stored units; no state/EOM access or centering"},
        {"eigenfunction_representation",extended?"analytic-polynomial":"linear-projected-dmdc"},
        {"eigenfunction_independent_variable","selected source snapshot vector q"},
        {"eigenfunctions_available",fit.hasEigenfunctions},{"full_row_rank",fit.fullRowRank},
        {"one_step_relative_error",fit.oneStepRelativeError},
        {"projection_relative_error",fit.projectionRelativeError},
        {"retained_condition_number",fit.retainedConditionNumber},
        {"eigenvector_condition_number",fit.eigenvectorConditionNumber},
        {"error_dataset_role","training; not held-out validation"},
        {"continuous_rates_available",continuous.size()==r},
        {"continuous_rate_convention","principal log(mu)/dt; not a continuous input map"},
        {"settings",settings},{"warnings",std::move(warnings)}};
    if(polynomial) {
        details["dictionary"]="polynomial";
        details["polynomial_degree"]=static_cast<double>(polynomial->dictionary.degree);
        details["include_constant"]=polynomial->dictionary.includeConstant;
        details["normalize"]=settings.at("normalize").boolean(true);
        details["dictionary_order"]="SINDy canonical graded total-degree monomials";
        details["operator_coordinates"]="z = P^T Theta(q); modes reconstruct dictionary features";
        details["one_step_error_coordinates"]="normalized polynomial features";
        details["snapshot_one_step_relative_error"]=polynomial->snapshotOneStepRelativeError;
        details["snapshot_decoder"]="center + scale * predicted degree-one dictionary entries";
    }
    const auto text=json::dump(json::Value(std::move(details)));
    auto meta=json::parse(metadata(source,method,label,text.substr(1,text.size()-2)));
    if(polynomial&&fit.hasEigenfunctions&&continuous.size()==r) {
        const auto sourceMeta=json::parse(source.metadataJson);
        // Only a native solution's values have an unambiguous physical-state
        // contract. Arbitrary measured/latent snapshots need their own provider.
        if(key=="values"&&sourceMeta.at("format").string()=="dynsyslab.solution/v1") {
            auto bank=polynomial->eigenfunctionSet();auto coordinates=bank.coordinates();
            coordinates.kind=core::EigenfunctionCoordinateKind::physicalState;
            core::KoopmanModel model(bank.withCoordinates(coordinates),systemConfiguration(source));
            const auto plan=planAutonomousKoopmanValidation(model,{},root);
            auto refs=core::koopmanDatasetReferencesFromMetadata(meta);
            for(const auto& ref:core::koopmanDatasetReferencesFromMetadata(plan))core::mergeKoopmanDatasetReference(refs,ref);
            meta.object()["dataset_references"]=core::koopmanDatasetReferencesJson(refs);
            meta.object()["autonomous_validation_plan"]=plan;
        } else meta.object()["autonomous_validation_plan"]=json::Value::Object{
            {"status","coordinate-provider-required"},{"simulation_executed",false},
            {"message","Export the KROM with its system and frozen coordinate provider, then run the reuse-first autonomous validation lookup. Do not assume measured snapshots are physical states."}};
    }
    archive.arrays["metadata_json"]=core::NpyArray::fromString(json::dump(meta));
    return publish(root,source,"koopman-"+method,label,std::move(archive));
}
} // namespace dynsyslab::workflow::postprocess_detail
