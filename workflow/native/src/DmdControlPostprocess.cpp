#include "PostprocessInternal.hpp"
#include <dynsyslab/core/SpectrumArchive.hpp>

#include <dynsyslab/core/DmdControl.hpp>
#include <dynsyslab/core/EdmdControl.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
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
    auto times=time.doubles();
    if(!std::all_of(times.begin(),times.end(),[](double t){return std::isfinite(t);}))
        throw std::invalid_argument("DMDc time must be finite");
    const std::size_t samples=time.shape[0];
    source.record.series=core::TimeSeries(times,data.doubles(),data.shape.size()==1?1:data.shape[1]);
    const bool interval=archive.contains("control_interval_values");
    if(!interval&&!archive.contains("control_values"))
        throw std::invalid_argument("DMDc requires recorded control_interval_values or control_values");
    const auto& input=archive.at(interval?"control_interval_values":"control_values");
    if(input.shape.empty()||input.shape.size()>2
        ||input.shape[0]!=samples-(interval?1:0))
        throw std::invalid_argument("DMDc inputs must have N-1 interval rows or N sample rows");
    const std::size_t inputDimension=input.shape.size()==1?1:input.shape[1];
    auto inputs=input.doubles();
    inputs.resize((samples-1)*inputDimension);  // one row per interval
    source.record.inputs=inputSeries(times,std::move(inputs),inputDimension,true);
    source.record.inputsOnIntervals=true;
    return source;
}
} // namespace

std::filesystem::path runDmdControl(const std::filesystem::path& path,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label,bool extended) {
    const auto input=core::readNpz(path);
    const auto key=settings.at("snapshot_array").string(input.contains("values")?"values":"states");
    const auto source=measuredSource(path,input,key);
    const auto n=source.dimension(),m=source.inputDimension();
    core::DmdControlOptions options;
    options.rank=rankSetting(settings,"rank");
    options.augmentedRank=rankSetting(settings,"augmented_rank");
    options.singularValueTolerance=settings.at("svd_tolerance").number(0);
    const double dt=uniformDeltaTime(source);
    // Tighter relative test than the legacy helper for very small physical dt.
    for(std::size_t i=1;i<source.samples();++i)
        if(source.time()[i]<=source.time()[i-1]||std::abs((source.time()[i]-source.time()[i-1])-dt)
            >1e-7*dt+16*std::numeric_limits<double>::epsilon()*std::abs(source.time()[i]))
            throw std::invalid_argument("DMDc requires uniform positive sample intervals");
    // Optional pair selection: time_window [t0, t1] and state_window [[lo, hi], ...] per snapshot
    // component. A pair is kept only when both of its snapshots lie inside both windows, so a fit can be
    // restricted to the neighbourhood of an operating point (for example an unstable target).
    // time_window may also name two of the source's metadata keys, {"from_metadata": ["drive_start",
    // "drive_stop"]}, so each record of a bank is fitted on its own active (e.g. forcing) interval.
    const auto& stateWindow=settings.at("state_window");
    if(!stateWindow.isNull()&&stateWindow.array().size()!=n)
        throw std::invalid_argument("DMDc state_window needs one [low, high] pair per snapshot component");
    // Appends the windowed snapshot pairs of one record (every record must match the first's
    // snapshot and input widths and its sample interval).
    auto collect=[&](const SourceData& record,std::vector<double>& xs,std::vector<double>& ys,std::vector<double>& us) {
        if(record.dimension()!=n||record.inputDimension()!=m||std::abs(uniformDeltaTime(record)-dt)>1e-9*dt)
            throw std::invalid_argument("DMDc pooled and validation records must match the source's widths and sample interval");
        auto timeWindow=settings.at("time_window");
        if(!timeWindow.isNull()&&!timeWindow.at("from_metadata").isNull()) {
            const auto meta=json::parse(record.metadataJson);
            json::Value::Array bounds;
            for(const auto& name:timeWindow.at("from_metadata").array()) {
                const auto& value=meta.at(name.string());
                if(!std::holds_alternative<double>(value.storage()))
                    throw std::invalid_argument("DMDc time_window: the source metadata has no number "+name.string());
                bounds.push_back(value);
            }
            if(bounds.size()!=2) throw std::invalid_argument("DMDc time_window from_metadata names two keys");
            timeWindow=json::Value(std::move(bounds));
        }
        auto inside=[&](std::size_t k) {
            if(!timeWindow.isNull()&&(record.time()[k]<timeWindow.array().at(0).number()
                    ||record.time()[k]>timeWindow.array().at(1).number())) return false;
            if(!stateWindow.isNull()) for(std::size_t j=0;j<n;++j) {
                const double v=record.values()[k*n+j];
                const auto& b=stateWindow.array().at(j).array();
                if(v<b.at(0).number()||v>b.at(1).number()) return false;
            }
            return true;
        };
        for(std::size_t k=0;k+1<record.samples();++k) {
            if(!inside(k)||!inside(k+1)) continue;
            xs.insert(xs.end(),record.values().begin()+k*n,record.values().begin()+(k+1)*n);
            ys.insert(ys.end(),record.values().begin()+(k+1)*n,record.values().begin()+(k+2)*n);
            us.insert(us.end(),record.inputs().begin()+k*m,record.inputs().begin()+(k+1)*m);
        }
    };
    auto load=[&](const json::Value& item) {
        const auto file=std::filesystem::u8path(item.string());
        return measuredSource(file,core::readNpz(file),key);
    };
    std::vector<double> x,y,u;
    collect(source,x,y,u);
    // pool_sources: further records fitted together with the source (pairs never cross records).
    json::Value::Array pooled;
    for(const auto& item:settings.at("pool_sources").array()) { collect(load(item),x,y,u); pooled.push_back(item); }
    const std::size_t pairs=x.size()/n;
    if(pairs<2) throw std::invalid_argument("DMDc time_window/state_window keep fewer than two snapshot pairs");
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
            x,y,u,pairs,n,m,dt,liftedOptions);
    } else linear=core::dynamicModeDecompositionWithControl(x,y,u,pairs,n,m,dt,options);
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
    const auto& spectrum=fit.spectrum;
    core::writeSpectrumArrays(archive,spectrum);  // a zero multiplier has no finite continuous rate
    const bool continuousRates=spectrum.hasContinuousRates();
    complex("modes",{features,r},fit.modes);
    complex("exact_modes",{features,r},fit.exactModes);
    if(polynomial) {
        const auto& d=polynomial->dictionary;
        archive.arrays["polynomial_exponents"]=core::NpyArray::fromInt64s(
            {features,n},{d.exponents.begin(),d.exponents.end()});
        real("state_center",{n},d.center);
        real("state_scale",{n},d.scale);
    }
    if(fit.hasEigenfunctions) {
        complex("eigenfunction_coefficients",{r,features},fit.eigenfunctionCoefficients);
        complex("modal_input_matrix",{r,m},fit.modalInputMatrix);
        if(settings.at("save_eigenfunction_values").boolean(true)) {
            complex("eigenfunction_values",{source.samples(),r},
                polynomial?core::evaluateEdmdControlEigenfunctions(*polynomial,source.values(),source.samples())
                    :core::evaluateDmdControlEigenfunctions(fit,source.values(),source.samples()));
            real("time",{source.samples()},source.time());
        }
    }
    json::Value::Array warnings;
    if(!fit.fullRowRank) warnings.emplace_back("[X;U] is row-rank deficient: A and B are not separately identifiable in the full supplied space.");
    if(fit.augmentedRank<fit.augmentedNumericalRank) warnings.emplace_back("Augmented SVD was truncated: discarded directions are not identified.");
    if(!fit.hasEigenfunctions) warnings.emplace_back("Ill-conditioned/defective eigenbasis: modal eigenfunctions are omitted; A and B remain usable.");
    if(extended) warnings.emplace_back("Polynomial lifting does not generally preserve additive constant input coupling; validate held-out forced rollouts. No continuous B_psi is inferred.");
    json::Value::Object details{
        {"schema",extended?"dynsyslab.edmdc/v1":"dynsyslab.dmdc/v1"},{"eigenvalue_domain","discrete-time"},
        {"delta_t",dt},{"rank",static_cast<double>(r)},
        {"pairs",static_cast<double>(pairs)},{"available_pairs",static_cast<double>(source.samples()-1)},
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
        {"error_dataset_role","training; not held-out validation (see validation)"},
        {"continuous_rates_available",continuousRates},
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
    if(!pooled.empty()) details["pool_sources"]=json::Value(pooled);
    // validation_sources: one-step prediction of held-out records' windowed pairs, with the fit
    // frozen: ||y_hat - y|| / ||y|| and ||y_hat - y|| / ||y - x|| (error over the increment).
    {
        json::Value::Array scores;
        double sumError=0,sumEnergy=0,sumIncrement=0;
        for(const auto& item:settings.at("validation_sources").array()) {
            std::vector<double> vx,vy,vu;
            collect(load(item),vx,vy,vu);
            const auto count=vx.size()/n;
            if(count<1) throw std::invalid_argument("DMDc validation record keeps no snapshot pairs: "+item.string());
            const auto predicted=polynomial?core::predictEdmdControl(*polynomial,vx,vu,count):core::predictDmdControl(linear,vx,vu,count);
            double e=0,energy=0,increment=0;
            for(std::size_t i=0;i<vy.size();++i){e+=std::pow(predicted[i]-vy[i],2);energy+=vy[i]*vy[i];increment+=std::pow(vy[i]-vx[i],2);}
            sumError+=e;sumEnergy+=energy;sumIncrement+=increment;
            scores.push_back(json::Value::Object{{"source",item},{"pairs",static_cast<double>(count)},
                {"one_step_relative_l2",energy>0?std::sqrt(e/energy):0.0},{"error_over_increment_l2",increment>0?std::sqrt(e/increment):0.0}});
        }
        if(!scores.empty()) details["validation"]=json::Value::Object{{"records",json::Value(std::move(scores))},
            {"one_step_relative_l2",sumEnergy>0?std::sqrt(sumError/sumEnergy):0.0},
            {"error_over_increment_l2",sumIncrement>0?std::sqrt(sumError/sumIncrement):0.0}};
    }
    const auto text=json::dump(json::Value(std::move(details)));
    auto meta=json::parse(metadata(source,method,label,text.substr(1,text.size()-2)));
    std::optional<core::KoopmanModel> exportable;
    if(fit.hasEigenfunctions&&continuousRates) {
        const auto sourceMeta=json::parse(source.metadataJson);
        // Only a native solution's values have an unambiguous physical-state
        // contract. Arbitrary measured/latent snapshots need their own provider.
        if(key=="values"&&sourceMeta.at("schema").string()=="dynsyslab.solution/v1") {
            auto bank=polynomial?polynomial->eigenfunctionSet():linear.eigenfunctionSet();
            auto coordinates=core::FunctionDomain::physicalState(bank.functions.coordinateDimension());
            coordinates.names=bank.functions.coordinates().names;
            core::KoopmanModel model(bank.functions.withCoordinates(coordinates),bank.spectrum,systemConfiguration(source));
            const auto plan=planAutonomousKoopmanValidation(model,{},root);
            auto refs=core::koopmanDatasetReferencesFromMetadata(meta);
            for(const auto& ref:core::koopmanDatasetReferencesFromMetadata(plan))core::mergeKoopmanDatasetReference(refs,ref);
            meta.object()["dataset_references"]=core::koopmanDatasetReferencesJson(refs);
            meta.object()["autonomous_validation_plan"]=plan;
            if(settings.at("export_krom").boolean(false)) exportable.emplace(std::move(model));
        } else meta.object()["autonomous_validation_plan"]=json::Value::Object{
            {"status","coordinate-provider-required"},{"simulation_executed",false},
            {"message","Export the KROM with its system and frozen coordinate provider, then run the reuse-first autonomous validation lookup. Do not assume measured snapshots are physical states."}};
    }
    archive.arrays["metadata_json"]=core::NpyArray::fromString(json::dump(meta));
    const auto fitPath=publishObject(root,source,dynsyslab::core::object_kind::dynamicalModel,"koopman-"+method,std::move(archive));  // the controlled linear fit
    if(!settings.at("export_krom").boolean(false)) return fitPath;
    if(!exportable)
        throw std::invalid_argument("export_krom needs a fit with eigenfunctions and finite "
            "continuous rates on a native solution's state values");
    // Export the eigenfunction bank as a krom so krom_extract, input_map, KLQE and KLQR can use it. The
    // selected (windowed) training states and the eigenfunctions' values there (its samples) are
    // attached so krom_extract can refit named readouts on exactly the fitted domain.
    const auto values=polynomial?core::evaluateEdmdControlEigenfunctions(*polynomial,x,pairs)
                                :core::evaluateDmdControlEigenfunctions(linear,x,pairs);
    std::map<std::string,core::NpyArray> extra;
    extra["training_states"]=core::NpyArray::fromDoubles({pairs,n},x);
    {
        core::EigenfunctionSamples samples;
        samples.modes=r;
        samples.values=values;
        for(std::size_t k=0;k<pairs;++k){samples.sampleIndices.push_back(k);samples.time.push_back(double(k));}
        exportable->setSamples(std::move(samples));
    }
    // Per-mode consistency on the fitted pairs: ||psi_j(x+) - mu_j psi_j(x) - b_j u|| / ||psi_j(x+) - psi_j(x)||
    // (ResDMD-style relative residual of the controlled model; Colbrook & Townsend, Comm. Pure Appl.
    // Math. 77, 2024, arXiv:2111.14889). krom_extract's principal rule uses it.
    const auto next=polynomial?core::evaluateEdmdControlEigenfunctions(*polynomial,y,pairs)
                              :core::evaluateDmdControlEigenfunctions(linear,y,pairs);
    std::vector<double> consistency(r,0.0);
    for(std::size_t j=0;j<r;++j) {
        double num=0,den=0;
        for(std::size_t k=0;k<pairs;++k) {
            std::complex<double> drive=0; for(std::size_t i=0;i<m;++i) drive+=fit.modalInputMatrix[j*m+i]*u[k*m+i];
            num+=std::norm(next[k*r+j]-fit.spectrum.values()[j]*values[k*r+j]-drive);
            den+=std::norm(next[k*r+j]-values[k*r+j]);
        }
        consistency[j]=den>0?std::sqrt(num/den):0.0;
    }
    extra["mode_consistency"]=core::NpyArray::fromDoubles({r},consistency);
    meta.object()["edmdc_fit_archive"]=fitPath.generic_u8string();
    meta.object()["extraction_algorithm"]=method;
    return publishKrom(root,source,*exportable,"koopman-"+method+"-krom",json::dump(meta),std::move(extra));
}
} // namespace dynsyslab::workflow::postprocess_detail
