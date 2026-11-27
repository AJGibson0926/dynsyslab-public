#include "LinearAnalysisModal.hpp"
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/core/EigenfunctionObservability.hpp>
#include <complex>
#include <dynsyslab/core/Controllability.hpp>
#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <cmath>
#include <optional>
#include <stdexcept>
namespace dynsyslab::workflow::postprocess_detail {
namespace {
using namespace dynsyslab::core;
ObservabilityOptions options(const json::Value& s) {
    ObservabilityOptions o;
    const auto method=s.at("analysis_method").string("continuous-gramian");
    if(method=="matrix")o.method=ObservabilityMethod::classicalMatrix;
    else if(method=="continuous-gramian")o.method=ObservabilityMethod::continuousFiniteGramian;
    else if(method=="continuous-infinite-gramian")o.method=ObservabilityMethod::continuousInfiniteGramian;
    else if(method=="discrete-gramian")o.method=ObservabilityMethod::discreteGramian;
    else if(method=="discrete-infinite-gramian")o.method=ObservabilityMethod::discreteInfiniteGramian;
    else throw std::invalid_argument("unknown linear analysis method");
    const double blocks=s.at("blocks").number();
    if(!std::isfinite(blocks)||blocks<0||std::floor(blocks)!=blocks||blocks>1e8)
        throw std::invalid_argument("analysis blocks must be a bounded nonnegative integer");
    o.blockCount=static_cast<std::size_t>(blocks);o.horizon=s.at("horizon").number(1);
    o.sampleTime=s.at("sample_time").number();o.relativeTolerance=s.at("relative_tolerance").number(1e-10);
    return o;
}
void spectralArrays(NpzRecord& record,const std::string& prefix,std::size_t n,
    const std::vector<double>& gramian,const std::vector<double>& strengths,
    const std::vector<double>& eigenvalues,const std::vector<double>& directions) {
    record.arrays[prefix+"_gramian"]=NpyArray::fromDoubles({n,n},gramian);
    record.arrays[prefix+"_singular_values"]=NpyArray::fromDoubles({strengths.size()},strengths);
    record.arrays[prefix+"_gramian_eigenvalues"]=NpyArray::fromDoubles({eigenvalues.size()},eigenvalues);
    record.arrays[prefix+"_directions"]=NpyArray::fromDoubles({n,n},directions);
}
json::Value finite(double value){return std::isfinite(value)?json::Value(value):json::Value{};}
template<class Result> json::Value summary(const Result& r) {
    return json::Value::Object{{"rank",double(r.rank)},{"state_dimension",double(r.stateDimension)},
        {"rank_tolerance",r.rankTolerance},{"condition_number",finite(r.conditionNumber)},
        {"resolved_condition_number",finite(r.resolvedConditionNumber)},
        {"effective_rank",r.effectiveRank},{"horizon",finite(r.horizon)}};
}
}
// analysis_method "modal-matrix": obsv in discrete modal coordinates, one complex coordinate
// per eigenfunction (M = diag(exp(lambda dt))), with its direct SVD.
std::filesystem::path runModalObservabilityMatrix(const std::filesystem::path& sourcePath,const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    using namespace dynsyslab::core;
    std::vector<std::size_t> selected,rows;
    for(const auto& v:settings.at("selected_modes").array())selected.push_back(static_cast<std::size_t>(v.number()));
    for(const auto& v:settings.at("output_rows").array())rows.push_back(static_cast<std::size_t>(v.number()));
    const auto model=loadKoopmanModel(sourcePath,selected);
    const auto measurement=settings.at("measurement").string();
    if(measurement.empty())throw std::invalid_argument("modal-matrix needs a measurement");
    const double blocks=settings.at("blocks").number(0);
    if(!(blocks>=0)||std::floor(blocks)!=blocks)throw std::invalid_argument("blocks must be a nonnegative integer");
    const auto r=modalObservabilityMatrix(model,measurement,settings.at("sample_time").number(0),
                                          static_cast<std::size_t>(blocks),rows);
    const auto complexArray=[](const std::vector<std::complex<double>>& v,std::size_t a,std::size_t b) {
        std::vector<double> flat;for(const auto& x:v){flat.push_back(x.real());flat.push_back(x.imag());}
        return NpyArray::fromDoubles({a,b,2},flat);
    };
    const std::size_t n=r.modes,k=r.singularValues.size();
    NpzRecord out;
    out.arrays["modal_observability_matrix"]=complexArray(r.matrix,r.blocks*r.outputs,n);
    out.arrays["modal_singular_values"]=NpyArray::fromDoubles({k},r.singularValues);
    out.arrays["modal_directions"]=complexArray(r.directions,n,k);
    std::vector<double> magnitude;for(const auto& x:r.directions)magnitude.push_back(std::abs(x));
    out.arrays["modal_direction_magnitudes"]=NpyArray::fromDoubles({n,k},magnitude);
    out.arrays["modal_multipliers"]=complexArray(r.multipliers,n,1);
    const auto lambda=model.continuousEigenvalues();
    out.arrays["continuous_eigenvalues"]=complexArray(lambda,n,1);
    json::Value::Array sv;for(double s:r.singularValues)sv.push_back(s);
    std::size_t rank=0;for(double s:r.singularValues)rank+=s>1e-10*r.singularValues.front();
    json::Value::Object details{{"schema","dynsyslab.linear-analysis/v1"},{"analysis_method","modal-matrix"},
        {"source_krom",sourcePath.generic_u8string()},{"measurement",measurement},{"modes",double(n)},
        {"blocks",double(r.blocks)},{"sample_time",r.sampleTime},{"singular_values",sv},{"numerical_rank",double(rank)},
        {"definition","O = [C; C M; ...; C M^(blocks-1)], M = diag(exp(lambda dt)), complex modal coordinates; direct SVD"}};
    auto serialized=json::dump(json::Value(details));
    out.arrays["metadata_json"]=NpyArray::fromString(metadata(source,"linear-analysis",label,
        serialized.substr(1,serialized.size()-2)));
    return publishObject(root, source, dynsyslab::core::object_kind::report, "linear-analysis",std::move(out));
}

// analysis_method "eigenfunction-observability": how well each eigenfunction can be determined
// from the measurement over a window (EigenfunctionObservability.hpp).
std::filesystem::path runEigenfunctionObservability(const std::filesystem::path& sourcePath,const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    using namespace dynsyslab::core;
    std::vector<std::size_t> selected;
    for(const auto& v:settings.at("selected_modes").array())selected.push_back(static_cast<std::size_t>(v.number()));
    const auto model=loadKoopmanModel(sourcePath,selected);
    const auto measurement=settings.at("measurement").string();
    if(measurement.empty())throw std::invalid_argument("eigenfunction-observability needs a measurement");
    EigenfunctionObservabilityOptions o;
    o.horizon=settings.at("horizon").number(1.0);
    o.noiseIntensity=settings.at("noise_intensity").number(1e-4);
    o.regularization=settings.at("regularization").number(1e-12);
    o.processCovariance=settings.at("process_covariance").number(0.0);
    const auto r=eigenfunctionObservability(model,measurement,o);
    const std::size_t n=r.modes;
    NpzRecord out;
    const auto put=[&](const char* name,const std::vector<double>& v){if(!v.empty())out.arrays[name]=NpyArray::fromDoubles({v.size()},v);};
    put("signal",r.signal);put("distinguishable",r.distinguishable);put("isolation",r.isolation);
    put("bound_std",r.boundStd);put("pair_bound_std",r.pairBoundStd);put("kalman_std",r.kalmanStd);
    std::vector<double> partner,lambda,gram;
    for(auto p:r.partner)partner.push_back(double(p));
    for(const auto& l:r.eigenvalues){lambda.push_back(l.real());lambda.push_back(l.imag());}
    for(const auto& g:r.gramian){gram.push_back(g.real());gram.push_back(g.imag());}
    put("partner",partner);
    out.arrays["continuous_eigenvalues"]=NpyArray::fromDoubles({n,2},lambda);
    out.arrays["modal_gramian"]=NpyArray::fromDoubles({n,n,2},gram);
    std::vector<double> index;for(std::size_t k=0;k<n;++k)index.push_back(double(k+1));
    put("eigenfunction_index",index);
    json::Value::Object details{{"schema","dynsyslab.linear-analysis/v1"},{"analysis_method","eigenfunction-observability"},
        {"source_krom",sourcePath.generic_u8string()},{"measurement",measurement},{"modes",double(n)},
        {"horizon",o.horizon},{"noise_intensity",o.noiseIntensity},{"regularization",o.regularization},
        {"process_covariance",o.processCovariance},
        {"definition","G = finite-horizon modal gramian over [0, horizon]; signal = sqrt(G_nn); "
                      "distinguishable = 1/sqrt((G^-1)_nn); bound_std = sqrt(V (G^-1)_nn); pair_bound_std = worst direction of the conjugate pair"}};
    auto serialized=json::dump(json::Value(details));
    out.arrays["metadata_json"]=NpyArray::fromString(metadata(source,"linear-analysis",label,
        serialized.substr(1,serialized.size()-2)));
    return publishObject(root, source, dynsyslab::core::object_kind::report, "linear-analysis",std::move(out));
}

std::filesystem::path runLinearAnalysis(const std::filesystem::path& sourcePath,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    using namespace dynsyslab::core;
    const auto sourceArchive=readNpz(sourcePath);
    SourceData source;source.path=sourcePath;source.metadataJson=sourceArchive.metadataJson();
    if(settings.at("analysis_method").string()=="modal-matrix")
        return runModalObservabilityMatrix(sourcePath,source,settings,root,label);
    if(settings.at("analysis_method").string()=="eigenfunction-observability")
        return runEigenfunctionObservability(sourcePath,source,settings,root,label);
    if(settings.at("analysis_method").string()=="modal-gramian")  // LinearAnalysisModal.cpp
        return runModalGramian(sourcePath,source,settings,root,label);
    if(settings.at("analysis_method").string()=="klqe-assessment")
        return runKlqeAssessment(sourcePath,source,settings,root,label);
    const auto o=options(settings);
    NpzRecord out;json::Value::Object details;
    details["schema"]="dynsyslab.linear-analysis/v1";
    details["analysis_method"]=observabilityMethodName(o.method);
    std::optional<KoopmanModel> model;
    const bool matrixSource=sourceArchive.contains("A");
    std::size_t n=0;std::vector<double> A;
    if(matrixSource) {
        const auto& a=sourceArchive.at("A");
        if(a.shape.size()!=2||a.shape[0]!=a.shape[1])throw std::invalid_argument("A must be square");
        n=a.shape[0];A=a.doubles();
    } else {
        std::vector<std::size_t> selected;
        for(const auto& v:settings.at("selected_modes").array()) {
            const auto x=v.number(-1);
            if(x<0||std::floor(x)!=x)throw std::invalid_argument("invalid selected mode");
            selected.push_back(static_cast<std::size_t>(x));
        }
        model=loadKoopmanModel(sourcePath,selected);n=model->splitDimension();
        details["source_krom"]=sourcePath.generic_u8string();
    }
    const auto measurement=settings.at("measurement").string();
    if((matrixSource&&sourceArchive.contains("C"))||(!matrixSource&&!measurement.empty())) {
        ObservabilityResult r;
        if(model)r=analyzeObservability(*model,measurement,o);
        else {
            const auto& c=sourceArchive.at("C");
            if(c.shape.size()!=2||c.shape[1]!=n)throw std::invalid_argument("C must be outputs x states");
            r=analyzeObservability(A,c.doubles(),n,c.shape[0],o);
        }
        spectralArrays(out,"observability",n,r.gramian,r.singularValues,r.gramianEigenvalues,r.rightDirections);
        if(!r.matrix.empty())out.arrays["observability_matrix"]=NpyArray::fromDoubles(
            {r.outputDimension*r.blockCount,n},r.matrix);
        details["observability"]=summary(r);
    }
    const auto input=settings.at("control_map").string();
    if((matrixSource&&sourceArchive.contains("B"))||(!matrixSource&&!input.empty())) {
        ControllabilityResult r;
        if(model)r=analyzeControllability(*model,input,o);
        else {
            const auto& b=sourceArchive.at("B");
            if(b.shape.size()!=2||b.shape[0]!=n)throw std::invalid_argument("B must be states x inputs");
            r=analyzeControllability(A,b.doubles(),n,b.shape[1],o);
        }
        spectralArrays(out,"controllability",n,r.gramian,r.singularValues,r.gramianEigenvalues,r.leftDirections);
        if(!r.matrix.empty())out.arrays["controllability_matrix"]=NpyArray::fromDoubles(
            {n,r.inputDimension*r.blockCount},r.matrix);
        details["controllability"]=summary(r);
    }
    if(out.arrays.empty())throw std::invalid_argument("analysis requires C/measurement or B/control_map");
    auto serialized=json::dump(json::Value(details));
    out.arrays["metadata_json"]=NpyArray::fromString(metadata(source,"linear-analysis",label,
        serialized.substr(1,serialized.size()-2)));
    return publishObject(root, source, dynsyslab::core::object_kind::report, "linear-analysis",std::move(out));
}
} // namespace dynsyslab::workflow::postprocess_detail
