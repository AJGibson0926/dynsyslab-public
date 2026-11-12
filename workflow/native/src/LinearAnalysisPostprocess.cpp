#include "PostprocessInternal.hpp"
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
std::filesystem::path runLinearAnalysis(const std::filesystem::path& sourcePath,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    using namespace dynsyslab::core;
    const auto sourceArchive=readNpz(sourcePath);
    SourceData source;source.path=sourcePath;source.metadataJson=sourceArchive.metadataJson();
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
    return publish(root,source,"linear-analysis",label,std::move(out));
}
} // namespace dynsyslab::workflow::postprocess_detail
