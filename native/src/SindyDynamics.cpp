#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
void aligned(const TimeSeries& a,const TimeSeries& b,bool sameWidth=true) {
    if(a.sampleCount()!=b.sampleCount()||(sameWidth&&a.componentCount()!=b.componentCount()))
        throw std::invalid_argument("SINDy coordinates and rates must align");
    for(std::size_t k=0;k<a.sampleCount();++k)
        if(a.time(k)!=b.time(k))throw std::invalid_argument("SINDy coordinate and rate clocks differ");
}
void finite(double x) {if(!std::isfinite(x))throw std::invalid_argument("SINDy requires finite data");}
}
SindyDynamicsModel::SindyDynamicsModel(Observable library,std::vector<double> coefficients,std::size_t outputs)
    :library_(std::move(library)),coefficients_(std::move(coefficients)),
     outputs_(outputs?outputs:library_.stateDimension()) {
    if(coefficients_.size()!=outputs_*library_.outputDimension())
        throw std::invalid_argument("SINDy coefficient shape mismatch");
    for(double v:coefficients_)finite(v);
}
std::vector<double> SindyDynamicsModel::evaluate(const double* q) const {
    const auto theta=library_.evaluateReal(q);
    std::vector<double> result(outputs_,0);
    for(std::size_t f=0;f<theta.size();++f)for(std::size_t j=0;j<outputs_;++j)
        result[j]+=theta[f]*coefficients_[f*outputs_+j];
    return result;
}
TimeSeries SindyDynamicsModel::evaluate(const TimeSeries& q) const {
    if(q.componentCount()!=dimension())throw std::invalid_argument("SINDy evaluation dimension mismatch");
    std::vector<double> result(q.sampleCount()*outputs_);
    for(std::size_t k=0;k<q.sampleCount();++k) {
        const auto rate=evaluate(q.valuesData()+k*dimension());
        std::copy(rate.begin(),rate.end(),result.begin()+k*outputs_);
    }
    return TimeSeries({q.timeData(),q.timeData()+q.sampleCount()},std::move(result),outputs_);
}
Observable SindyDynamicsModel::vectorField() const {
    const auto library=library_.expressions();
    if(library.empty())throw std::logic_error("Library has no portable expressions");
    if(!isVectorField())throw std::logic_error("A cross-coordinate rate law is not a closed vector field");
    std::vector<AnalyticExpression> rhs;
    for(std::size_t j=0;j<dimension();++j) {
        auto value=AnalyticExpression::constant(0,dimension());
        for(std::size_t f=0;f<library_.outputDimension();++f)
            if(coefficients_[f*dimension()+j]!=0)
                value=value+library[f]*coefficients_[f*dimension()+j];
        rhs.push_back(value);
    }
    return Observable::fromExpressions({"sindy_dynamics","SINDy vector field","","model",{},false},rhs);
}
json::Value SindyDynamicsModel::definition() const {
    json::Value::Array c;for(auto v:coefficients_)c.emplace_back(v);
    json::Value::Object out{{"schema","dynsyslab.sindy-dynamics/v1"},
        {"library",library_.expressionDefinition()},{"coefficients",c}};
    if(!isVectorField())out["outputs"]=double(outputs_);
    return out;
}
SindyDynamicsModel SindyDynamicsModel::fromDefinition(const json::Value& v) {
    if(v.at("schema").string()!="dynsyslab.sindy-dynamics/v1")
        throw std::invalid_argument("Unknown SINDy dynamics model schema");
    std::vector<double> c;for(const auto& n:v.at("coefficients").array()) {
        if(!std::holds_alternative<double>(n.storage()))throw std::invalid_argument("SINDy coefficients must be numbers");
        c.push_back(n.number());
    }
    return {Observable::fromExpressionDefinition(v.at("library")),c,std::size_t(v.at("outputs").number(0))};
}
Observable polynomialDynamicsLibrary(std::size_t d,std::size_t degree,
    const std::vector<double>& center,const std::vector<double>& scale) {
    if(!d||degree>32)throw std::invalid_argument("Invalid polynomial dictionary size");
    long double count=1;
    for(std::size_t k=1;k<=degree;++k) {
        count*=(static_cast<long double>(d)+k)/k;
        if(count>4096)throw std::invalid_argument("Polynomial dictionary exceeds 4096 features");
    }
    PolynomialDictionary dictionary{d,polynomialExponentTable(d,degree,true),
        center.empty()?std::vector<double>(d,0):center,scale.empty()?std::vector<double>(d,1):scale};
    return dictionary.asObservable();
}
std::vector<double> trapezoidWeights(const TimeSeries& series) {
    const auto n=series.sampleCount();std::vector<double> w(n,0);
    for(std::size_t k=1;k<n;++k){const double h=series.time(k)-series.time(k-1);w[k-1]+=h/2;w[k]+=h/2;}
    return w;
}
namespace {
// Regression rows: theta (rows x p) and targets (rows x d), each row with its weight.
struct Rows { std::vector<double> theta,targets,weights;std::size_t count{}; };
Rows instantaneousRows(const TimeSeries& q,const TimeSeries& rates,const Observable& library,const SindyDynamicsOptions& o) {
    const auto n=q.sampleCount(),p=library.outputDimension(),inputs=q.componentCount(),d=rates.componentCount();
    if(!o.sampleWeights.empty()&&o.sampleWeights.size()!=n)throw std::invalid_argument("SINDy sample weights need one value per sample");
    Rows r;r.count=n;r.theta.resize(n*p);r.targets.resize(n*d);r.weights.assign(n,1.0);
    for(std::size_t k=0;k<n;++k) {
        for(std::size_t j=0;j<inputs;++j)finite(q.value(k,j));
        const auto row=library.evaluateReal(q.valuesData()+k*inputs);
        for(std::size_t f=0;f<p;++f){finite(row[f]);r.theta[k*p+f]=row[f];}
        for(std::size_t j=0;j<d;++j){finite(rates.value(k,j));r.targets[k*d+j]=rates.value(k,j);}
        if(!o.sampleWeights.empty()) {
            if(!(o.sampleWeights[k]>=0)||!std::isfinite(o.sampleWeights[k]))throw std::invalid_argument("SINDy sample weights must be nonnegative");
            r.weights[k]=o.sampleWeights[k];
        }
    }
    return r;
}
// Window averages of Theta over [t, t + D] and the mean rate (q(t + D) - q(t)) / D, exact for the
// piecewise-linear interpolant of both.
Rows weakRows(const TimeSeries& q,const Observable& library,const SindyDynamicsOptions& o) {
    const auto n=q.sampleCount(),p=library.outputDimension(),d=q.componentCount();
    if(!o.weakStride)throw std::invalid_argument("SINDy weak_stride must be positive");
    std::vector<double> theta(n*p);
    for(std::size_t k=0;k<n;++k){const auto row=library.evaluateReal(q.valuesData()+k*d);
        for(std::size_t f=0;f<p;++f){finite(row[f]);theta[k*p+f]=row[f];}}
    Rows r;
    for(const double D:o.weakDurations) {
        if(!(D>0)||!std::isfinite(D))throw std::invalid_argument("SINDy weak durations must be positive");
        for(std::size_t i=0;i<n;i+=o.weakStride) {
            const double end=q.time(i)+D;
            if(end>q.time(n-1)*(1+1e-12)+1e-15)break;  // the window would cross the record's end
            std::vector<double> average(p,0),last(d);
            std::size_t k=i;
            for(;k+1<n&&q.time(k+1)<=end;++k){const double h=q.time(k+1)-q.time(k);
                for(std::size_t f=0;f<p;++f)average[f]+=h*(theta[k*p+f]+theta[(k+1)*p+f])/2;}
            double fraction=0;
            if(k+1<n&&end>q.time(k)){  // the partial interval up to the window end
                const double h=q.time(k+1)-q.time(k),part=end-q.time(k);fraction=part/h;
                for(std::size_t f=0;f<p;++f){const double atEnd=theta[k*p+f]+fraction*(theta[(k+1)*p+f]-theta[k*p+f]);
                    average[f]+=part*(theta[k*p+f]+atEnd)/2;}
            }
            for(std::size_t j=0;j<d;++j)last[j]=q.value(k,j)+(fraction>0?fraction*(q.value(k+1,j)-q.value(k,j)):0.0);
            for(std::size_t f=0;f<p;++f)r.theta.push_back(average[f]/D);
            for(std::size_t j=0;j<d;++j)r.targets.push_back((last[j]-q.value(i,j))/D);
            r.weights.push_back(1.0);++r.count;
        }
    }
    if(!r.count)throw std::invalid_argument("SINDy weak windows are longer than the record");
    return r;
}
}  // namespace
SindyDynamicsModel fitSindyDynamics(const TimeSeries& q,const TimeSeries& rates,
    const Observable& library,const SindyDynamicsOptions& options,SindyDynamicsFitDiagnostics* diagnostics) {
    aligned(q,rates,false);const auto p=library.outputDimension();
    const auto inputs=q.componentCount(),d=rates.componentCount();  // d: fitted outputs
    if(library.stateDimension()!=inputs||!d||options.threshold<0||!std::isfinite(options.threshold)
        ||options.singularValueTolerance<0||!std::isfinite(options.singularValueTolerance)||!options.maximumIterations
        ||!(options.ridge>=0)||!std::isfinite(options.ridge))
        throw std::invalid_argument("Invalid SINDy dimensions, sample count or options");
    const bool weak=!options.weakDurations.empty();
    if(weak&&d!=inputs)throw std::invalid_argument("weak-form SINDy fits a vector field (the predictors' own rates)");
    const auto rows=weak?weakRows(q,library,options):instantaneousRows(q,rates,library,options);
    const auto n=rows.count;
    if(n<p)throw std::invalid_argument("Invalid SINDy dimensions, sample count or options");
    // RMS scales of features and responses over the weighted rows.
    std::vector<double> scales(p,0),outputs(d,0),coefficients(p*d,0),root(n);
    double total=0;for(const double w:rows.weights)total+=w;
    if(!(total>0))throw std::invalid_argument("SINDy sample weights are all zero");
    for(std::size_t k=0;k<n;++k) {
        root[k]=std::sqrt(rows.weights[k]);
        for(std::size_t f=0;f<p;++f)scales[f]+=rows.weights[k]*rows.theta[k*p+f]*rows.theta[k*p+f]/total;
        for(std::size_t j=0;j<d;++j)outputs[j]+=rows.weights[k]*rows.targets[k*d+j]*rows.targets[k*d+j]/total;
    }
    for(auto& s:scales)s=std::sqrt(s);
    for(auto& s:outputs)s=std::sqrt(s);
    numerics::LeastSquaresOptions ls;ls.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    ls.singularValueTolerance=options.singularValueTolerance;
    if(diagnostics){diagnostics->relativeResidual.assign(d,{});diagnostics->supportSize.assign(d,{});}
    const double ridge=std::sqrt(options.ridge*total);  // the weighted row sum carries the sample count
    for(std::size_t j=0;j<d;++j) {
        if(outputs[j]==0)continue;
        std::vector<std::size_t> support;
        for(std::size_t f=0;f<p;++f)if(scales[f]>0)support.push_back(f);
        std::vector<double> solution;
        for(std::size_t iteration=0;iteration<=options.maximumIterations&&!support.empty();++iteration) {
            const auto m=support.size();const std::size_t extra=options.ridge>0?m:0;
            std::vector<double> matrix((n+extra)*m,0.0),rhs(n+extra,0.0);
            for(std::size_t k=0;k<n;++k) {
                rhs[k]=root[k]*rows.targets[k*d+j]/outputs[j];
                for(std::size_t f=0;f<m;++f)matrix[k*m+f]=root[k]*rows.theta[k*p+support[f]]/scales[support[f]];
            }
            for(std::size_t f=0;f<extra;++f)matrix[(n+f)*m+f]=ridge;  // ridge as augmented rows
            solution=numerics::solveLeastSquares(matrix,n+extra,m,rhs,1,ls);
            if(diagnostics) {
                double error=0;
                for(std::size_t k=0;k<n;++k) {
                    double predicted=0;for(std::size_t f=0;f<m;++f)predicted+=matrix[k*m+f]*solution[f];
                    error+=std::pow(predicted-rhs[k],2)/total;
                }
                diagnostics->relativeResidual[j].push_back(std::sqrt(error));
                diagnostics->supportSize[j].push_back(m);
            }
            std::vector<std::size_t> retained;
            for(std::size_t f=0;f<m;++f)if(std::abs(solution[f])>=options.threshold)retained.push_back(support[f]);
            if(retained==support||iteration==options.maximumIterations)break;
            support=std::move(retained);
        }
        for(std::size_t f=0;f<support.size();++f)
            coefficients[support[f]*d+j]=solution[f]*outputs[j]/scales[support[f]];
    }
    return {library,std::move(coefficients),d};
}
json::Value SindyDynamicsFitDiagnostics::toJson() const {
    json::Value::Array outputs;
    for(std::size_t j=0;j<relativeResidual.size();++j) {
        json::Value::Array steps;
        for(std::size_t k=0;k<relativeResidual[j].size();++k)
            steps.emplace_back(json::Value::Object{{"relative_residual",relativeResidual[j][k]},
                {"active_columns",double(supportSize[j][k])}});
        outputs.emplace_back(steps);
    }
    return outputs;
}
DynamicsError dynamicsError(const TimeSeries& expected,const TimeSeries& predicted) {
    aligned(expected,predicted);const auto d=expected.componentCount(),n=expected.sampleCount();
    DynamicsError e{std::vector<double>(d),std::vector<double>(d),std::vector<double>(d),std::vector<double>(d)};
    for(std::size_t k=0;k<n;++k)for(std::size_t j=0;j<d;++j) {
        const double truth=expected.value(k,j),delta=predicted.value(k,j)-truth;
        finite(truth);finite(delta);
        e.referenceRms[j]+=truth*truth/n;e.rmse[j]+=delta*delta/n;
        e.maximumAbsolute[j]=std::max(e.maximumAbsolute[j],std::abs(delta));
    }
    for(std::size_t j=0;j<d;++j) {
        e.referenceRms[j]=std::sqrt(e.referenceRms[j]);e.rmse[j]=std::sqrt(e.rmse[j]);
        e.relativeRms[j]=e.referenceRms[j]>0?e.rmse[j]/e.referenceRms[j]:0;
    }
    return e;
}
json::Value DynamicsError::toJson() const {
    json::Value::Object result;auto array=[](const std::vector<double>& x){
        json::Value::Array a;for(double v:x)a.emplace_back(v);return json::Value(a);};
    result["reference_rms"]=array(referenceRms);result["rmse"]=array(rmse);
    auto rel=array(relativeRms);for(std::size_t j=0;j<referenceRms.size();++j)
        if(referenceRms[j]==0)rel.array()[j]=json::Value{};
    result["relative_rms"]=rel;result["maximum_absolute"]=array(maximumAbsolute);return result;
}
}
