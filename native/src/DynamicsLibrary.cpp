#include <dynsyslab/core/DynamicsLibrary.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
DynamicsLibrary buildDynamicsLibrary(const TimeSeries& q,const DynamicsLibraryOptions& o) {
    const auto d=q.componentCount(),n=q.sampleCount();
    if(o.inverseDegree>16)throw std::invalid_argument("Inverse library degree exceeds 16");
    std::vector<double> center(d,0),scale(d,1),lo(d),hi(d);
    for(std::size_t j=0;j<d;++j) {
        lo[j]=hi[j]=q.value(0,j);
        for(std::size_t k=0;k<n;++k) {
            const double v=q.value(k,j);
            if(!std::isfinite(v))throw std::invalid_argument("Library training coordinates must be finite");
            lo[j]=std::min(lo[j],v);hi[j]=std::max(hi[j],v);
            if(o.normalize)center[j]+=v/n;
        }
        if(o.normalize) {
            double variance=0;for(std::size_t k=0;k<n;++k)variance+=std::pow(q.value(k,j)-center[j],2)/n;
            scale[j]=std::sqrt(variance);if(scale[j]<=1e-15)scale[j]=1;
        }
    }
    const auto polynomial=polynomialDynamicsLibrary(d,o.polynomialDegree,center,scale);
    auto terms=polynomial.expressions();auto names=polynomial.descriptor().outputNames;
    std::vector<std::string> excluded;
    std::vector<AnalyticExpression> raw,z;
    for(std::size_t j=0;j<d;++j) {
        raw.push_back(AnalyticExpression::state(d,j));z.push_back((raw.back()-center[j])/scale[j]);
    }
    const auto add=[&](const AnalyticExpression& term,const std::string& name,std::size_t axis) {
        terms.push_back(term);names.push_back(name);
        if(o.mixedTerms)for(std::size_t other=0;other<d;++other)if(other!=axis)
            for(std::size_t degree=1;degree<=2;++degree) {
                // Raw factors retain familiar rational terms such as v/r^2.
                terms.push_back(term*raw[other].power(double(degree)));
                names.push_back(name+" * q"+std::to_string(other+1)+"^"+std::to_string(degree));
            }
        if(terms.size()>4096)throw std::invalid_argument("Dynamics candidate library exceeds 4096 terms");
    };
    for(std::size_t j=0;j<d;++j) {
        const auto label="z"+std::to_string(j+1),original="q"+std::to_string(j+1);
        if(o.trigonometric)for(double w:o.frequencies) {
            if(!std::isfinite(w)||w<=0)throw std::invalid_argument("Library frequencies must be finite and positive");
            add((z[j]*w).unary("sin"),"sin("+std::to_string(w)+"*"+label+")",j);
            add((z[j]*w).unary("cos"),"cos("+std::to_string(w)+"*"+label+")",j);
        }
        if(o.reciprocal) {
            if(lo[j]>0||hi[j]<0)for(std::size_t k=1;k<=o.inverseDegree;++k)
                add(raw[j].power(-double(k)),original+"^-"+std::to_string(k),j);
            else excluded.push_back(original+" negative powers: training interval includes zero");
            add((1+z[j].power(2)).power(-1),"1/(1+"+label+"^2)",j);
        }
        if(o.logarithmic) {
            if(lo[j]>0)add(raw[j].unary("log"),"log("+original+")",j);
            else excluded.push_back("log("+original+"): training interval is not strictly positive");
            add((1+z[j].power(2)).unary("log"),"log(1+"+label+"^2)",j);
        }
        if(o.exponential)for(double rate:o.exponentialRates) {
            if(!std::isfinite(rate)||rate==0)throw std::invalid_argument("Library exponential rates must be finite and nonzero");
            const double argument=std::max(rate*(lo[j]-center[j])/scale[j],rate*(hi[j]-center[j])/scale[j]);
            const auto name="exp("+std::to_string(rate)+"*"+label+")";
            if(argument>40){excluded.push_back(name+": training argument exceeds 40");continue;}
            add((z[j]*rate).unary("exp"),name,j);
        }
        if(o.hyperbolic) {
            add(z[j].unary("tanh"),"tanh("+label+")",j);
            add((1+z[j].power(2)).power(.5),"sqrt(1+"+label+"^2)",j);
        }
    }
    return {Observable::fromExpressions({"dynamics_library","Mixed dynamics candidate library","","model",names,false},terms),excluded};
}
}
