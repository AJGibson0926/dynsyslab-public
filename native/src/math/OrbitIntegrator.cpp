#include <dynsyslab/math/OrbitIntegrator.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
namespace dynsyslab::math {
namespace {
constexpr int kLevels=8;
using V=std::vector<double>;
// Gragg modified midpoint with the symmetric smoothing endpoint.
void midpoint(const FlowField& f,const V& x,double H,int n,V& out){
    const auto d=x.size();const double h=H/n;V a=x,b(d),k(d);
    f(a,k);for(std::size_t i=0;i<d;++i)b[i]=a[i]+h*k[i];
    for(int m=1;m<n;++m){f(b,k);for(std::size_t i=0;i<d;++i){const double c=a[i]+2*h*k[i];a[i]=b[i];b[i]=c;}}
    f(b,k);out.resize(d);for(std::size_t i=0;i<d;++i)out[i]=.5*(a[i]+b[i]+h*k[i]);
}
std::size_t stencil(const V& x,double at,std::size_t count){
    const auto up=std::size_t(std::upper_bound(x.begin(),x.end(),at)-x.begin());
    const auto left=up==0?0:up-1;const auto back=(count-1)/2;
    return std::min(left>back?left-back:0,x.size()-count);
}
double lagrange(const V& x,const V& f,std::size_t first,std::size_t count,double at){
    double sum=0;
    for(std::size_t i=first;i<first+count;++i){double w=1;
        for(std::size_t j=first;j<first+count;++j)if(j!=i)w*=(at-x[j])/(x[i]-x[j]);
        sum+=w*f[i];}
    return sum;
}
}
void requireIncreasing(const V& x,const char* name){
    if(x.size()<2)throw std::invalid_argument(std::string(name)+" needs at least two values");
    for(std::size_t i=0;i<x.size();++i){
        if(!std::isfinite(x[i]))throw std::invalid_argument(std::string(name)+" must be finite");
        if(i&&x[i]<=x[i-1])throw std::invalid_argument(std::string(name)+" must be strictly increasing");
    }
}
void extrapolatedStep(const FlowField& f,const V& x,double H,V& out,V* lower){
    std::vector<V> row(kLevels);V current;
    for(int k=0;k<kLevels;++k){
        midpoint(f,x,H,2*(k+1),current);
        V previous=row[0];row[0]=current;
        for(int j=1;j<=k;++j){
            const double ratio=double(k+1)/double(k+1-j),denominator=ratio*ratio-1;
            V next(x.size());
            for(std::size_t i=0;i<x.size();++i)next[i]=row[j-1][i]+(row[j-1][i]-previous[i])/denominator;
            previous=row[j];row[j]=std::move(next);
        }
    }
    out=row[kLevels-1];if(lower)*lower=row[kLevels-2];
}
double adaptiveFlowStep(const FlowField& f,V& x,double& h,const FlowOptions& o,double limit){
    if(!(o.relativeTolerance>0)||!(o.absoluteTolerance>0))throw std::invalid_argument("Flow tolerances must be positive");
    if(!(h>0)||!std::isfinite(h))h=o.initialStep>0?o.initialStep:1e-2;
    V high,low;
    for(int attempt=0;attempt<80;++attempt){
        double trial=h;if(o.maximumStep>0)trial=std::min(trial,o.maximumStep);
        const bool truncated=limit<trial;if(truncated)trial=limit;
        if(!(trial>0))throw std::invalid_argument("Flow step limit must be positive");
        extrapolatedStep(f,x,trial,high,&low);
        double error=0;bool finite=true;
        for(std::size_t i=0;i<x.size();++i){
            if(!std::isfinite(high[i])){finite=false;break;}
            const double scale=o.absoluteTolerance+o.relativeTolerance*std::max(std::abs(x[i]),std::abs(high[i]));
            error=std::max(error,std::abs(high[i]-low[i])/scale);
        }
        if(finite&&error<=1){
            const double growth=error>0?std::clamp(.9*std::pow(error,-1./(2*kLevels-1)),.2,3.):3.;
            h=truncated?std::max(h,trial*growth):trial*growth;x=high;return trial;
        }
        h=trial*(finite?std::clamp(.9*std::pow(error,-1./(2*kLevels-1)),.1,.9):.25);
        if(h<std::numeric_limits<double>::min()*1e3)break;
    }
    throw std::runtime_error("Adaptive flow step underflow");
}
void integrateFlow(const FlowField& f,V& x,double duration,const FlowOptions& o){
    if(!(duration>=0)||!std::isfinite(duration))throw std::invalid_argument("Flow duration must be finite and nonnegative");
    double done=0,h=o.initialStep;std::size_t steps=0;
    while(done<duration){
        const double left=duration-done;
        if(left<=1e-15*duration)break;
        done+=adaptiveFlowStep(f,x,h,o,left);
        if(++steps>o.maximumSteps)throw std::runtime_error("Flow exceeded maximum steps");
    }
}
double interpolateTable(const V& x,const V& f,double at){
    if(x.size()!=f.size()||x.empty())throw std::invalid_argument("Table sizes differ");
    if(x.size()==1)return f[0];
    const auto count=std::min<std::size_t>(4,x.size());
    return lagrange(x,f,stencil(x,at,count),count,at);
}
V cumulativeTableIntegral(const V& x,const V& f){
    requireIncreasing(x,"Table abscissa");
    if(x.size()!=f.size())throw std::invalid_argument("Table sizes differ");
    V out(x.size());const double g=.5/std::sqrt(3.);
    const auto count=std::min<std::size_t>(4,x.size());
    for(std::size_t i=0;i+1<x.size();++i){
        const double mid=.5*(x[i]+x[i+1]),w=x[i+1]-x[i];
        const auto first=stencil(x,mid,count);
        out[i+1]=out[i]+.5*w*(lagrange(x,f,first,count,mid-g*w)+lagrange(x,f,first,count,mid+g*w));
    }
    return out;
}
V tableDerivative(const V& x,const V& f){
    requireIncreasing(x,"Table abscissa");
    if(x.size()!=f.size())throw std::invalid_argument("Table sizes differ");
    const auto n=x.size();V d(n);
    if(n==2){d[0]=d[1]=(f[1]-f[0])/(x[1]-x[0]);return d;}
    const auto three=[&](std::size_t a,double at){
        const double x0=x[a],x1=x[a+1],x2=x[a+2];
        return f[a]*(2*at-x1-x2)/((x0-x1)*(x0-x2))+f[a+1]*(2*at-x0-x2)/((x1-x0)*(x1-x2))
            +f[a+2]*(2*at-x0-x1)/((x2-x0)*(x2-x1));};
    for(std::size_t i=0;i<n;++i)d[i]=three(i==0?0:(i==n-1?n-3:i-1),x[i]);
    return d;
}
}
