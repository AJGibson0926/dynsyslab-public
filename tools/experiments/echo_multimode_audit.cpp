#include "echo_multimode.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace echo_multi {
namespace {
Json readJson(const fs::path& path){std::ifstream f(path);std::stringstream s;s<<f.rdbuf();return core::json::parse(s.str());}
double ratio(double a,double b){return std::sqrt(a/std::max(b,1e-300));}
Json selfTests() {
    const std::vector<double>A{-.1,-6.,6.,-.1},B{.2,.3},Q{1.,0.,0.,1.};
    const auto single=core::solveContinuousLQR(A,B,Q,{1.},2,1);
    double duplicationError=0,scalingError=0;
    for(std::size_t m=2;m<=3;++m) {
        const auto n=2*m;std::vector<double>a(n*n),b(n),q(n*n);
        for(std::size_t j=0;j<m;++j)for(std::size_t k=0;k<2;++k) {
            b[2*j+k]=B[k];q[(2*j+k)*n+2*j+k]=1./m;
            for(std::size_t l=0;l<2;++l)a[(2*j+k)*n+2*j+l]=A[k*2+l];
        }
        const auto fit=core::solveContinuousLQR(a,b,q,{1.},n,1);
        for(std::size_t k=0;k<2;++k){double gain=0;for(std::size_t j=0;j<m;++j)gain+=fit.gain[2*j+k];
            duplicationError=std::max(duplicationError,std::abs(gain-single.gain[k]));}
    }
    for(double s:{.01,100.}) {
        auto q=Q,b=B;for(auto&v:q)v/=s*s;for(auto&v:b)v/=s;
        const auto raw=core::solveContinuousLQR(A,B,q,{1.},2,1);
        const auto normalized=core::solveContinuousLQR(A,b,Q,{1.},2,1);
        for(std::size_t k=0;k<2;++k)scalingError=std::max(scalingError,
            std::abs(raw.gain[k]-normalized.gain[k]/s)/std::max(1.,std::abs(raw.gain[k])));
    }
    if(duplicationError>1e-9||scalingError>1e-9)throw std::runtime_error("multi-mode normalization tests failed");
    return Json::Object{{"duplicate_mode_weight_invariance_max_error",duplicationError},
        {"normalization_equivalence_max_relative_error",scalingError},
        {"online_offline_projection_and_canonical_map_parity","asserted during all four target/window runs"}};
}
Json truthScore(const fs::path& actual,const fs::path& target) {
    const auto a=core::readNpz(actual),b=core::readNpz(target);
    const auto t=a.at("time").doubles(),tt=b.at("time").doubles();
    const auto x=a.at("states").doubles(),y=b.at("states").doubles();
    if(t!=tt || x.size()!=y.size())throw std::runtime_error("truth time/shape mismatch");
    double re=0,rr=0,ve=0,vv=0,rmin=1e300,rmax=-1e300;
    for(std::size_t k=0;k<t.size();++k)if(t[k]>=30) {
        re+=std::pow(x[2*k]-y[2*k],2);rr+=std::pow(y[2*k]-1.,2);
        ve+=std::pow(x[2*k+1]-y[2*k+1],2);vv+=std::pow(y[2*k+1],2);
        rmin=std::min(rmin,x[2*k]);rmax=std::max(rmax,x[2*k]);
    }
    return Json::Object{{"radius_oscillation_nrmse",ratio(re,rr)},
        {"velocity_nrmse",ratio(ve,vv)},{"minimum_radius",rmin},{"maximum_radius",rmax}};
}
Json causal(const fs::path& path,const Frozen& f,const core::KoopmanModel& model) {
    const auto a=core::readNpz(path);
    if(a.contains("states"))throw std::runtime_error("causal diagnostic rejects state access");
    const auto time=a.at("time").doubles(),hist=a.at("echo_history").doubles();
    const auto tails=a.at("echo_tail").doubles(),input=a.at("input").doubles();
    const auto d=hist.size()/time.size(),stride=tails.size()/(time.size()*8),m=model.liftedDimension();
    std::vector<double> drift(m),act(m),berr(m),bnorm(m);
    for(std::size_t k=0;k<time.size();++k) {
        auto value=[&](int branch){core::CausalHistoryProjector p(f.basis);
            for(std::size_t j=0;j<d;++j)p.append(hist[k*d+j]);
            if(branch>=0)for(std::size_t j=0;j<stride;++j)p.append(tails[(k*8+branch)*stride+j]);
            const auto xi=p.projectLatest();return std::make_pair(*xi,model.evaluateCoordinates(xi->data()));};
        const auto base=value(-1),zero=value(0),plus=value(1),minus=value(2),forced=value(7);
        const auto bfit=model.controlMapCoordinates(base.first.data(),"fitted");
        const double delta=input[k*8+1];
        for(std::size_t j=0;j<m;++j) {
            const auto l=model.continuousEigenvalues()[j],mu=std::exp(l*.01);
            const auto q0=zero.second[j]-mu*base.second[j];
            const auto b=(plus.second[j]-minus.second[j])/(2*delta);
            drift[j]+=std::norm(q0);act[j]+=std::norm(forced.second[j]-zero.second[j]);
            berr[j]+=std::norm(bfit[j]*(mu-1.)/l-b);bnorm[j]+=std::norm(b);
        }
    }
    Json::Array rows;for(std::size_t j=0;j<m;++j)rows.emplace_back(Json::Object{{"mode",double(j+1)},
        {"zero_input_drift_over_actual_input_effect",ratio(drift[j],act[j])},
        {"regressed_B_relative_error_to_causal_B",ratio(berr[j],bnorm[j])}});
    return rows;
}
}
void auditResults(const fs::path& root,const fs::path& out) {
    const auto inputs=root/"archives/experiments/echo-input-regression-20261006-v1";
    Json::Array results,causalRows;int failed=0,completed=0;
    for(const auto& entry:fs::directory_iterator(out)) {
        const auto name=entry.path().filename().string();
        if(entry.path().extension()!=".json" || name.find("-r")==std::string::npos)continue;
        auto row=readJson(entry.path());const auto id=row.at("id").string();
        if(row.object().count("failure")){++failed;results.push_back(row);continue;}
        ++completed;const auto target=id.find("-mild-")!=std::string::npos?"mild":"nonlinear";
        row.object()["truth_late"]=truthScore(out/(id+".truth.npz"),inputs/(std::string("target-")+target+".truth.npz"));
        const auto recomputed=metrics(readMeasurements(out/(id+".measurements.npz")),
            readMeasurements(inputs/(std::string("target-")+target+".measurements.npz")),30);
        for(const auto key:{"echo_nrmse","driver_nrmse","saturation_fraction"})
            if(std::abs(recomputed.at(key).number()-row.at("late").at(key).number())>1e-12)
                throw std::runtime_error("saved control metric parity failed");
        results.push_back(row);
    }
    for(const std::string window:{"long","short"}) {
        const auto source=inputs/(window=="long"?"refit-xi2.krom.npz":"short-history-0p05/refit-xi2.krom.npz");
        const auto old=core::loadKoopmanModel(source);const auto f=frozen(core::readNpz(source),old.eigenfunctions());
        for(const std::string target:{"mild","nonlinear"}) {
            const auto path=out/(window+"-"+target+"-m123.krom.npz");
            auto model=core::loadKoopmanModel(path);
            // Custom map names are loaded through the existing independent-map API.
            model.setControlMap("fitted",core::loadApproximateControlMap(path,"control_fitted",3,1));
            causalRows.emplace_back(Json::Object{{"window",window},{"target",target},
                {"modes",causal(root/"archives/experiments/echo-control-diagnosis-20261007-v2"/(target+"-probes.npz"),f,model)}});
        }
    }
    const auto result=Json::Object{{"schema","dynsyslab.echo-multimode-audit/v1"},{"checks",selfTests()},
        {"completed",double(completed)},{"failed",double(failed)},{"results",results},{"causal",causalRows},
        {"causal_probes_role","previous target-history interventions, diagnostic only, never used for fitting"}};
    saveJson(out/"audit.json",result);std::cout<<"AUDIT completed="<<completed<<" failed="<<failed<<std::endl;
}
}
