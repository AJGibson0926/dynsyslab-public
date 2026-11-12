#include "echo_multimode.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace echo_multi {
namespace {
std::vector<Complex> evaluate(const core::KoopmanModel& model,const Coordinates& x) {
    std::vector<Complex> out;
    for(std::size_t k=0;k<x.time.size();++k) {
        const auto z=model.evaluateCoordinates(x.x.data()+2*k);out.insert(out.end(),z.begin(),z.end());
    }return out;
}
std::vector<core::HistoryInputSample> pairs(const Coordinates& x,const std::vector<Complex>& z,
    std::size_t mode,std::size_t modes) {
    std::vector<core::HistoryInputSample> out;
    for(std::size_t k=0;k+1<x.time.size();++k)
        out.push_back({{x.x[2*k],x.x[2*k+1]},z[k*modes+mode],z[(k+1)*modes+mode],x.input[k]});
    return out;
}
Json inverse(const core::KoopmanModel& model,const Coordinates& x,const std::vector<double>& scales,
             const fs::path& path) {
    double ue=0,uu=0,perp=0,qq=0,modelerr=0;
    std::vector<double> time,known,recovered;
    for(std::size_t k=0;k+1<x.time.size();++k)if(x.time[k]>=30) {
        const auto z=model.evaluateCoordinates(x.x.data()+2*k),next=model.evaluateCoordinates(x.x.data()+2*k+2);
        auto b=model.controlMapCoordinates(x.x.data()+2*k,"fitted");
        std::vector<Complex> q(z.size());double dot=0,bb=0;
        for(std::size_t j=0;j<z.size();++j) {
            const auto l=model.continuousEigenvalues()[j],mu=std::exp(l*.01);
            b[j]*=(mu-1.)/l; q[j]=(next[j]-mu*z[j])/scales[j];b[j]/=scales[j];
            dot+=(std::conj(b[j])*q[j]).real();bb+=std::norm(b[j]);
        }
        const double u=dot/std::max(bb,1e-300);
        for(std::size_t j=0;j<z.size();++j) {
            perp+=std::norm(q[j]-b[j]*u);qq+=std::norm(q[j]);
            modelerr+=std::norm(q[j]-b[j]*x.input[k]);
        }
        ue+=(u-x.input[k])*(u-x.input[k]);uu+=x.input[k]*x.input[k];
        time.push_back(x.time[k]);known.push_back(x.input[k]);recovered.push_back(u);
    }
    core::NpzRecord a;
    a.arrays["time"]=core::NpyArray::fromDoubles({time.size()},time);
    a.arrays["known_input"]=core::NpyArray::fromDoubles({time.size()},known);
    a.arrays["inverse_input"]=core::NpyArray::fromDoubles({time.size()},recovered);
    core::writeNpzAtomicNew(path,a);
    return Json::Object{{"driver_nrmse",std::sqrt(ue/std::max(uu,1e-300))},
        {"unachievable_residual_fraction",std::sqrt(perp/std::max(qq,1e-300))},
        {"known_driver_model_residual",std::sqrt(modelerr/std::max(qq,1e-300))}};
}
Json trial(const fs::path& out,const std::string& id,const Frozen& f,
    const core::KoopmanModel& model,const core::SystemConfiguration& config,
    const Measurements& target,const Coordinates& reference,const std::vector<double>& scales,double penalty) {
    const auto n=model.splitDimension(),m=model.liftedDimension();
    if(n!=2*m)throw std::runtime_error("this experiment selects complex modes only");
    // Equivalent w_j=z_j/targetRMS_j scaling avoids an ill-conditioned CARE.
    // A is unchanged: both real components in each eigenmode share a scale.
    std::vector<double> Q(n*n);for(std::size_t j=0;j<n;++j)Q[j*n+j]=1./m;
    const auto A=model.generatorSplit();
    core::CausalHistoryProjector projector(f.basis);
    std::size_t tick=0;double held=0,maxCare=0;
    std::vector<double> traceTime,modeErrors;
    EchoController feedback=[&](double t,double e) {
        projector.append(e);const auto now=tick++;
        if(now%f.stride || t<6.-1e-8)return held;
        const auto xi=projector.projectLatest();if(!xi)return 0.;
        const auto index=static_cast<std::size_t>(std::llround((t-reference.time.front())/.01));
        if(index+1>=reference.time.size())return held;
        const auto z=model.evaluateCoordinatesSplit(xi->data());
        const auto ref=model.evaluateCoordinatesSplit(reference.x.data()+2*index);
        const auto bc=model.controlMapCoordinates(xi->data(),"fitted");
        std::vector<double> B;for(std::size_t j=0;j<m;++j){B.push_back(bc[j].real()/scales[j]);B.push_back(bc[j].imag()/scales[j]);}
        const auto lqr=core::solveContinuousLQR(A,B,Q,{penalty/.01},n,1);
        maxCare=std::max(maxCare,lqr.relativeResidual);
        double value=0;for(std::size_t j=0;j<n;++j)value-=lqr.gain[j]*(z[j]-ref[j])/scales[j/2];
        if(!std::isfinite(value))throw std::runtime_error("nonfinite KLQR input");
        held=std::clamp(value,-.6,.6);traceTime.push_back(t);
        for(std::size_t j=0;j<m;++j)modeErrors.push_back(std::hypot(z[2*j]-ref[2*j],z[2*j+1]-ref[2*j+1])/scales[j]);
        return held;
    };
    // Unchanged acquisition adapter alone sees plant state. No target input enters feedback.
    simulate(config,out,id,f,target.time.back(),1.,feedback,Json::Object{
        {"role","multimode closed-loop test"},{"modes",double(m)},{"latent_coordinates",2.},
        {"R_multiplier",penalty},{"Q","w=psi/targetRMS; Q=I/m; algebraically equivalent to raw Q_j=1/(m RMS_j^2)"},
        {"law","existing native continuous CARE; u=-K(z-z_target); sample/hold 0.01; no feedforward"}});
    core::NpzRecord trace;
    trace.arrays["time"]=core::NpyArray::fromDoubles({traceTime.size()},traceTime);
    trace.arrays["normalized_mode_error"]=core::NpyArray::fromDoubles({traceTime.size(),m},modeErrors);
    core::writeNpzAtomicNew(out/(id+".feedback.npz"),trace);
    const auto actual=readMeasurements(out/(id+".measurements.npz"));
    return Json::Object{{"id",id},{"late",metrics(actual,target,30)},
        {"whole_driven",metrics(actual,target,6)},{"final_10tau",metrics(actual,target,50)},
        {"maximum_care_relative_residual",maxCare},{"R_multiplier",penalty},{"modes",double(m)}};
}
}
Json runModelTrial(const fs::path& out,const std::string& id,const Frozen& f,
    const core::KoopmanModel& model,const core::SystemConfiguration& config,
    const Measurements& target,const Coordinates& reference,const std::vector<double>& scales,double penalty) {
    return trial(out,id,f,model,config,target,reference,scales,penalty);
}
void runControls(const fs::path& root,const fs::path& out,const std::string& window,
    const Frozen& f,const core::KoopmanModel& bank,const core::SystemConfiguration& config,bool joint) {
    const auto inputs=root/"archives/experiments/echo-input-regression-20261006-v1";
    Json::Array fits,jointFits;
    for(const std::string targetName:{"mild","nonlinear"}) {
        const auto target=readMeasurements(inputs/("target-"+targetName+".measurements.npz"));
        const auto reference=project(target,f);const auto targetZ=evaluate(bank,reference);
        // Check online/offline projection parity once per target, using the same samples.
        core::CausalHistoryProjector checkProjector(f.basis);
        double projectionError=0,projectionEnergy=0;
        for(std::size_t k=0;k<target.echo.size();++k) {
            checkProjector.append(target.echo[k]);
            if(k%2000 || target.time[k]<reference.time.front())continue;
            const auto xi=checkProjector.projectLatest();
            const auto row=std::size_t(std::llround((target.time[k]-reference.time.front())/.01));
            for(std::size_t j=0;j<2;++j) {
                projectionError+=std::pow((*xi)[j]-reference.x[2*row+j],2);
                projectionEnergy+=std::pow(reference.x[2*row+j],2);
            }
        }
        if(std::sqrt(projectionError/std::max(projectionEnergy,1e-300))>1e-10)
            throw std::runtime_error("online/offline projection parity failed");
        std::vector<std::vector<core::HistoryInputSample>> rows(3),validation(3);
        for(int k=1;k<=8;++k) {
            const auto rec=readMeasurements(inputs/("local-"+targetName+"-"+std::to_string(k)+".measurements.npz"));
            const auto x=project(rec,f,12,rec.time.back()-6.02);
            const auto z=evaluate(bank,x);
            for(std::size_t j=0;j<3;++j) {
                const auto p=pairs(x,z,j,3);auto& dst=k<=6?rows[j]:validation[j];
                dst.insert(dst.end(),p.begin(),p.end());
            }
        }
        std::vector<core::HistoryInputMap> maps;std::vector<double> rms(3);
        ModeSamples targetRows(3);
        for(std::size_t j=0;j<3;++j) {
            targetRows[j]=pairs(reference,targetZ,j,3);
            const auto map=core::fitHistoryInputMap(rows[j],std::exp(bank.continuousEigenvalues()[j]*.01),.01,1,{1,1});
            maps.push_back(map);
            const auto train=core::scoreHistoryInputMap(map,rows[j]),val=core::scoreHistoryInputMap(map,validation[j]);
            const auto test=core::scoreHistoryInputMap(map,pairs(reference,targetZ,j,3));
            fits.emplace_back(Json::Object{{"target",targetName},{"mode",double(j+1)},
                {"train_residual",train.relativeResidual},{"heldout_waveform_residual",val.relativeResidual},
                {"target_residual",test.relativeResidual},{"condition",map.conditionNumber}});
            // Same target-only amplitude normalization for a mode in every subset.
            for(std::size_t k=0;k<reference.time.size();++k)rms[j]+=std::norm(targetZ[k*3+j]);
            rms[j]=std::sqrt(rms[j]/reference.time.size());
            if(!(rms[j]>1e-15))throw std::runtime_error("vanishing target eigenfunction");
        }
        // All seven nonempty subsets: fixed rows for historical reproduction;
        // --joint refits every output from the predictors in this precise subset.
        for(int mask=1;mask<8;++mask) {
            std::vector<std::size_t> selected;std::string tag;
            for(std::size_t j=0;j<3;++j)if(mask&(1<<j)){selected.push_back(j);tag+=std::to_string(j+1);}
            auto model=bank.extract(selected);const auto m=selected.size();
            core::BilinearControlMap b;b.modes=m;b.inputs=1;b.operators.resize(m*m);b.conjugateOperators.resize(m*m);
            b.provenance="isolated nested experiment; existing per-mode bilinear fit; fixed diagonal rows across subsets";
            std::vector<double> scales;
            for(std::size_t j=0;j<m;++j) {
                const auto one=std::get<core::BilinearControlMap>(maps[selected[j]].canonicalRepresentation());
                b.offsets.push_back(one.offsets[0]);b.operators[j*m+j]=one.operators[0];
                b.conjugateOperators[j*m+j]=one.conjugateOperators[0];scales.push_back(rms[selected[j]]);
            }
            Json jointReport;
            std::map<std::string,core::NpyArray> extra;
            if(joint) {
                const auto fit=fitJointInput(rows,selected,bank.continuousEigenvalues());
                b=fit.map;
                jointReport=Json::Object{{"target",targetName},{"subset",tag},{"modes",double(m)},
                    {"columns",double(fit.columns)},{"rank",double(fit.rank)},{"condition",fit.condition},
                    {"train",scoreJointInput(fit,rows,selected,bank.continuousEigenvalues())},
                    {"heldout_waveforms",scoreJointInput(fit,validation,selected,bank.continuousEigenvalues())},
                    {"target_score_only",scoreJointInput(fit,targetRows,selected,bank.continuousEigenvalues())}};
                jointFits.push_back(jointReport);
                extra["sampled_joint_coefficients"]=core::NpyArray::fromComplex({m,fit.columns},fit.coefficients);
            }
            model.setControlMap("fitted",b);
            // Canonical multi-mode map must equal the same independently fitted rows.
            for(std::size_t k=0;!joint && k<reference.time.size();k+=100) {
                const auto z=model.evaluateCoordinates(reference.x.data()+2*k);
                const auto bs=model.controlMapCoordinates(reference.x.data()+2*k,"fitted");
                for(std::size_t j=0;j<m;++j) {
                    const auto l=model.continuousEigenvalues()[j];
                    const auto expected=maps[selected[j]].evaluate(z[j])/((std::exp(l*.01)-1.)/l);
                    if(std::abs(bs[j]-expected)>1e-10*std::max(1.,std::abs(expected)))
                        throw std::runtime_error("canonical map assembly parity failed");
                }
            }
            const auto prefix=window+"-"+targetName+"-m"+tag;
            if(joint && !fs::exists(out/(prefix+".fit.json")))saveJson(out/(prefix+".fit.json"),jointReport);
            Json::Array savedScales;for(double s:scales)savedScales.emplace_back(s);
            if(!fs::exists(out/(prefix+".krom.npz")))core::saveKoopmanModel(out/(prefix+".krom.npz"),model,
                core::json::dump(Json::Object{{"schema","dynsyslab.echo-multimode-control/v1"},
                    {"basis_reference",joint?(root/"archives/experiments/echo-multimode-20261009-v2"/(window+"-bank.krom.npz")).generic_string():window+"-bank.krom.npz"},
                    {"B_policy",joint?"joint refit for this subset; all rows and cross-mode predictors":"fixed independently fitted rows"},
                    {"calibration","local records 1-6; records 7-8 withheld"},{"target_rms_scales",savedScales}}),extra);
            const auto invPath=out/(prefix+".inverse.json");
            if(!fs::exists(invPath))saveJson(invPath,inverse(model,reference,scales,out/(prefix+".inverse.npz")));
            // Additional m1 cases isolate the 1/m reduction of fundamental weight:
            // m12 at R=1 vs m1 at R=2; m123 at R=1 vs m1 at R=3 (also R=100).
            const auto penalties=mask==1?std::vector<double>{100.,1.,.01,2.,3.,200.,300.}
                                       :std::vector<double>{100.,1.,.01};
            for(double penalty:penalties) {
                const auto id=prefix+"-r"+std::to_string(penalty);const auto resultPath=out/(id+".json");
                if(fs::exists(resultPath))continue;
                std::cout<<"RUN "<<id<<std::endl;
                try {const auto result=trial(out,id,f,model,config,target,reference,scales,penalty);
                    saveJson(resultPath,result);std::cout<<core::json::dump(result)<<std::endl;
                }catch(const std::exception& e){saveJson(resultPath,Json::Object{{"id",id},{"failure",e.what()},
                    {"R_multiplier",penalty},{"modes",double(m)}});std::cout<<"FAILED "<<e.what()<<std::endl;}
            }
        }
    }
    if(!fs::exists(out/(window+"-maps.json")))saveJson(out/(window+"-maps.json"),fits);
    if(joint && !fs::exists(out/(window+"-joint-maps.json")))saveJson(out/(window+"-joint-maps.json"),jointFits);
}
}
