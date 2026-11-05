// Measurement-only: no state archives, plant evaluators or equations of motion.
#include "echo_pulse.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace echo_pulse {
Json analyze(const fs::path& output,const Preparation& p,const std::string& window,
             const Frozen& f,const core::KoopmanModel& model) {
    const auto path=output/(window+"-"+p.id+".response.npz"),json=output/(window+"-"+p.id+".response.json");
    if(fs::exists(json))return readJson(json);
    const auto a=core::readNpz(output/(p.id+".pulse.npz"));
    if(a.contains("states")||a.contains("values"))throw std::runtime_error("state-bearing pulse rejected");
    const auto hist=a.at("echo_history").doubles(),echo=a.at("echo_continuations").doubles();
    const auto time=a.at("time_since_pulse").doubles(),amp=a.at("input_amplitudes").doubles();
    const auto width=a.at("pulse_widths").doubles();
    const auto branches=amp.size(),n=(time.size()-1)/f.stride+1,m=model.liftedDimension();
    std::vector<Complex> z(branches*n*m);std::vector<double> xi(branches*n*2),t(n),e(branches*n);
    double parity=0;
    for(std::size_t b=0;b<branches;++b) {
        std::vector<double> observations=hist;
        observations.insert(observations.end(),echo.begin()+b*time.size()+1,echo.begin()+(b+1)*time.size());
        const auto projected=core::projectHistory(observations,f.basis);
        for(std::size_t k=0;k<n;++k) {
            const auto row=hist.size()-f.basis.delays+k*f.stride;
            const auto coords=projected.coordinates.data()+row*f.basis.rank;
            const auto value=model.evaluateCoordinates(coords);
            for(std::size_t j=0;j<m;++j)z[(b*n+k)*m+j]=value[j];
            for(std::size_t j=0;j<2;++j)xi[(b*n+k)*2+j]=coords[j];
            t[k]=time[k*f.stride];e[b*n+k]=echo[b*time.size()+k*f.stride];
        }
        if(b<2) {
            core::CausalHistoryProjector online(f.basis);
            for(std::size_t k=0;k<observations.size();++k) {
                online.append(observations[k]);
                if(k+1<hist.size() || (k+1-hist.size())%4000)continue;
                const auto coords=online.projectLatest();const auto row=(k+1-hist.size())/f.stride;
                for(std::size_t j=0;j<2;++j)parity=std::max(parity,std::abs((*coords)[j]-xi[(b*n+row)*2+j]));
            }
        }
    }
    if(parity>1e-10)throw std::runtime_error("online/offline pulse projection mismatch");
    auto Z=[&](std::size_t b,std::size_t k,std::size_t j){return z[(b*n+k)*m+j];};
    auto G=[&](std::size_t plus,std::size_t k,std::size_t j){return (Z(plus,k,j)-Z(plus+1,k,j))/(2*amp[plus]);};
    const double span=(f.basis.delays-1)*f.sensorStep;
    const auto clean=std::size_t(std::ceil((span+.01)/f.controlStep-1e-9));
    Json::Array modes,widthChecks;
    std::vector<Complex> response(n*m),predicted(n*m),postClean(n*m);
    for(std::size_t j=0;j<m;++j) {
        const auto l=model.continuousEigenvalues()[j],mu=std::exp(l*f.controlStep),g1=G(1,1,j);
        double peak=0,peakTime=0,halfErr=0,doubleErr=0,gEnergy=0,evenErr=0,convErr=0;
        for(std::size_t k=0;k<n;++k) {
            const auto g=G(1,k,j);response[k*m+j]=g;
            predicted[k*m+j]=k?std::exp(l*(t[k]-f.controlStep))*g1:Complex{};
            postClean[k*m+j]=k>=clean?std::exp(l*(t[k]-t[clean]))*G(1,clean,j):Complex{};
            if(std::abs(g)>peak){peak=std::abs(g);peakTime=t[k];}
            halfErr+=std::norm(G(3,k,j)-g);doubleErr+=std::norm(G(5,k,j)-g);gEnergy+=std::norm(g);
            evenErr+=std::norm((Z(1,k,j)+Z(2,k,j)-2.*Z(0,k,j))/(2*p.delta));
            if(p.refine)convErr+=std::norm(G(branches-2,k,j)-g);
        }
        Json::Array bands;
        for(auto band:std::vector<std::pair<double,double>>{{.01,.05},{.05,.5},{.5,4.34},{4.34,8.67},{8.67,14.01}}) {
            double er=0,en=0,drift=0,znorm=0,after=0,afterEn=0;
            for(std::size_t k=1;k<n;++k)if(t[k]>=band.first-1e-9&&t[k]<band.second-1e-9) {
                const auto g=G(1,k,j);er+=std::norm(g-predicted[k*m+j]);en+=std::norm(g);
                drift+=std::norm(Z(0,k,j)-mu*Z(0,k-1,j));znorm+=std::norm(Z(0,k-1,j));
                if(k>=clean){after+=std::norm(g-postClean[k*m+j]);afterEn+=std::norm(g);}
            }
            bands.emplace_back(Json::Object{{"from",band.first},{"to",band.second},
                {"pulse_semigroup_relative_error",relative(er,en)},
                {"zero_branch_pair_error",relative(drift,znorm)},
                {"post_clean_response_semigroup_error",afterEn>0?Json(relative(after,afterEn)):Json()}});
        }
        modes.emplace_back(Json::Object{{"mode",double(j+1)},{"immediate_sensitivity_magnitude",std::abs(g1)},
            {"peak_sensitivity_magnitude",peak},{"time_to_peak",peakTime},{"peak_over_immediate",peak/std::max(std::abs(g1),1e-300)},
            {"initial_zero_drift_over_pulse_effect",std::abs(Z(0,1,j)-mu*Z(0,0,j))/std::max(std::abs(g1*p.delta),1e-300)},
            {"half_amplitude_sensitivity_change",relative(halfErr,gEnergy)},
            {"double_amplitude_sensitivity_change",relative(doubleErr,gEnergy)},
            {"even_nonlinear_over_odd_response",relative(evenErr,gEnergy)},
            {"half_integration_step_sensitivity_change",p.refine?Json(relative(convErr,gEnergy)):Json()},
            {"bands",bands}});
        if(p.extended)for(std::size_t b:{7u,9u}) {
            const auto end=std::size_t(std::llround(width[b]/f.controlStep));
            double error=0,energy=0;const auto g0=G(b,std::max(std::size_t(1),end),j);
            for(std::size_t k=std::max(std::size_t(1),end);k<n;++k) {
                const auto pred=std::exp(l*(t[k]-t[std::max(std::size_t(1),end)]))*g0;
                error+=std::norm(G(b,k,j)-pred);energy+=std::norm(G(b,k,j));
            }
            widthChecks.emplace_back(Json::Object{{"mode",double(j+1)},{"width",width[b]},
                {"post_pulse_semigroup_relative_error",relative(error,energy)}});
        }
    }
    Json result=Json::Object{{"schema","dynsyslab.echo-pulse-response/v1"},{"id",p.id},{"window",window},
        {"regime",p.regime},{"split",p.split},{"source_time",p.time},{"modes",modes},{"width_checks",widthChecks},
        {"history_span",span},{"first_clean_control_time",t[clean]},{"projection_max_absolute_error",parity}};
    core::NpzRecord out;
    out.arrays["time"]=core::NpyArray::fromDoubles({n},t);
    out.arrays["xi"]=core::NpyArray::fromDoubles({branches,n,2},xi);
    out.arrays["psi"]=core::NpyArray::fromComplex({branches,n,m},z);
    out.arrays["echo"]=core::NpyArray::fromDoubles({branches,n},e);
    out.arrays["causal_response"]=core::NpyArray::fromComplex({n,m},response);
    out.arrays["predicted_response"]=core::NpyArray::fromComplex({n,m},predicted);
    out.arrays["post_clean_prediction"]=core::NpyArray::fromComplex({n,m},postClean);
    out.arrays["input_amplitudes"]=a.at("input_amplitudes");out.arrays["pulse_widths"]=a.at("pulse_widths");
    out.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(result));
    if(!fs::exists(path))core::writeNpzAtomicNew(path,out);
    saveJson(json,result);return result;
}
}
