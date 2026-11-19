#include "SindyDynamicsPayload.hpp"
#include <algorithm>
#include <cstdio>
#include <map>
#include <stdexcept>
namespace dsl {
namespace {
using dynsyslab::core::NpzRecord;
Trace trace(const NpzRecord& record,const std::string& key,const std::string& timeKey,
    std::size_t column,const std::string& label,const std::string& color,
    const std::string& style,double stop,bool phase=false) {
    const auto& a=record.at(key);const auto t=record.at(timeKey).doubles(),v=a.doubles();
    if(a.shape.size()!=2||column>=a.shape[1]||a.shape[0]!=t.size()||(phase&&a.shape[1]<2))
        throw std::invalid_argument("SINDy plot samples must align");
    Trace result;result.label=label;result.color=color;result.lineStyle=style;
    for(std::size_t k=0;k<t.size();++k)if(t[k]<=stop) {
        result.x.push_back(phase?v[k*a.shape[1]]:t[k]);
        result.y.push_back(v[k*a.shape[1]+column]);
    }
    return result;
}
// Display labels (TeX) when recorded, else identifiers, else q_k.
std::vector<std::string> labelsOf(const json::Value& m,const char* labels,const char* names,
                                  std::size_t count,const char* fallback) {
    const auto& shown=m.at(labels).array().size()==count?m.at(labels):m.at(names);
    std::vector<std::string> out;
    for(std::size_t j=0;j<count;++j)
        out.push_back(j<shown.array().size()?shown.array()[j].string():std::string(fallback)+std::to_string(j+1));
    return out;
}
}
Payload buildSindyDynamicsPayload(const std::filesystem::path& source,const std::string& view,
    const NpzRecord& a,const json::Value& m) {
    Payload p{source.u8string(),view,m.at("label").string(),2};
    const auto d=std::size_t(m.at("coordinate_dimension").number());
    const double stop=view=="sindy_dynamics_transient"?8.:1e300;
    const auto names=labelsOf(m,"coordinate_labels","coordinate_names",d,"q_");
    // A cross-coordinate law w' = g(z) predicts the rates of the response coordinates.
    const auto r=std::size_t(m.at("response_dimension").number(double(d)));
    const bool cross=!m.at("response_source").isNull();
    const auto responses=cross?labelsOf(m,"response_coordinate_labels","response_coordinate_names",r,"w_"):names;
    const auto count=1+m.at("validation").array().size();
    // One page per kind of plot, training and each held-out record side by side.
    std::map<std::string,int> pages;
    const auto page=[&](const std::string& kind){return pages.emplace(kind,int(pages.size())).first->second;};
    for(std::size_t dataset=0;dataset<count;++dataset) {
        const auto prefix=dataset?"validation_"+std::to_string(dataset-1)+"_":"training_";
        const std::string role=dataset?"Held-out "+std::to_string(dataset):"Training";
        const auto add=[&](Figure f,const std::string& kind){f.dashboardPage=page(kind);p.figures.push_back(std::move(f));};
        const auto& detail=dataset?m.at("validation").array()[dataset-1]:m.at("training");
        // Free rollouts first: the trajectory the identified law produces on its own.
        std::size_t index=0;bool plotted=false;
        const auto& runs=detail.at("rollouts");
        for(const auto& run:(std::holds_alternative<json::Value::Array>(runs.storage())?runs.array():json::Value::Array{})) {
            const auto status=run.at("status").string();
            if(status=="diverged")  // plotted up to where the law left the solver's reach
                p.warnings.push_back(role+" rollout diverged at tau = "+std::to_string(run.at("valid_until").number())+": "+run.at("reason").string());
            else if(status!="complete") {
                p.warnings.push_back(role+" rollout failed at "+std::to_string(run.at("start").number())+": "+run.at("reason").string());
                ++index;continue;
            }
            if(plotted){++index;continue;}
            const auto key=prefix+"rollout_"+std::to_string(index)+"_";
            // Phase portraits over the whole record; time panels to `stop` (8 tau in the transient view).
            const auto overlay=[&](Figure f,const std::string& kind,const std::string& stem,std::size_t j,bool phase,
                                   const char* color,const std::string& label){
                const double until=phase?1e300:stop;
                f.traces.push_back(trace(a,key+stem+"actual",key+"time",j,"observed","#202a35","solid",until,phase));
                f.traces.push_back(trace(a,key+stem+"predicted",key+"time",j,label,color,"dash",until,phase));
                add(f,kind);
            };
            if(d>=2)overlay(Figure{"phase",role+": free SINDy rollout",names[0],names[1]},"phase","",1,true,"#009b99","free SINDy rollout");
            if(a.contains(key+"mapped_actual")) {  // the same rollout mapped, e.g. (r,v) -> (e,e')
                const auto map=run.at("map").string();
                const auto& ml=run.at("map_labels").array();
                const std::string x=ml.size()>0?ml[0].string():"e",y=ml.size()>1?ml[1].string():"e'";
                overlay(Figure{"phase",role+": rollout mapped to ("+x+", "+y+")",x,y},"mapped phase","mapped_",1,true,"#e59626","rollout, "+map);
                overlay(Figure{"time",role+": rollout mapped to "+x,"\\tau",x},"mapped time","mapped_",0,false,"#e59626","rollout, "+map);
            }
            for(std::size_t j=0;j<d;++j)
                overlay(Figure{"time",role+": rollout "+names[j],"\\tau",names[j]},"time "+std::to_string(j),"",j,false,"#009b99","free SINDy rollout");
            plotted=true;++index;
        }
        for(std::size_t j=0;j<r;++j) {
            Figure f{"time",role+": derivative "+responses[j],"\\tau","d("+responses[j]+")/d\\tau"};
            if(a.contains(prefix+"reference_derivatives"))
                f.traces.push_back(trace(a,prefix+"reference_derivatives",prefix+"time",j,"chain-rule reference","#202a35","solid",stop));
            f.traces.push_back(trace(a,prefix+"numerical_derivatives",prefix+"time",j,"sample derivative","#e59626","dash",stop));
            f.traces.push_back(trace(a,prefix+"predicted_derivatives",prefix+"time",j,"SINDy","#009b99","dash",stop));
            add(f,"derivative "+std::to_string(j));
        }
        if(d>=2&&!cross&&!plotted) {
            Figure f{"phase",role+": observed coordinates",names[0],names[1]};
            f.traces.push_back(trace(a,prefix+"coordinates",prefix+"time",1,"observed","#202a35","solid",stop,true));
            add(f,"observed phase");
        }
    }
    p.warnings.push_back("Derivative fit and free rollout are different tests. Plot thinning does not affect archived full-resolution error metrics.");
    return p;
}
}
