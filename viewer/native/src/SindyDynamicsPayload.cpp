#include "SindyDynamicsPayload.hpp"
#include <algorithm>
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
}
Payload buildSindyDynamicsPayload(const std::filesystem::path& source,const std::string& view,
    const NpzRecord& a,const json::Value& m) {
    Payload p{source.u8string(),view,m.at("label").string(),2};
    const auto d=std::size_t(m.at("coordinate_dimension").number());
    const double stop=view=="sindy_dynamics_transient"?8.:1e300;
    std::vector<std::string> names;
    for(std::size_t j=0;j<d;++j) {
        const auto& labels=m.at("coordinate_names").array();
        names.push_back(j<labels.size()?labels[j].string():"q_"+std::to_string(j+1));
    }
    const auto count=1+m.at("validation").array().size();
    for(std::size_t dataset=0;dataset<count;++dataset) {
        const auto prefix=dataset?"validation_"+std::to_string(dataset-1)+"_":"training_";
        const std::string role=dataset?"Validation "+std::to_string(dataset):"Training";
        const auto add=[&](Figure f){f.dashboardPage=int(dataset);p.figures.push_back(std::move(f));};
        for(std::size_t j=0;j<d;++j) {
            Figure f{"time",role+": derivative "+names[j],"\\tau","d("+names[j]+")/d\\tau"};
            if(a.contains(prefix+"reference_derivatives"))
                f.traces.push_back(trace(a,prefix+"reference_derivatives",prefix+"time",j,"chain-rule reference","#202a35","solid",stop));
            f.traces.push_back(trace(a,prefix+"numerical_derivatives",prefix+"time",j,"sample derivative","#e59626","dash",stop));
            f.traces.push_back(trace(a,prefix+"predicted_derivatives",prefix+"time",j,"SINDy","#009b99","dash",stop));
            add(f);
        }
        if(d==2) {
            Figure f{"phase",role+": observed coordinates",names[0],names[1]};
            f.traces.push_back(trace(a,prefix+"coordinates",prefix+"time",1,"observed","#202a35","solid",stop,true));
            add(f);
        }
        const auto& detail=dataset?m.at("validation").array()[dataset-1]:m.at("training");
        std::size_t index=0;bool plotted=false;
        for(const auto& run:detail.at("rollouts").array()) {
            if(run.at("status").string()!="complete") {
                p.warnings.push_back(role+" rollout failed at "+std::to_string(run.at("start").number())+": "+run.at("reason").string());
                ++index;continue;
            }
            if(plotted){++index;continue;}
            const auto key=prefix+"rollout_"+std::to_string(index)+"_";
            const auto rolloutTitle=role+": free rollout from tau="+std::to_string(run.at("start").number());
            const auto overlay=[&](Figure f,std::size_t j,bool phase){
                f.traces.push_back(trace(a,key+"actual",key+"time",j,"observed","#202a35","solid",1e300,phase));
                f.traces.push_back(trace(a,key+"predicted",key+"time",j,"free SINDy rollout","#009b99","dash",1e300,phase));
                add(f);
            };
            if(d==2)overlay(Figure{"phase",rolloutTitle,names[0],names[1]},1,true);
            for(std::size_t j=0;j<d;++j)
                overlay(Figure{"time",role+": rollout "+names[j],"\\tau",names[j]},j,false);
            plotted=true;++index;
        }
    }
    p.warnings.push_back("Derivative fit and free rollout are different tests. Plot thinning does not affect archived full-resolution error metrics.");
    return p;
}
}
