#include <dynsyslab/core/KoopmanArchive.hpp>
#include <chrono>
#include <iostream>
#include <stdexcept>

namespace {
namespace c=dynsyslab::core;
void check(bool ok,const char* why) { if(!ok)throw std::runtime_error(why); }
template<class F> void rejects(F f,const char* why) {
    try { f(); } catch(const std::logic_error&) { return; }
    throw std::runtime_error(why);
}
c::KoopmanModel makeModel(const c::SystemConfiguration& s) {
    return c::KoopmanModel::fromPolynomialCoefficients(s,2,{1,0,0,1},
        {1.,{0,1},2.,{0,-3}},{{-.1,2},{-.2,4}});
}
}
int main() {
    std::vector<std::filesystem::path> paths;
    const auto id=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    auto path=[&](const char* suffix) {
        paths.push_back(std::filesystem::temp_directory_path()/("dsl-binding-"+id+suffix+".npz"));
        return paths.back();
    };
    auto cleanup=[&] { for(const auto& p:paths){std::error_code ec;std::filesystem::remove(p,ec);} };
    try {
        c::SystemConfiguration request;request.identifier="rpe";request.preset="air_water_20c";
        const auto system=c::createSystem(request)->configuration();
        auto changed=system;changed.parameters["vapor_pressure_pa"]=system.parameters.at("vapor_pressure_pa")==0?2330:0;
        check(!c::systemCompatibilityDifferences(system,changed).empty(),"Vapor pressure mismatch was missed");
        auto same=system;same.initialState={1.1,.2};same.driver.amplitude=10;
        c::requireCompatibleSystems(system,same);
        auto equation=system;equation.equationsVersion="different-equation";
        rejects([&]{c::requireCompatibleSystems(system,equation);},"Equation mismatch accepted");
        auto model=makeModel(system);
        model.setControlMap("constant",c::ConstantControlMap{{1.,2.},2,1,"test"});
        const auto saved=path("model");c::saveKoopmanModel(saved,model);
        const auto loaded=c::loadKoopmanModel(saved,system);
        check(c::koopmanModelIdentity(loaded)==c::koopmanModelIdentity(model),"Identity changed on round trip");
        rejects([&]{(void)c::loadKoopmanModel(saved,changed);},"Archive silently rebound physics");
        c::ModelBindingOptions transfer;transfer.systemTransferReason="intentional transfer experiment";
        const auto transferred=c::loadKoopmanModel(saved,changed,{},transfer);
        check(transferred.systemConfiguration().parameters.at("vapor_pressure_pa")==system.parameters.at("vapor_pressure_pa"),"Transfer rewrote training physics");
        check(transferred.bindingEvidenceJson().find("explicit-system-transfer")!=std::string::npos,"Transfer not recorded");
        auto reordered=model.extract({1,0});
        rejects([&]{reordered.setControlMap("wrong",model.storedControlMap("constant"));},"Mode ordering identity ignored");
        rejects([&]{auto other=makeModel(changed);other.setControlMap("wrong",model.storedControlMap("constant"));},"Physical identity ignored");
        const auto external=path("map");c::NpzRecord map;
        map.arrays["matrix"]=c::NpyArray::fromComplex({2,1},{1.,2.});
        map.arrays["model_identity"]=c::NpyArray::fromString(c::koopmanModelIdentity(model));
        c::writeNpzAtomic(external,map);
        (void)c::loadApproximateControlMap(external,"",model);
        rejects([&]{(void)c::loadApproximateControlMap(external,"",reordered);},"External same-sized mismatched map accepted");
        const auto oldMap=path("old-map");map.arrays.erase("model_identity");c::writeNpzAtomic(oldMap,map);
        rejects([&]{(void)c::loadApproximateControlMap(oldMap,"",model);},"Unbound map accepted");
        const auto vouched=c::loadApproximateControlMap(oldMap,"",model,true);
        check(std::visit([](const auto& map){return map.provenance.find("external map")!=std::string::npos;},vouched),
                "external map provenance not recorded");
        const auto dtPath=path("legacy-dt");c::saveKoopmanModel(dtPath,model,"{\"delta_t\":0.01}");
        (void)c::loadKoopmanModel(dtPath); // map identity must reflect the serialized dt
        c::requireContinuousControlModel(model,"constant");
        model.setDiscreteInputMap("sampled",{{1.,2.},2,1,.01,"test"});
        rejects([&]{c::requireContinuousControlModel(model,"sampled");},"Sampled B implicitly converted to continuous B");
        cleanup();std::cout<<"Scientific compatibility tests passed\n";return 0;
    }catch(const std::exception& e){cleanup();std::cerr<<e.what()<<'\n';return 1;}
}
