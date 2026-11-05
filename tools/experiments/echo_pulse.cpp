#include "echo_pulse.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace echo_pulse;
int main(int argc,char** argv) {
    try {
        if(argc<4)throw std::invalid_argument("usage: echo_pulse ROOT OUTPUT acquire|long|short|fit-long|fit-short|control-long|control-short [id prefix]");
        const fs::path root=argv[1],out=argv[2];std::string command=argv[3];fs::create_directories(out);
        const bool dense=command.rfind("dense-",0)==0;if(dense)command=command.substr(6);
        const bool shortWindow=command.find("short")!=std::string::npos;
        const std::string window=shortWindow?"short":"long";
        const auto source=root/"archives/experiments/echo-input-regression-20261006-v1"/
            (shortWindow?"short-history-0p05/refit-xi2.krom.npz":"refit-xi2.krom.npz");
        const auto old=core::loadKoopmanModel(source);const auto f=frozen(core::readNpz(source),old.eigenfunctions());
        const auto bank=core::loadKoopmanModel(root/"archives/experiments/echo-multimode-20261009-v2"/(window+"-bank.krom.npz"));
        if(command=="acquire"&&dense){acquireDense(root,out,f,old.systemConfiguration());return 0;}
        if(command=="verify") {
            const auto folder=out/"resolution-checks";fs::create_directories(folder);Json::Array checks;
            for(auto p:preparations())if(p.id=="free-train-1-p3"||p.id=="target-mild-p0"||p.id=="target-nonlinear-p0") {
                p.refine=true;p.extended=true;acquire(root,folder,p,f,old.systemConfiguration());
                for(const std::string w:{"long","short"}) {
                    const auto src=root/"archives/experiments/echo-input-regression-20261006-v1"/(w=="long"?"refit-xi2.krom.npz":"short-history-0p05/refit-xi2.krom.npz");
                    const auto original=core::loadKoopmanModel(src);const auto basis=frozen(core::readNpz(src),original.eigenfunctions());
                    const auto b=core::loadKoopmanModel(root/"archives/experiments/echo-multimode-20261009-v2"/(w+"-bank.krom.npz"));
                    auto result=analyze(folder,p,w,basis,b);
                    const auto response=core::readNpz(folder/(w+"-"+p.id+".response.npz"));
                    const auto z=response.at("psi").complexDoubles();
                    const auto n=response.at("time").doubles().size();
                    const auto branches=response.at("input_amplitudes").doubles().size();
                    Json::Array immediate;
                    for(std::size_t j=0;j<3;++j) {
                        const auto g=(z[(n+1)*3+j]-z[(2*n+1)*3+j])/(2*p.delta);
                        const auto refined=(z[((branches-2)*n+1)*3+j]-z[((branches-1)*n+1)*3+j])/(2*p.delta);
                        const auto half=(z[(3*n+1)*3+j]-z[(4*n+1)*3+j])/p.delta;
                        const auto twice=(z[(5*n+1)*3+j]-z[(6*n+1)*3+j])/(4*p.delta);
                        immediate.emplace_back(Json::Object{{"mode",double(j+1)},
                            {"refined_relative_change",std::abs(g-refined)/std::abs(g)},
                            {"half_amplitude_relative_change",std::abs(g-half)/std::abs(g)},
                            {"double_amplitude_relative_change",std::abs(g-twice)/std::abs(g)}});
                    }
                    result.object()["immediate_sensitivity_checks"]=immediate;
                    checks.push_back(result);
                }
            }
            const auto summary=folder/"summary-with-immediate.json";
            if(!fs::exists(summary))saveJson(summary,Json::Object{{"checks",checks}});return 0;
        }
        if(command.rfind("fit-",0)==0||command.rfind("control-",0)==0) {
            fitAndControl(root,out,window,f,bank,command.rfind("control-",0)==0,dense);return 0;
        }
        Json::Array rows;const auto all=preparations();
        for(const auto& p:all) {
            if(argc>=5&&p.id.rfind(argv[4],0)!=0)continue;
            std::cout<<command<<" "<<p.id<<std::endl;
            rows.push_back(command=="acquire"?acquire(root,out,p,f,old.systemConfiguration()):analyze(out,p,window,f,bank));
        }
        if(argc<5)saveJson(out/(command+"-summary.json"),Json::Object{{"preparations",rows},
            {"contract","frozen 2 latent coordinates and 3-mode banks; pulse-trained B only; no core changes"}});
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
