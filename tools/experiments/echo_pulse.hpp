#pragma once
// Isolated pulse experiment: acquisition is the only state-aware boundary.
#include "echo_multimode.hpp"
namespace echo_pulse {
using namespace echo_multi;
struct Preparation {
    std::string id,source,regime,split;
    double time{},delta{};
    bool extended{},refine{};
};
std::vector<Preparation> preparations();
Json acquire(const fs::path& root,const fs::path& output,const Preparation&,
             const Frozen&,const core::SystemConfiguration&);
Json analyze(const fs::path& output,const Preparation&,const std::string& window,
             const Frozen&,const core::KoopmanModel&);
void fitAndControl(const fs::path& root,const fs::path& output,const std::string& window,
                   const Frozen&,const core::KoopmanModel&,bool controls,bool dense=false);
void acquireDense(const fs::path&,const fs::path&,const Frozen&,const core::SystemConfiguration&);
std::map<std::string,ModeSamples> denseSamples(const fs::path&,const fs::path&,const std::string&,
                                             const Frozen&,const core::KoopmanModel&);
Json readJson(const fs::path&);
inline double relative(double error,double energy){return std::sqrt(error/std::max(energy,1e-300));}
}
