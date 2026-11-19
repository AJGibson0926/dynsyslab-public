#include <dynsyslab_viewer/NativePayload.hpp>
#include <iostream>
#include <stdexcept>
int main(){try{
    namespace c=dynsyslab::core;
    c::NpzRecord a;const std::vector<double> t{0,1,2},q{1,0,0,1,-1,0};
    a.arrays["training_time"]=c::NpyArray::fromDoubles({3},t);
    for(const char* key:{"coordinates","reference_derivatives","numerical_derivatives","predicted_derivatives"})
        a.arrays[std::string("training_")+key]=c::NpyArray::fromDoubles({3,2},q);
    a.arrays["training_rollout_0_time"]=c::NpyArray::fromDoubles({3},t);
    for(const char* key:{"actual","predicted"})
        a.arrays[std::string("training_rollout_0_")+key]=c::NpyArray::fromDoubles({3,2},q);
    const auto meta=c::json::parse(R"({"label":"Dynamics fit","coordinate_dimension":2,"coordinate_names":["\\xi_1","\\xi_2"],"training":{"rollouts":[{"start":0,"status":"complete"}]}})");
    const auto payload=dsl::buildNativePayload("fixture.npz","sindy_dynamics",a,meta);
    // Rollout first (phase, then each coordinate), then the derivative fits; one
    // page per kind of plot. The observed-only phase panel is dropped when a
    // rollout already shows the observed path.
    if(payload.figures.size()!=5||payload.figures[0].kind!="phase"||payload.figures[0].traces.size()!=2
        ||payload.figures[3].traces.size()!=3||payload.figures[0].dashboardPage==payload.figures[3].dashboardPage)
        throw std::runtime_error("Dynamics derivative/rollout layout");
    if(payload.figures[0].xLabel!="\\xi_1"||payload.figures[1].yLabel!="\\xi_1")
        throw std::runtime_error("Coordinate labels lost");
    std::cout<<"Dynamics Viewer payload passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
