#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace c=dynsyslab::core;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F> void rejects(F f){try{f();}catch(const std::logic_error&){return;}throw std::runtime_error("Invalid history contract accepted");}
}
int main(){std::filesystem::path path;try{
    c::HankelLeftBasis basis;basis.rows=3;basis.rank=2;
    basis.leftSingularVectors={1,0,0,1,0,0};basis.singularValues={2,4};basis.observationScale=7;
    c::FrozenHistoryTransform transform(basis,.1,{1,0},"echo");
    basis.leftSingularVectors[0]=999;
    c::CausalHistoryProjector first(transform),second(transform);
    check(&first.basis()==&second.basis(),"Frozen transform storage was duplicated");
    check(!first.push(0,2)&&!first.push(.1,8),"History became available before warmup");
    auto projected=first.push(.2,16);
    check(projected&&(*projected)[0]==2&&(*projected)[1]==1,"History selection/order/scaling changed");
    check(second.samplesSeen()==0&&!second.projectLatest(),"Two streams share a live ring buffer");
    rejects([&]{first.append(.2,1);});rejects([&]{first.append(.4,1);});rejects([&]{first.append(1);});
    check(first.samplesSeen()==3,"Rejected sample changed the ring buffer");
    const auto offline=c::projectHistory({2,8,16,32},transform.basis());
    projected=first.push(.3,32);
    check(std::abs((*projected)[0]-offline.coordinates[3])<1e-12&&
          std::abs((*projected)[1]-offline.coordinates[2])<1e-12,"FFT and online projections differ");
    first.reset();check(first.samplesSeen()==0&&!first.push(3,1),"Reset did not reset clock/warmup");
    c::NpzRecord a;c::writeFrozenHistoryTransform(a,transform);
    check(c::readFrozenHistoryTransform(a).identity()==transform.identity(),"Transform archive identity changed");
    a.arrays["singular_values"]=c::NpyArray::fromDoubles({2},{2,5});
    rejects([&]{(void)c::readFrozenHistoryTransform(a);});
    path=std::filesystem::temp_directory_path()/("dsl-transform-"+std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    c::saveFrozenHistoryTransform(path,transform);
    check(c::loadFrozenHistoryTransform(path).identity()==transform.identity(),"Standalone transform changed");
    c::FunctionDomain q;q.dimension=2;
    q.providerIdentity=transform.identity();
    c::EigenfunctionFitResult polynomial{c::FunctionExpansion({2,{1,0,0,1},{0,0},{1,1}}, {1.,{0,1}}, q), c::Spectrum::fromContinuousRates({{-1,2}}), {}};
    auto model=std::make_shared<c::KoopmanModel>(polynomial,1);
    model->setControlMap("constant",c::ConstantControlMap{{1.},1,1,"regression"});
    c::KoopmanLQRController controller(model,{1,0,0,1},{1},{},"constant");
    check(std::isfinite(controller.valueFromCoordinates(.3,projected->data())[0]),"Runtime-only history control failed");
    std::filesystem::remove(path);std::cout<<"Runtime observation contracts passed\n";return 0;
}catch(const std::exception& e){if(!path.empty()){std::error_code ec;std::filesystem::remove(path,ec);}std::cerr<<e.what()<<'\n';return 1;}}
