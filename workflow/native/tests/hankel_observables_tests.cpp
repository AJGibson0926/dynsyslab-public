#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/ObservableQuery.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Hankel.hpp>
#include <dynsyslab/core/System.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace dynsyslab;
namespace {
void require(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
void near(double a,double b) {
    require(std::abs(a-b)<1.e-8*(1+std::abs(b)),"unexpected numeric value");
}
template<class F> void rejects(F f) {
    bool failed=false;
    try {f();} catch (const std::exception&) {failed=true;}
    require(failed,"invalid case unexpectedly succeeded");
}
void test(const std::filesystem::path& root) {
    core::SystemConfiguration spec; spec.identifier="rpe"; auto system=core::createSystem(spec);
    const auto catalog=system->observableCatalog();
    require(catalog.timeDerivative("echo")->target=="far_field_acoustic_echo_dot","echo link");
    require(catalog.timeDerivative("echo_dot")->zeroInputOnly,"echo rate scope");
    require(catalog.timeDerivative("r")->target=="radial_velocity","radius link");
    const auto list=workflow::evaluateObservableQuery(core::json::parse(
        R"({"schema":"dynsyslab.observable-query/v1","system":{"name":"rpe"},"operation":"list"})"));
    require(!list.at("observables").array().empty(),"catalog discovery");
    std::vector<double> time,x;
    for (std::size_t i=0;i<52;++i) {
        const double t=.01*i;
        time.push_back(t); x.push_back(1.15+.06*std::sin(5*t)); x.push_back(.3*std::cos(5*t));
    }
    core::NpzRecord source;
    source.arrays["time"]=core::NpyArray::fromDoubles({52},time);
    source.arrays["values"]=core::NpyArray::fromDoubles({52,2},x);
    source.arrays["metadata_json"]=core::NpyArray::fromString(R"({"system":"rpe","archive_coordinate_system":"dimensionless"})");
    auto path=root/"source.npz";core::writeNpzAtomic(path,source);
    auto channels=core::json::parse(R"({"observables":[{"observable":"echo","taps":[0,2]},{"observable":"echo","derivative_order":1}],"window":{"delays":4,"columns":18,"sample_stride":2},"channel_scales":[2,3,4],"save_matrix":true})");
    auto job=core::json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"matched hankel","methods":[]})");
    job.at("sources")=core::json::Value::Array{path.generic_u8string()};
    auto pod=channels;pod.at("center")=true;
    auto dmd=channels;dmd.at("rank")=2.;
    job.at("methods")=core::json::Value::Array{
        core::json::Value::Object{{"name","hankel_svd"},{"settings",pod}},
        core::json::Value::Object{{"name","hdmd"},{"settings",dmd}}};
    auto result=workflow::executePostprocess(job,root);
    auto a=core::readNpz(result.archives.at(0)),b=core::readNpz(result.archives.at(1));
    auto H=a.at("hankel_matrix").doubles(),Y=b.at("shifted_hankel_matrix").doubles();
    require(H==b.at("hankel_matrix").doubles(),"SVD and DMD inputs differ");
    const auto e=catalog.at("echo").evaluateRealBatch(x,52),q=catalog.at("echo_dot").evaluateRealBatch(x,52);
    for (std::size_t j=0;j<18;++j) for (std::size_t d=0;d<4;++d) {
        const auto i=4+2*(j+d),k=3*d*18+j;
        near(H[k],e[i]/2);near(H[k+18],e[i-4]/3);near(H[k+36],q[i]/4);
        near(Y[k],e[i+2]/2);near(Y[k+18],e[i-2]/3);near(Y[k+36],q[i+2]/4);
        if(j+1<18) near(Y[k],H[k+1]);
    }
    auto t=a.at("time").doubles();near(t.front(),time[10]);near(t.back(),time[44]);
    near(b.at("training_time").doubles().back(),time[46]);
    const auto modes=a.at("pod_modes").doubles(),scores=a.at("scores").doubles();
    const auto mean=a.at("mean").doubles();const auto rank=a.at("singular_values").shape[0];
    for(std::size_t i=0;i<12;++i) for(std::size_t j=0;j<18;++j) {
        double v=mean[i];for(std::size_t k=0;k<rank;++k)v+=modes[i*rank+k]*scores[j*rank+k];
        near(v,H[i*18+j]);
    }
    auto replay=job;
    replay.at("methods")=core::json::Value::Array{
        core::json::Value::Object{{"name","hankel_matrix"},{"settings",channels}},
        core::json::Value::Object{{"name","hankel_svd"},{"settings",core::json::parse(R"({"use_previous_hankel":true,"center":true,"save_matrix":true})")}},
        core::json::Value::Object{{"name","hdmd"},{"settings",core::json::parse(R"({"use_previous_hankel":true,"rank":2,"save_matrix":true})")}}};
    const auto stored=workflow::executePostprocess(replay,root);
    for(const auto& archive:stored.archives) require(core::readNpz(archive).at("hankel_matrix").doubles()==H,"matrix replay changed H");
    const auto savedMatrix=core::readNpz(stored.archives.front());
    require(savedMatrix.at("time").doubles()==t,"stored matrix alignment changed");
    near(savedMatrix.at("successor_time").doubles().back(),time[46]);
    const auto repeated=core::readNpz(stored.archives.back());
    require(repeated.at("shifted_hankel_matrix").doubles()==Y,"stored successor changed");
    const auto eig=repeated.at("discrete_eigenvalues").complexDoubles(),expected=b.at("discrete_eigenvalues").complexDoubles();
    for(std::size_t i=0;i<eig.size();++i) near(std::abs(eig[i]-expected[i]),0.);
    source.arrays["control_values"]=core::NpyArray::fromDoubles({52,1},std::vector<double>(52,.4));
    path=root/"forced.npz";core::writeNpzAtomic(path,source);job.at("sources")=core::json::Value::Array{path.generic_u8string()};
    job.at("methods").array().resize(1);
    auto forced=core::readNpz(workflow::executePostprocess(job,root).archives.at(0));
    const auto forcedH=forced.at("hankel_matrix").doubles();
    double f[2];system->rhsWithInput(time[4],x.data()+8,.4,f);
    const auto gradient=catalog.at("echo").jacobian(x.data()+8);
    near(forcedH[36],(gradient[0].real()*f[0]+gradient[1].real()*f[1])/4);
    require(std::abs(forcedH[36]-H[36])>1e-6,"forced rate used zero-input echo-dot");
    auto& opts=job.at("methods").array().front().at("settings");
    opts=core::json::Value::Object{{"hankel_archive",stored.archives.front().generic_u8string()}};
    rejects([&]{(void)workflow::executePostprocess(job,root);}); // different source bytes
    opts=pod; opts.at("observables").array().push_back("echo");
    rejects([&]{(void)workflow::executePostprocess(job,root);});
    opts=pod; opts.at("observables")=core::json::Value::Array{core::json::Value::Object{{"observable","echo"},{"derivative_order",2.}}};
    opts.object().erase("channel_scales");
    std::vector<double> varying(52,.4);varying.back()=.3;
    source.arrays["control_values"]=core::NpyArray::fromDoubles({52,1},varying);
    path=root/"varying.npz";core::writeNpzAtomic(path,source);job.at("sources")=core::json::Value::Array{path.generic_u8string()};
    rejects([&]{(void)workflow::executePostprocess(job,root);});
    const auto tooLong=std::numeric_limits<std::size_t>::max();
    rejects([&]{core::buildHankelMatrix(x,52,2,{{0,0,1.}},{tooLong,1,2,1,0},true);});
    rejects([&]{core::buildHankelMatrix(x,52,2,{{0,52,1.}},{1},false);});
}
}
int main() {
    const auto root=std::filesystem::temp_directory_path()/("dynsyslab-hankel-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {std::filesystem::create_directories(root);test(root);std::filesystem::remove_all(root);
        std::cout<<"Hankel observable pipeline passed\n";return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
