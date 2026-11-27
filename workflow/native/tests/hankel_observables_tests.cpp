#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/CoordinateInverse.hpp>
#include <algorithm>
#include <dynsyslab/workflow/ObservableQuery.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Hankel.hpp>
#include <dynsyslab/core/HankelOperator.hpp>
#include <dynsyslab/core/PathCrossings.hpp>
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
// Archives a postprocess job published for one method, in order.
std::vector<std::filesystem::path> byMethod(const workflow::PostprocessResult& r,const std::string& method) {
    std::vector<std::filesystem::path> out;
    for(std::size_t i=0;i<r.methods.size();++i) if(r.methods[i]==method) out.push_back(r.archives[i]);
    return out;
}
// The one object of `kind` a method published (methods publish data objects, primary first).
std::filesystem::path objectOf(const workflow::PostprocessResult& r,const std::string& method,const std::string& kind) {
    std::vector<std::filesystem::path> out;
    for(const auto& path:byMethod(r,method))
        if(core::dataObjectKind(core::readNpzIndex(path).metadataJson)==kind) out.push_back(path);
    require(out.size()==1,"expected one object of the kind");
    return out.front();
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
    // Named observables are a coordinate change: one shared step precedes both Hankel methods.
    require(byMethod(result,"coordinate_transform").size()==1,"shorthand expands to one coordinate change");
    require(result.methods.front()=="coordinate_transform","coordinate change runs before the Hankel module");
    auto a=core::readNpz(byMethod(result,"hankel_svd").at(0)),b=core::readNpz(byMethod(result,"hdmd").at(0));
    const auto bm=core::readNpz(objectOf(result,"hdmd","matrix"));
    auto H=a.at("hankel_matrix").doubles(),Y=bm.at("shifted_hankel_matrix").doubles();
    require(H==bm.at("hankel_matrix").doubles(),"SVD and DMD inputs differ");
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
    for(const auto* method:{"hankel_matrix","hankel_svd"})
        require(core::readNpz(byMethod(stored,method).at(0)).at("hankel_matrix").doubles()==H,"matrix replay changed H");
    require(core::readNpz(objectOf(stored,"hdmd","matrix")).at("hankel_matrix").doubles()==H,"matrix replay changed H");
    const auto savedMatrix=core::readNpz(byMethod(stored,"hankel_matrix").at(0));
    require(savedMatrix.at("time").doubles()==t,"stored matrix alignment changed");
    near(savedMatrix.at("successor_time").doubles().back(),time[46]);
    const auto repeated=core::readNpz(byMethod(stored,"hdmd").at(0));
    require(core::readNpz(objectOf(stored,"hdmd","matrix")).at("shifted_hankel_matrix").doubles()==Y,"stored successor changed");
    const auto eig=repeated.at("discrete_eigenvalues").complexDoubles(),expected=b.at("discrete_eigenvalues").complexDoubles();
    for(std::size_t i=0;i<eig.size();++i) near(std::abs(eig[i]-expected[i]),0.);
    source.arrays["control_values"]=core::NpyArray::fromDoubles({52,1},std::vector<double>(52,.4));
    path=root/"forced.npz";core::writeNpzAtomic(path,source);job.at("sources")=core::json::Value::Array{path.generic_u8string()};
    job.at("methods").array().resize(1);
    auto forced=core::readNpz(byMethod(workflow::executePostprocess(job,root),"hankel_svd").at(0));
    const auto forcedH=forced.at("hankel_matrix").doubles();
    double f[2];system->rhsWithInput(time[4],x.data()+8,.4,f);
    const auto gradient=catalog.at("echo").jacobian(x.data()+8);
    near(forcedH[36],(gradient[0].real()*f[0]+gradient[1].real()*f[1])/4);
    require(std::abs(forcedH[36]-H[36])>1e-6,"forced rate used zero-input echo-dot");
    auto& opts=job.at("methods").array().front().at("settings");
    opts=core::json::Value::Object{{"hankel_archive",byMethod(stored,"hankel_matrix").at(0).generic_u8string()}};
    rejects([&]{(void)workflow::executePostprocess(job,root);}); // different source bytes
    opts=pod; opts.at("observables").array().push_back("echo");
    rejects([&]{(void)workflow::executePostprocess(job,root);});
    opts=pod; opts.at("observables")=core::json::Value::Array{core::json::Value::Object{{"observable","echo"},{"derivative_order",2.}}};
    opts.object().erase("channel_scales");
    std::vector<double> varying(52,.4);varying.back()=.3;
    source.arrays["control_values"]=core::NpyArray::fromDoubles({52,1},varying);
    path=root/"varying.npz";core::writeNpzAtomic(path,source);job.at("sources")=core::json::Value::Array{path.generic_u8string()};
    rejects([&]{(void)workflow::executePostprocess(job,root);});
    // Explicit form of the same study: coordinate change, then Hankel on its columns.
    auto explicitJob=job;
    explicitJob.at("sources")=core::json::Value::Array{(root/"source.npz").generic_u8string()};
    auto columns=dmd;
    columns.at("observables")=core::json::parse(R"([{"observable":0,"taps":[0,2]},{"observable":1}])");
    columns.at("use_previous_output")=true;
    explicitJob.at("methods")=core::json::Value::Array{
        core::json::Value::Object{{"name","coordinate_transform"},{"settings",core::json::parse(
            R"({"kind":"channels","channels":["echo",{"observable":"echo","derivative_order":1}]})")}},
        core::json::Value::Object{{"name","hdmd"},{"settings",columns}}};
    const auto explicitRun=workflow::executePostprocess(explicitJob,root);
    require(core::readNpz(objectOf(explicitRun,"hdmd","matrix")).at("hankel_matrix").doubles()==H,"explicit form changed H");
    const auto coordinates=core::readNpz(byMethod(explicitRun,"coordinate_transform").at(0));
    require(coordinates.at("source_states").doubles()==x,"coordinate change keeps the physical states");
    {   // The eigenfunctions' interpolation stores the physical state at each training sample.
        const auto fit=core::readNpz(byMethod(explicitRun,"hdmd").at(0));
        require(!fit.contains("training_states"),"H-DMD must not copy the physical state");
        require(!fit.contains("training_coordinates")&&!fit.contains("domain__input"),"H-DMD must not copy its input");
        const auto set=core::readKoopmanModel(fit);
        const auto& interpolation=set.interpolation();
        require(interpolation.coordinates().isPhysicalState(),"the default interpolation domain is the full state");
        const auto tt=fit.at("training_time").doubles();const auto ts=interpolation.points();
        require(interpolation.coordinates().dimension==2,"physical state width");
        for(std::size_t k=0;k<tt.size();++k) {
            const auto i=static_cast<std::size_t>(std::lround(tt[k]/.01));
            near(ts[2*k],x[2*i]);near(ts[2*k+1],x[2*i+1]);
        }
    }
    // The Hankel module selects columns only; a physics name is refused there.
    auto physics=explicitJob;
    physics.at("methods").array()[1].at("settings").at("observables")=core::json::Value::Array{"echo"};
    rejects([&]{(void)workflow::executePostprocess(physics,root);});
    // Chained coordinate changes: (r,v) -> (e,e') -> xi from the echo column.
    {
        auto free=source; free.arrays.erase("control_values");
        free.arrays["metadata_json"]=core::NpyArray::fromString(R"({"system":"rpe","archive_coordinate_system":"dimensionless","forcing_type":"free"})");
        const auto freePath=root/"free.npz";core::writeNpzAtomic(freePath,free);
        const auto history=core::json::parse(R"({"kind":"history","measurement":"far_field_acoustic_echo","delays":6,"hankel_rank":3,"coordinate_rank":2,"seed":11})");
        auto chained=core::json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"chain","methods":[]})");
        chained.at("sources")=core::json::Value::Array{freePath.generic_u8string()};
        chained.at("methods")=core::json::Value::Array{core::json::Value::Object{{"name","coordinate_transform"},
            {"settings",core::json::Value::Object{{"kind","chain"},{"steps",core::json::Value::Array{
                core::json::parse(R"({"kind":"observable","observable":"echo_jet2"})"),history}}}}}};
        const auto chainRun=workflow::executePostprocess(chained,root);
        require(chainRun.archives.size()==1,"a chain publishes its last coordinates");
        const auto xi=core::readNpz(chainRun.archives.front());
        const auto meta=core::json::parse(xi.metadataJson());
        require(meta.at("coordinate_chain").array().size()==1&&meta.at("value_space").string()=="history-coordinates","chain provenance");
        auto direct=chained;
        direct.at("methods").array()[0].at("settings")=history;
        const auto once=core::readNpz(workflow::executePostprocess(direct,root).archives.front());
        const auto a1=xi.at("values").doubles(),b1=once.at("values").doubles();
        require(a1.size()==b1.size(),"chain and direct history shapes");
        for(std::size_t i=0;i<a1.size();++i)near(a1[i],b1[i]);
        const auto states=xi.at("source_states").doubles();
        require(xi.at("source_states").shape==std::vector<std::size_t>({a1.size()/2,2}),"chain carries physical states");
        near(states.front(),x[2*5]);near(states.back(),x.back());
        // The retained inverse: xi back to (e,e') by readout, callable from the CLI path.
        const auto xiMeta=core::json::parse(xi.metadataJson());
        require(xiMeta.at("coordinate_inverse").at("kind").string()=="readout","history keeps a readout inverse");
        const auto decoded=root/"decoded.npz";
        workflow::invertCoordinates(chainRun.archives.front(),chainRun.archives.front(),decoded);
        const auto back=core::readNpz(decoded).at("values").doubles();
        const auto echoes=core::readNpz(std::filesystem::u8path(meta.at("coordinate_chain").array()[0].string()));
        const auto e=echoes.at("values").doubles();const auto et=echoes.at("time").doubles(),xt=xi.at("time").doubles();
        for(std::size_t k=0;k<xt.size();k+=7){  // training pairs decode to their own targets
            const auto row=std::size_t(std::lower_bound(et.begin(),et.end(),xt[k])-et.begin());
            near(back[2*k]+1.,e[2*row]+1.);
        }
        require(core::json::parse(echoes.metadataJson()).at("coordinate_inverse").at("kind").string()=="none","echo_jet2 records no inverse");
        auto catalogMap=chained;  // a catalog map after the first step is refused
        catalogMap.at("methods").array()[0].at("settings").at("steps").array().push_back(
            core::json::parse(R"({"kind":"observable","observable":"echo_jet2"})"));
        rejects([&]{(void)workflow::executePostprocess(catalogMap,root);});
    }
    // Path crossings: a figure-eight crosses itself once, at the origin.
    {
        std::vector<double> fx,fy;
        for(int k=0;k<400;++k){const double s=.3+6.283185307179586*k/400.;fx.push_back(std::sin(s));fy.push_back(std::sin(s)*std::cos(s));}
        const auto hits=core::findPathCrossings(fx,fy);
        require(hits.size()==1,"figure-eight has one crossing");
        const auto& h=hits.front();
        // The polyline crossing is a chord approximation of the origin.
        require(std::abs(fx[h.first]+h.u*(fx[h.first+1]-fx[h.first]))<1e-3
                &&std::abs(fy[h.second]+h.v*(fy[h.second+1]-fy[h.second]))<1e-3,"crossing at the origin");
        require(core::findPathCrossings({0,1,2,3},{0,1,0,1}).empty(),"zigzag without self-intersection");
    }
    // Matrix-free realization equals the explicit matrix, across pooled segments.
    {
        const std::vector<core::HankelSegment> segments{{x,52,2},{std::vector<double>(x.begin(),x.begin()+60),30,2}};
        const std::vector<core::HankelChannel> lanes{{0,0,2.},{1,3,.5}};
        const core::HankelOperator op(segments,lanes,5);
        std::vector<double> w(op.rows()*2),v(op.columns()*2);
        for(std::size_t i=0;i<w.size();++i)w[i]=std::sin(.37*double(i)+.1);
        for(std::size_t i=0;i<v.size();++i)v[i]=std::cos(.23*double(i));
        const auto adjoint=op.applyAdjoint(w,2),applied=op.apply(v,2);
        std::vector<double> expectedApply(op.rows()*2,0.);
        std::size_t offset=0;double energy=0;
        for(const auto& segment:segments) {
            const auto m=core::buildHankelMatrix(segment.values,segment.samples,2,lanes,{5,1,1,1,0},false);
            for(std::size_t j=0;j<m.columns;++j) for(std::size_t k=0;k<2;++k) {
                double sum=0;for(std::size_t r=0;r<m.rows;++r)sum+=m.x[r*m.columns+j]*w[r*2+k];
                near(adjoint[(offset+j)*2+k],sum);
            }
            for(std::size_t r=0;r<m.rows;++r) for(std::size_t k=0;k<2;++k)
                for(std::size_t j=0;j<m.columns;++j)expectedApply[r*2+k]+=m.x[r*m.columns+j]*v[(offset+j)*2+k];
            for(double value:m.x)energy+=value*value;
            offset+=m.columns;
        }
        require(offset==op.columns(),"pooled column count");
        for(std::size_t i=0;i<applied.size();++i)near(applied[i],expectedApply[i]);
        near(op.energy(),energy);
        const auto basis=core::fitHankelLeftBasis(op,{3,4,16,7});
        for(std::size_t a2=0;a2<3;++a2) for(std::size_t b2=0;b2<3;++b2) {
            double dot=0;for(std::size_t r=0;r<op.rows();++r)dot+=basis.leftSingularVectors[r*3+a2]*basis.leftSingularVectors[r*3+b2];
            near(dot,a2==b2?1.:0.);
        }
        require(basis.singularValues[0]>=basis.singularValues[1]&&basis.cumulativeEnergy[2]<=1+1e-12,"ordered basis");
    }
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
