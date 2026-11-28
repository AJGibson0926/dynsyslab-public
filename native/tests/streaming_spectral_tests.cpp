#include <dynsyslab/numerics/StreamingQr.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace {
namespace n=dynsyslab::numerics;namespace c=dynsyslab::core;namespace m=dynsyslab::math;
void check(bool b,const char*message){if(!b)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"Expected rejection");}
double relative(const std::vector<double>&a,const std::vector<double>&b){check(a.size()==b.size(),"Array size");double e=0,s=0;
    for(std::size_t i=0;i<a.size();++i){e+=std::pow(a[i]-b[i],2);s+=b[i]*b[i];}return std::sqrt(e/std::max(s,1e-30));}
void qr(){
    const std::size_t rows=91,columns=5;std::vector<double> a,w,weighted;
    for(std::size_t i=0;i<rows;++i){w.push_back(i%9==0?0:.2+i*.031);
        a.insert(a.end(),{1,std::sin(.31*i),std::cos(.53*i),std::sin(.31*i)+1e-9*std::cos(.71*i),.01*i});}
    n::WeightedQrAccumulator q(columns);
    for(std::size_t begin=0;begin<rows;){const auto count=std::min<std::size_t>(begin==0?2:17,rows-begin);
        q.append({a.begin()+begin*columns,a.begin()+(begin+count)*columns},count,{w.begin()+begin,w.begin()+begin+count});begin+=count;}
    check(q.rows()==columns&&q.samples()==rows,"QR shape bookkeeping");
    double weight=0;for(auto v:w)weight+=v;check(std::abs(q.totalWeight()-weight)<1e-12,"QR weight bookkeeping");
    for(std::size_t i=0;i<columns;++i)for(std::size_t j=0;j<columns;++j){double lhs=0,rhs=0;
        for(std::size_t k=0;k<rows;++k)lhs+=w[k]*a[k*columns+i]*a[k*columns+j];
        for(std::size_t k=0;k<q.rows();++k)rhs+=q.factor()[k*columns+i]*q.factor()[k*columns+j];
        check(std::abs(lhs-rhs)<2e-12*(1+std::abs(lhs)),"QR preserves weighted joint products");}
    weighted=a;for(std::size_t i=0;i<rows;++i)for(std::size_t j=0;j<columns;++j)weighted[i*columns+j]*=std::sqrt(w[i]);
    const auto raw=n::singularValueDecomposition(weighted,rows,columns).singularValues;
    const auto compressed=n::singularValueDecomposition(q.factor(),q.rows(),columns).singularValues;
    check(raw.back()>0&&std::abs(raw.back()-compressed.back())/raw.back()<2e-6,"QR retains a small singular direction without Gram squaring");
    const auto saved=q.factor();rejects([&]{q.append({1,2},1);});
    rejects([&]{q.append(std::vector<double>(columns,1),1,{-1});});
    rejects([&]{q.append(std::vector<double>(columns,std::numeric_limits<double>::quiet_NaN()),1);});
    check(q.factor()==saved,"Rejected QR append leaves factor intact");
}
void scalarPrediction(){
    std::vector<m::PeriodicSnapshotOrbit> bank;
    for(const auto& item:std::vector<std::pair<double,double>>{{.7,1},{1.3,3}})
        bank.push_back({{m::FourierSeries{item.first,{{-1,{.5,0}},{1,{.5,0}}}}},item.second});
    const double dt=.4,correlation=(std::cos(.7*dt)+3*std::cos(1.3*dt))/4;
    for(const std::string method:{"edmd","resdmd","mpedmd"})for(const std::string storage:{"dense","streaming_qr"}){
        c::KoopmanSpectralOptions options;options.method=method;options.sampleInterval=dt;options.acknowledgeMeasurePreserving=true;
        const auto model=c::fitPeriodicDelayModel(bank,{0},{0},32,options,true,storage);
        const auto prediction=c::rolloutCompression(*model.compression(),"coordinates",{2},5);
        for(std::size_t k=0;k<prediction.size();++k)
            check(std::abs(prediction[k]-2*std::pow(method=="mpedmd"?1:correlation,double(k)))<1e-12,
                "Independent weighted cosine regression versus unitary scalar fit");
        if(method!="mpedmd")rejects([&]{c::fitPeriodicDelaySpectrum(bank,{0},{0},32,options);});
    }
    // Fixed-space short-step limit: the cosine-only generator projection is
    // zero, although its exact evolution is nonconstant. This independently
    // checks sparse integer powers and prevents conflating a small local error
    // with an accurate fixed-time rollout (Galerkin generator: Klus et al. 2020).
    const double horizonCorrelation=(std::cos(.7)+3*std::cos(1.3))/4;
    for(std::size_t count:{5,10,20,40,100})for(const std::string method:{"edmd","mpedmd"}){
        const double h=1./count,corr=(std::cos(.7*h)+3*std::cos(1.3*h))/4;
        c::KoopmanSpectralOptions options;options.method=method;options.sampleInterval=h;options.acknowledgeMeasurePreserving=true;
        const auto model=c::fitPeriodicDelayModel(bank,{0},{0},32,options);
        const auto audit=c::auditPeriodicDelayForecast(model,bank,{count,count},{{0,1,count},{0,1,count}});
        const double power=method=="edmd"?std::pow(corr,double(count)):1.;
        const double relativeSquared=(audit[0].errorMeanSquare.back()+3*audit[1].errorMeanSquare.back())/2;
        check(std::abs(relativeSquared-(1+power*power-2*power*horizonCorrelation))<2e-12,"Exact scalar fixed-time forecast formula");
        if(count==100)check(relativeSquared>1.1,"Small native step does not fix a noninvariant trial space");
    }
}
void periodic(){
    std::vector<m::PeriodicSnapshotOrbit> bank;
    for(int j=0;j<3;++j){const double w=.7+.26*j;
        m::FourierSeries a{w,{{-2,{.1,.06}},{-1,{.5,0}},{0,{.2*j,0}},{1,{.5,0}},{2,{.1,-.06}}}};
        m::FourierSeries b{w,{{-2,{.03,0}},{-1,{0,.4}},{0,{.1,0}},{1,{0,-.4}},{2,{.03,0}}}};
        bank.push_back({{a,b},double(j*j+1)});}
    c::KoopmanSpectralOptions options;options.method="mpedmd";options.acknowledgeMeasurePreserving=true;
    options.sampleInterval=.17;options.rank=6;
    for(const auto&offsets:std::vector<std::vector<double>>{{-.34,-.17,0,.17,.34},{-.41,-.113,0,.207,.419}}){
        const auto dense=c::fitPeriodicDelaySpectrum(bank,{0,1},offsets,64,options,true,"dense");
        const auto stream=c::fitPeriodicDelaySpectrum(bank,{0,1},offsets,64,options,true,"streaming_qr");
        check(stream.compression()->samples==192&&stream.compression()->rank==dense.compression()->rank,"Original sample count and retained rank");
        const auto contract=c::json::parse(stream.compression()->coordinateContract);check(contract.at("fit_rows").number()<=20,"Bounded fitting rows");
        const auto ref=c::periodicSpectralReference(bank,bank,0,.17,{-.7,-.2,0,.17,.4},{.1,.3},.08,4);
        const auto a=c::comparePeriodicDelaySpectrum(dense,{0,1},offsets,bank,ref);
        const auto b=c::comparePeriodicDelaySpectrum(stream,{0,1},offsets,bank,ref);
        check(relative(a.density,b.density)<2e-9&&relative(a.packetRms,b.packetRms)<2e-9,"Dense/QR spectral agreement");
        check(relative(a.packetErrorRms,b.packetErrorRms)<2e-8,"Dense/QR full complex error agreement");
        const auto residuals=[](const c::KoopmanModel& model){std::vector<double> out;
            for(std::size_t j=0;j<model.spectrum().size();++j)out.push_back(model.spectrum().residual(j)->value);return out;};
        auto ar=residuals(dense),br=residuals(stream);
        std::sort(ar.begin(),ar.end());std::sort(br.begin(),br.end());check(relative(ar,br)<2e-8,"Training residual equivalence");
        const auto sampled=m::periodicSnapshotPairs(bank,64,offsets,.17);
        const auto x=m::periodicDelaySamples(bank.front().channels,64,offsets);
        check(std::equal(x.begin(),x.end(),sampled.x.begin()),"Streaming sampler and dense sampler agree");
        const std::vector<double> initial(x.begin(),x.begin()+10);
        check(relative(c::rolloutCompression(*dense.compression(),"coordinates",initial,30),c::rolloutCompression(*stream.compression(),"coordinates",initial,30))<1e-8,"Dense/QR free rollout equivalence");
        const auto forecast=c::auditPeriodicDelayForecast(stream,bank,{5,5,5});
        const std::vector<std::vector<std::size_t>> selectedSteps{{0,1,5,19,69},{0,2,69},{0,7,35}};
        const auto sparseForecast=c::auditPeriodicDelayForecast(stream,bank,{69,69,35},selectedSteps);
        const auto denseForecast=c::auditPeriodicDelayForecast(stream,bank,{69,69,35});
        for(std::size_t j=0;j<bank.size();++j){
            check(sparseForecast[j].stepIndices==selectedSteps[j],"Sparse forecast step identity");
            for(std::size_t k=0;k<selectedSteps[j].size();++k)for(std::size_t ch=0;ch<2;++ch){
                const auto index=k*2+ch,refIndex=selectedSteps[j][k]*2+ch;
                check(std::abs(sparseForecast[j].errorMeanSquare[index]-denseForecast[j].errorMeanSquare[refIndex])<2e-10,
                      "Sparse binary propagation matches ordinary repeated powers");
            }
        }
        rejects([&]{c::auditPeriodicDelayForecast(stream,bank,{5,5,5},{{0,5},{0,1,1,5},{0,5}});});
        rejects([&]{c::auditPeriodicDelayForecast(stream,bank,{5,5,5},{{0,5},{0,5},{1,5}});});
        auto past=offsets;for(auto&t:past)t-=*std::max_element(offsets.begin(),offsets.end());
        const auto latest=std::size_t(std::max_element(offsets.begin(),offsets.end())-offsets.begin());
        for(std::size_t orbit=0;orbit<bank.size();++orbit){
            const auto histories=m::periodicDelaySamples(bank[orbit].channels,64,past);
            const auto predictions=c::rolloutCompressionBatch(*stream.compression(),"coordinates",histories,5);
            std::vector<double> sampledError(12),sampledPower(12);
            for(std::size_t phase=0;phase<64;++phase)for(std::size_t k=0;k<6;++k)for(std::size_t ch=0;ch<2;++ch){
                const double theta=2*std::acos(-1.)*phase/64+bank[orbit].channels[ch].angularFrequency*.17*k;
                const double truth=m::evaluate(bank[orbit].channels[ch],theta).real();
                const double predicted=predictions[(phase*6+k)*10+latest*2+ch];
                sampledError[k*2+ch]+=std::pow(predicted-truth,2)/64;
                sampledPower[k*2+ch]+=predicted*predicted/64;
            }
            check(relative(forecast[orbit].errorMeanSquare,sampledError)<2e-8,"Parseval forecast error agrees with independent phase quadrature");
            check(relative(forecast[orbit].predictionMeanSquare,sampledPower)<2e-8,"Parseval forecast norm agrees with independent phase quadrature");
        }
    }
    rejects([&]{c::fitPeriodicDelaySpectrum(bank,{0},{0},64,options,true,"typo");});
    rejects([&]{c::fitPeriodicDelaySpectrum(bank,{0},{0},3,options,true,"streaming_qr");});
}
}
int main(){try{qr();scalarPrediction();periodic();std::cout<<"streaming spectral tests passed\n";return 0;}
catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
