#include <dynsyslab/math/SpectralFilters.hpp>
#include <dynsyslab/math/Query.hpp>
#include <dynsyslab/math/AnalyticFourier.hpp>
#include <dynsyslab/math/PeriodicSnapshots.hpp>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <dynsyslab/core/System.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::math;
using C=std::complex<double>;
void check(bool x){if(!x)throw std::runtime_error("Spectral filter test failed");}
void near(C a,C b,double tol=2e-10){check(std::abs(a-b)<tol*(1+std::abs(b)));}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught);}
int main(){try{
    const double pi=std::acos(-1.);
    for(double eccentricity:{0.,.2,.6,.8}) {
        dynsyslab::core::SystemConfiguration configuration;configuration.identifier="kepler_planar";
        auto system=dynsyslab::core::createSystem(configuration);
        const std::vector<double> initial{1-eccentricity,0,0,std::sqrt((1+eccentricity)/(1-eccentricity))};
        const auto orbit=system->analyticSolutions()[0].instantiate(initial,0);
        check(bool(orbit));
        const auto cosine=keplerAnomalyFourier(false,eccentricity,511);
        const auto sine=keplerAnomalyFourier(true,eccentricity,511);
        for(double mean:{-.7,0.,.13,1.7,3.14,5.6}) {
            const auto state=orbit.evaluate(mean);
            near(evaluate(cosine,mean),state[0]+eccentricity,2e-12);
            near(evaluate(sine,mean),state[1]/std::sqrt(1-eccentricity*eccentricity),2e-12);
        }
    }
    rejects([]{keplerAnomalyFourier(false,1,10);});
    for(double m : {0.,1e-18,.2,.8,.999})for(const auto kind : {"sn","cn","dn"})for(int N : {0,4,40}){
        const auto f=jacobiFourier(kind,m,N,1.3);
        for(double phi : {-.8,0.,1.2,pi,2*pi-.2}){
            const auto j=dynsyslab::core::special::jacobi(phi*f.parameterPeriod/(2*pi),m);
            const double exact=std::string(kind)=="sn"?j.sn:std::string(kind)=="cn"?j.cn:j.dn;
            check(std::abs(evaluate(f.series,phi)-exact)<=f.uniformTailBound+2e-13);
        }
    }
    for(double k : {-.999,-.3,0.,.9})for(int N : {0,3,50}){
        const auto f=logarithmicCosineFourier(k,N);
        for(double phi : {0.,.7,2.,pi})check(std::abs(evaluate(f.series,phi)-std::log1p(k*std::cos(phi)))<=f.uniformTailBound+3e-13);
    }
    rejects([]{jacobiFourier("sn",1,10);});rejects([]{logarithmicCosineFourier(1,10);});
    FourierSeries fast{1.7,{{-21,.5},{21,.5},{0,.2}}};
    const auto grid=periodicGridValues(fast,8,.3);
    for(unsigned j=0;j<8;++j)near(grid[j],evaluate(fast,.3+2*pi*j/8));
    auto pairs=periodicSnapshotPairs({{{fast},2}},8,{-.4,0,.7},.13,false,.3);
    check(pairs.dimension==3&&pairs.samples==8);
    for(unsigned j=0;j<8;++j)for(unsigned k=0;k<3;++k){const double a=.3+2*pi*j/8+1.7*std::vector<double>{-.4,0,.7}[k];
        near(pairs.x[j*3+k],evaluate(fast,a));near(pairs.y[j*3+k],evaluate(fast,a+1.7*.13));near(pairs.weights[j],.25);}
    FourierSeries exponential{1,{{1,C(1)}}};
    rejects([&]{periodicSnapshotPairs({{{exponential},1}},8,{0},.1);});
    auto split=periodicSnapshotPairs({{{exponential},1}},8,{0},.1,true);
    near(C(split.y[2],split.y[3]),std::exp(C(0,2*pi/8+.1)));
    const std::vector<double> readoutOffsets{-.37,0,.81};
    FourierSeries second{1.7,{{-3,C(.2,.1)},{3,C(.2,-.1)}}};
    const std::vector<C> readoutCoefficients{{.3,.2},{-.2,0},{1,0},{0,.7},{.1,0},{.4,-.2}};
    const auto readout=periodicDelayReadout({fast,second},readoutOffsets,readoutCoefficients);
    double sampleMass=0;
    for(unsigned j=0;j<512;++j){const double p=2*pi*j/512;C direct{};
        for(unsigned k=0;k<3;++k)direct+=readoutCoefficients[2*k]*evaluate(fast,p+1.7*readoutOffsets[k])+
            readoutCoefficients[2*k+1]*evaluate(second,p+1.7*readoutOffsets[k]);
        near(evaluate(readout,p),direct);sampleMass+=std::norm(direct)/512;
    }
    near(meanSquare(readout),sampleMass);
    std::vector<double> uniformOffsets;std::vector<C> uniformCoefficients;
    for(unsigned j=0;j<128;++j){uniformOffsets.push_back((int(j)-64)*.1);uniformCoefficients.emplace_back(std::sin(.7*j),std::cos(.3*j));}
    const auto uniformReadout=periodicDelayReadout({fast},uniformOffsets,uniformCoefficients);
    for(double p:{-.7,.23,2.91}){C direct{};for(unsigned j=0;j<128;++j)
        direct+=uniformCoefficients[j]*evaluate(fast,p+1.7*uniformOffsets[j]);near(evaluate(uniformReadout,p),direct);}
    rejects([&]{periodicDelayReadout({fast,exponential},{0},{1,1});});
    for(unsigned order=1;order<=8;++order){
        for(double width:{1e-9,.05,1.,1000.})for(double delta:{0.,width/3,.17,2.9}){
            RationalCircleKernel candidate(width,order);C direct{};
            for(unsigned j=0;j<order;++j){const auto exponent=C(0,-delta)+C(0,width)*candidate.poles()[j];
                const auto rho=std::exp(exponent);const auto a=exponent.real(),b=std::remainder(exponent.imag(),2*pi);
                const C denominator(-std::expm1(a)+2*std::exp(a)*std::pow(std::sin(b/2),2),-std::exp(a)*std::sin(b));
                direct+=candidate.residues()[j]*(1.+rho)/denominator;}
            near(candidate(delta),direct.real()/(2*pi),3e-9);
        }
        RationalCircleKernel k(.08,order);
        for(unsigned n=0;n<order;++n){C moment{};
            for(unsigned j=0;j<order;++j)moment+=k.residues()[j]*std::pow(k.poles()[j],double(n));
            near(moment,n==0?1.:0.,2e-9);
        }
        C mass{},c3{};const int count=32768;
        for(int j=0;j<count;++j){double a=-pi+2*pi*j/count;
            mass+=k(a)*2*pi/count;c3+=k(a)*std::exp(C(0,3*a))*(2*pi/count);
            if(j%1024==0){near(k(a),k(-a));near(k(a),k(a+2*pi));}
        }
        near(mass,1,2e-9);near(c3,k.timeMultiplier(3),2e-9);
        near(k.timeMultiplier(0),1,2e-9);
        double tail=0;for(int lag=21;lag<2000;++lag)tail+=std::abs(k.timeMultiplier(lag))/pi;
        check(tail<=k.symmetricTimeTailBound(20)*(1+1e-12));
    }
    RationalCircleKernel poisson(.05);
    for(double a:{0.,.1,1.3})near(poisson(a),std::sinh(.05)/(2*pi*(std::cosh(.05)-std::cos(a))));
    near(poisson.timeMultiplier(-12),std::exp(-.6));
    near(poisson.symmetricTimeTailBound(20),std::exp(-.05*21)/(pi*(-std::expm1(-.05))));
    RationalCircleKernel narrow(1e-9);near(narrow(0),1/(pi*1e-9),1e-9);
    RationalCircleKernel broad(1000);near(broad(.8),1/(2*pi));
    near(RationalCircleKernel(1e308).timeMultiplier(2),0);
    FourierSeries s{.7,{{0,.3},{1,C(.4,.2)},{-1,C(.4,-.2)}}};
    auto packet=spectralPacket(s,.4,.9,poisson);
    const double phase=.2;C direct{};
    // Independent two-sided time-domain Poisson projector.
    for(int lag=-800;lag<=800;++lag)direct+=std::exp(-.05*std::abs(lag))*std::exp(C(0,-lag*.9))*evaluate(s,phase+.7*.4*lag)/(2*pi);
    near(evaluate(packet,phase),direct);
    auto measure=smoothOrbitSpectra({{s,2}},.4,{.9},poisson);
    near(measure.totalMass,2*meanSquare(s));near(measure.zeroHarmonicMass,.18);
    auto exact=predictFrequencyBands(s,.2,{0.,1.,20.},0.);
    auto quantized=predictFrequencyBands(s,.2,{0.,1.,20.},.3);
    for(int i=0;i<3;++i){double t=std::vector<double>{0.,1.,20.}[i];
        near(exact.values[i],evaluate(s,.2+.7*t));
        check(std::abs(exact.values[i]-quantized.values[i])<=quantized.knownCoefficientErrorBound[i]+1e-12);
    }
    rejects([]{RationalCircleKernel k(0);});rejects([]{RationalCircleKernel k(.1,9);});
    rejects([&]{spectralPacket(s,0,.2,poisson);});
    rejects([&]{smoothOrbitSpectra({{s,-1}},1,{0},poisson);});
    auto response=evaluateQuery(dynsyslab::core::json::parse(R"({"schema":"dynsyslab.math-query/v1","operation":"spectral_filter","action":"kernel","bandwidth":0.05,"order":1,"angles":[0],"lags":[2]})"));
    near(response.at("result").at("time_multipliers").array()[0].number(),std::exp(-.1));
    std::cout<<"Spectral filter identities passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
