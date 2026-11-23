#ifdef NDEBUG
#undef NDEBUG
#endif
#include <dynsyslab/core/PeriodicFeatureIdentification.hpp>
#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
using namespace dynsyslab;
int main(){
    const double pi=std::acos(-1.);std::vector<math::PeriodicAtlasKnot> knots;
    for(double p:{1.,1.25,1.5,1.75,2.})knots.push_back({p,{{p,{{1,p/2},{-1,p/2}}},{p,{{1,{0,-p/2}},{-1,{0,p/2}}}}}});
    const math::PeriodicAtlas atlas(knots);core::PeriodicHistoryFeatureMap identity;
    identity.channels={0};identity.timeOffsets={-1.3,-.4,-.1,0};identity.outputs=4;
    identity.coefficients.assign(16,0);for(int i=0;i<4;++i)identity.coefficients[5*i]=1;
    identity.noiseStandardDeviations={.1,.2,.1,.3};
    std::vector<double> histories;for(auto pair:{std::pair<double,double>{1.42,2.3},{1.83,-1.2}})
        for(double t:identity.timeOffsets)histories.push_back(pair.first*std::cos(pair.second+pair.first*t));
    core::PeriodicFeatureIdentificationOptions options;options.phaseSeeds=64;options.starts=8;options.parameterSubdivisions=2;
    options.optimization.gradientTolerance=1e-10;
    auto plain=core::identifyPeriodicHistoryFeatures(atlas,identity,histories,{0,.7,3.},options);
    assert(plain.rank==4);assert(std::abs(plain.estimates[0].fit.parameters[0]-1.42)<2e-6);
    assert(std::abs(plain.estimates[0].fit.parameters[1]-2.3)<2e-6);
    assert(std::abs(plain.estimates[1].fit.parameters[0]-1.83)<2e-6);
    for(int k=0;k<3;++k){double t=std::vector<double>{0,.7,3.}[k];
        assert(std::abs(plain.estimates[0].forecast[2*k]-1.42*std::cos(2.3+1.42*t))<3e-6);}
    // An invertible embedding with redundant, correlated, differently scaled
    // columns preserves Gaussian information after full-covariance whitening.
    auto mixed=identity;mixed.outputs=6;mixed.coefficients={2,0,0,0,2,7, 1,3,0,0,1,0, 0,2,4,0,0,0, 0,0,1,5,0,0};
    auto projected=core::identifyPeriodicHistoryFeatures(atlas,mixed,histories,{0,.7,3.},options);
    assert(projected.rank==4);
    for(std::size_t b=0;b<plain.estimates.size();++b)for(int j=0;j<2;++j)
        assert(std::abs(plain.estimates[b].fit.parameters[j]-projected.estimates[b].fit.parameters[j])<2e-8);
    for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j){double dot=0;
        for(std::size_t k=0;k<4;++k)dot+=projected.effectiveCoefficients[i*4+k]*projected.effectiveCoefficients[j*4+k];
        const double expected=i==j?1/std::pow(identity.noiseStandardDeviations[i],2):0;
        assert(std::abs(dot-expected)<1e-10);}
    mixed.coefficients.assign(24,0);bool rejected=false;
    try{core::identifyPeriodicHistoryFeatures(atlas,mixed,histories);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    mixed=identity;mixed.timeOffsets.back()=.01;rejected=false;
    try{core::identifyPeriodicHistoryFeatures(atlas,mixed,histories);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    // Exact rotating circle: exported newest-coordinate packets equal the
    // independent analytic Fourier kernel after shifting all offsets causal.
    const double w=1.2,dt=.05;math::FourierSeries cosine{w,{{1,.5},{-1,.5}}},sine{w,{{1,{0,-.5}},{-1,{0,.5}}}};
    core::KoopmanSpectralOptions fit;fit.method="mpedmd";fit.acknowledgeMeasurePreserving=true;fit.sampleInterval=dt;
    const auto model=core::fitPeriodicDelaySpectrum({{{cosine,sine},1}},{0,1},{-.4,0,.3},64,fit);
    const std::vector<double> centers={.04,.12,.3};const auto packet=core::periodicDelayPacketFeatures(model,0,centers,.1,2);
    assert(packet.timeOffsets.back()==0);assert(packet.outputs==6);
    const math::RationalCircleKernel kernel(.1,2);
    for(double phase:{.1,2.,-2.4})for(std::size_t k=0;k<centers.size();++k){std::complex<double> value{};
        for(std::size_t i=0;i<packet.timeOffsets.size();++i)for(std::size_t c=0;c<2;++c){const auto index=2*i+c;
            const double h=c?std::sin(phase+w*packet.timeOffsets[i]):std::cos(phase+w*packet.timeOffsets[i]);
            value+=h*std::complex<double>(packet.coefficients[index*6+2*k],packet.coefficients[index*6+2*k+1]);}
        const auto expected=math::evaluate(math::spectralPacket(cosine,dt,centers[k],kernel),phase);
        assert(std::abs(value-expected)<1e-9);}
    std::cout<<"periodic feature identification tests passed\n";
}
