#include <dynsyslab/core/PeriodicTransport.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Standard coordinate transformation of control-affine fields. See Kaiser,
// Kutz and Brunton, arXiv:1707.01146; the input term generally depends on state.
// Here a periodic-family chart replaces a finite linear eigenfunction lift.
namespace dynsyslab::core {
PeriodicTransportSystem::PeriodicTransportSystem(const math::PeriodicAtlas& atlas,
    PlanarInputField field,PeriodicTransportOptions o):chart_(atlas,o.phaseSamples),options_(o){
    if(!field||o.stateChannels[0]==o.stateChannels[1]||o.stateChannels[0]>=chart_.channels()||o.stateChannels[1]>=chart_.channels()
        ||o.inputMaximumHarmonic<0||2LL*o.inputMaximumHarmonic>=static_cast<long long>(o.phaseSamples)
        ||!std::isfinite(o.frozenFrequency)||o.frozenFrequency<0)throw std::invalid_argument("Invalid periodic transport options");
    descriptor_.identifier="periodic_transport";descriptor_.label="Periodic chart transport";descriptor_.stateDimension=2;
    descriptor_.stateNames={"orbit_parameter","phase"};descriptor_.stateLabels=descriptor_.stateNames;
    configuration_.identifier=descriptor_.identifier;
    std::vector<math::PeriodicAtlasKnot> inputKnots;const double twoPi=2*std::acos(-1.);
    auto condition=[&](const auto& pulled){
        diagnostics_.minimumAbsoluteNormalizedDeterminant=std::min(diagnostics_.minimumAbsoluteNormalizedDeterminant,std::abs(pulled.normalizedDeterminant));
        diagnostics_.minimumNormalizedDeterminant=std::min(diagnostics_.minimumNormalizedDeterminant,pulled.normalizedDeterminant);
        diagnostics_.maximumNormalizedDeterminant=std::max(diagnostics_.maximumNormalizedDeterminant,pulled.normalizedDeterminant);
    };
    for(double p:chart_.parameters()){
        std::array<std::vector<std::complex<double>>,2> b;
        for(std::size_t j=0;j<o.phaseSamples;++j){
            const auto x=chart_.evaluate(p,twoPi*j/o.phaseSamples);
            const auto pulled=math::pullbackPlanarField(x,o.stateChannels,field({x.value[o.stateChannels[0]],x.value[o.stateChannels[1]]}),o.minimumNormalizedDeterminant);
            condition(pulled);
            for(int c=0;c<2;++c)b[c].emplace_back(pulled.coordinates[c]);
        }
        inputKnots.push_back({p,{math::periodicFourierSeries(b[0],chart_.frequency(p).first,o.inputMaximumHarmonic),
            math::periodicFourierSeries(b[1],chart_.frequency(p).first,o.inputMaximumHarmonic)}});
    }
    input_=std::make_unique<math::PeriodicChart>(math::PeriodicAtlas(std::move(inputKnots)),o.phaseSamples);
    diagnostics_.chartRealityProjectionL1=chart_.maximumRealityProjectionL1();
    diagnostics_.inputRealityProjectionL1=input_->maximumRealityProjectionL1();
    // Disjoint phase points at knots and between knots detect projection and
    // interpolation errors; this finite diagnostic is not a uniform bound.
    double error2=0,field2=0;
    for(std::size_t i=0;i+1<chart_.parameters().size();++i)for(double fraction:{0.,.5,1.}){
        const double p=(1-fraction)*chart_.parameters()[i]+fraction*chart_.parameters()[i+1];
        for(int j=0;j<127;++j){const double phase=twoPi*(j+.37)/127;
            const auto x=chart_.evaluate(p,phase),b=input_->evaluate(p,phase);
            const auto wanted=field({x.value[o.stateChannels[0]],x.value[o.stateChannels[1]]});
            const auto pulled=math::pullbackPlanarField(x,o.stateChannels,wanted,o.minimumNormalizedDeterminant);
            condition(pulled);
            for(int c=0;c<2;++c){const auto k=o.stateChannels[c];
                const double error=x.parameterDerivative[k]*b.value[0]+x.phaseDerivative[k]*b.value[1]-wanted[c];
                error2+=error*error;field2+=wanted[c]*wanted[c];
                if(std::abs(error)>diagnostics_.maximumInputProjectionError){diagnostics_.maximumInputProjectionError=std::abs(error);diagnostics_.worstInputLocation={p,phase};}
            }
        }
    }
    diagnostics_.inputProjectionRelativeL2=std::sqrt(error2/std::max(field2,1e-300));
    if(diagnostics_.minimumNormalizedDeterminant<0&&diagnostics_.maximumNormalizedDeterminant>0)
        throw std::domain_error("Periodic chart orientation reverses on sampled points; refine the Fourier/parameter representation or use a valid chart domain. Worst input error at parameter="
            +std::to_string(diagnostics_.worstInputLocation[0])+", phase="+std::to_string(diagnostics_.worstInputLocation[1]));
}
std::vector<double> PeriodicTransportSystem::initialState()const{return {chart_.parameters().front(),0};}
void PeriodicTransportSystem::rhs(double,const double* x,double* dx)const{dx[0]=0;dx[1]=options_.frozenFrequency>0?options_.frozenFrequency:chart_.frequency(x[0]).first;}
void PeriodicTransportSystem::rhsWithInput(double t,const double* x,double u,double* dx)const{
    if(!std::isfinite(u))throw std::invalid_argument("Nonfinite periodic input");rhs(t,x,dx);const auto b=input_->evaluate(x[0],x[1]);dx[0]+=b.value[0]*u;dx[1]+=b.value[1]*u;
}
std::vector<double> PeriodicTransportSystem::jacobian(double,const double* x)const{return {0,0,options_.frozenFrequency>0?0:chart_.frequency(x[0]).second,0};}
std::vector<double> PeriodicTransportSystem::inputMatrix(const double* x)const{return input_->evaluate(x[0],x[1]).value;}
bool PeriodicTransportSystem::invalid(double,const double* x)const{return !std::isfinite(x[0])||!std::isfinite(x[1])||x[0]<chart_.parameters().front()||x[0]>chart_.parameters().back();}
math::PeriodicChartValue PeriodicTransportSystem::readout(const double* x)const{return chart_.evaluate(x[0],x[1]);}
std::array<double,2> PeriodicTransportSystem::physicalState(const double* x)const{const auto a=readout(x);return {a.value[options_.stateChannels[0]],a.value[options_.stateChannels[1]]};}
}
