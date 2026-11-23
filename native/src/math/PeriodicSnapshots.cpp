#include <dynsyslab/math/PeriodicSnapshots.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Uniform Haar quadrature on each circle, explicit transverse measure weights.
// Arnold, Mathematical Methods of Classical Mechanics, chapter 10;
// delay dictionaries: Arbabi/Mezic https://doi.org/10.1137/16M1105230.
namespace dynsyslab::math {
std::vector<double> periodicDelaySamples(const std::vector<FourierSeries>&channels,
    std::size_t N,const std::vector<double>&offsets,bool split,double phase){
    const std::size_t factor=split?2:1;
    if(channels.empty()||channels.size()>1024||N<3||N>1000000||offsets.empty()
        ||offsets.size()>10000||!std::isfinite(phase))throw std::invalid_argument("Invalid periodic delay sample dimensions");
    const auto dimension=offsets.size()*channels.size()*factor;
    if(dimension>50000000/N)throw std::invalid_argument("Periodic delay sample block exceeds fifty million entries");
    const double w=channels.front().angularFrequency;
    for(const auto&s:channels)if(s.angularFrequency!=w)throw std::invalid_argument("Channels on one orbit need the same physical frequency");
    std::vector<double> out(N*dimension);
    for(std::size_t tap=0;tap<offsets.size();++tap){
        if(!std::isfinite(offsets[tap]))throw std::invalid_argument("Finite time offsets required");
        for(std::size_t c=0;c<channels.size();++c){const auto values=periodicGridValues(channels[c],N,phase+w*offsets[tap]);
            for(std::size_t i=0;i<N;++i){const auto j=i*dimension+(tap*channels.size()+c)*factor;
                if(!split&&std::abs(values[i].imag())>1e-10*(1+std::abs(values[i].real())))throw std::invalid_argument("Complex channel requires split_complex=true");
                out[j]=values[i].real();if(split)out[j+1]=values[i].imag();}
        }
    }return out;
}
PeriodicQuadratureAudit auditPeriodicQuadrature(const std::vector<PeriodicSnapshotOrbit>&orbits,
    const std::vector<std::size_t>&counts,double phase){
    if(orbits.empty()||!std::isfinite(phase))throw std::invalid_argument("Finite periodic orbit bank required");
    PeriodicQuadratureAudit out;out.channels=orbits.front().channels.size();out.phaseCounts=counts;
    if(!out.channels||out.channels>1024||counts.size()>100)throw std::invalid_argument("Quadrature audit dimensions exceed limits");
    for(auto n:counts)if(n<3||n>1000000)throw std::invalid_argument("Quadrature phase count must be 3..1000000");
    out.exactMass.resize(out.channels);out.sampledMass.resize(counts.size()*out.channels);
    for(const auto&o:orbits){
        if(o.channels.size()!=out.channels||!std::isfinite(o.weight)||o.weight<0)throw std::invalid_argument("Orbit bank shape/weight mismatch");
        out.suppliedWeight+=o.weight;
        for(std::size_t c=0;c<out.channels;++c){const auto&s=o.channels[c];
            if(s.angularFrequency!=o.channels.front().angularFrequency)throw std::invalid_argument("Channels must share their orbit frequency");
            out.exactMass[c]+=o.weight*meanSquare(s);
            for(const auto&[n,a]:s.coefficients)out.maximumHarmonic=std::max(out.maximumHarmonic,std::size_t(std::abs(n)));
            for(std::size_t k=0;k<counts.size();++k){double mass=0;for(auto v:periodicGridValues(s,counts[k],phase))mass+=std::norm(v)/counts[k];
                out.sampledMass[k*out.channels+c]+=o.weight*mass;}
        }
    }
    if(!std::isfinite(out.suppliedWeight)||out.suppliedWeight<=0)throw std::invalid_argument("Positive finite total orbit weight required");
    for(auto&v:out.exactMass)v/=out.suppliedWeight;for(auto&v:out.sampledMass)v/=out.suppliedWeight;
    out.aliasFreePhaseCount=std::max<std::size_t>(3,2*out.maximumHarmonic+1);return out;
}
FourierSeries periodicDelayReadout(const std::vector<FourierSeries>&channels,
    const std::vector<double>&offsets,const std::vector<std::complex<double>>&coefficients){
    if(channels.empty()||offsets.empty()||coefficients.size()!=channels.size()*offsets.size())
        throw std::invalid_argument("Delay readout shape mismatch");
    const auto w=channels.front().angularFrequency;FourierSeries out{w,{}};
    for(double t:offsets)if(!std::isfinite(t))throw std::invalid_argument("Finite offsets required");
    for(auto c:coefficients)if(!std::isfinite(c.real())||!std::isfinite(c.imag()))throw std::invalid_argument("Finite readout required");
    const double step=offsets.size()>1?offsets[1]-offsets[0]:0;
    bool uniform=true;for(std::size_t k=0;k<offsets.size();++k)
        if(std::abs(offsets[k]-(offsets[0]+k*step))>1e-12*(1+std::abs(offsets[k])))uniform=false;
    // Fourier translation theorem (e.g. Katznelson, Introduction to Harmonic
    // Analysis, Ch. I); Parseval then gives the exact finite-series RMS.
    for(std::size_t c=0;c<channels.size();++c){const auto&s=channels[c];meanSquare(s);
        if(s.angularFrequency!=w)throw std::invalid_argument("Readout channels need one orbit frequency");
        for(const auto&[n,a]:s.coefficients){std::complex<double> b{};
            const auto increment=std::exp(std::complex<double>(0,n*w*step));std::complex<double> phase{};
            for(std::size_t tap=0;tap<offsets.size();++tap){
                // Restart every 32 terms to bound recurrence roundoff.
                if(!uniform||tap%32==0)phase=std::exp(std::complex<double>(0,n*w*offsets[tap]));else phase*=increment;
                b+=coefficients[tap*channels.size()+c]*phase;
            }
            out.coefficients[n]+=a*b;
        }
    }
    meanSquare(out);return out;
}
PeriodicSnapshots periodicSnapshotPairs(const std::vector<PeriodicSnapshotOrbit>&orbits,
    std::size_t N,const std::vector<double>&offset,double dt,bool split,double phase){
    if(orbits.empty()||orbits.size()>100000||N<3||N>1000000||offset.empty()||offset.size()>10000
        ||!std::isfinite(dt)||dt<=0||!std::isfinite(phase))throw std::invalid_argument("Invalid periodic snapshot request");
    const auto channels=orbits[0].channels.size();const std::size_t factor=split?2:1;
    if(!channels||channels>1024||offset.size()>1000000/(channels*factor))throw std::invalid_argument("Snapshot channel count exceeds limit");
    PeriodicSnapshots out;out.dimension=offset.size()*channels*factor;
    if(orbits.size()>1000000/N)throw std::invalid_argument("Snapshot row count exceeds one million");out.samples=N*orbits.size();
    if(out.dimension>50000000/out.samples)throw std::invalid_argument("Snapshot matrix exceeds fifty million entries");
    for(double t:offset)if(!std::isfinite(t))throw std::invalid_argument("Finite time offsets required");
    out.x.resize(out.samples*out.dimension);out.y.resize(out.x.size());out.weights.resize(out.samples);double mass=0;
    for(std::size_t orbit=0;orbit<orbits.size();++orbit){const auto&o=orbits[orbit];
        if(o.channels.size()!=channels||!std::isfinite(o.weight)||o.weight<0)throw std::invalid_argument("Orbit channel/weight mismatch");mass+=o.weight;
        for(std::size_t i=0;i<N;++i)out.weights[orbit*N+i]=o.weight/double(N);
        for(unsigned step=0;step<2;++step){auto times=offset;for(auto&t:times)t+=step*dt;
            const auto values=periodicDelaySamples(o.channels,N,times,split,phase);auto&matrix=step?out.y:out.x;
            std::copy(values.begin(),values.end(),matrix.begin()+orbit*N*out.dimension);}
    }
    if(!std::isfinite(mass)||mass<=0)throw std::invalid_argument("Orbit weights must have positive finite total mass");return out;
}
}
