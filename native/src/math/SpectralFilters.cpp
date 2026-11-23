#include <dynsyslab/math/SpectralFilters.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
void positive(double x,const char* message){if(!std::isfinite(x)||x<=0)throw std::invalid_argument(message);}
void finite(double x){if(!std::isfinite(x))throw std::invalid_argument("Finite argument required");}
}
RationalCircleKernel::RationalCircleKernel(double epsilon,std::size_t order):bandwidth_(epsilon){
    positive(epsilon,"Positive finite bandwidth required");
    if(order<1||order>8)throw std::invalid_argument("Kernel order must be 1..8");
    for(std::size_t j=0;j<order;++j)poles_.emplace_back(2.*(j+1)/(order+1)-1,1);
    // Lagrange interpolation weights at zero: sum_j alpha_j*p_j^k=delta_k0.
    // Equivalent to the Vandermonde moment equations, without a linear solver.
    for(std::size_t j=0;j<order;++j){C a=1;
        for(std::size_t k=0;k<order;++k)if(j!=k)a*=-poles_[k]/(poles_[j]-poles_[k]);
        residues_.push_back(a);
    }
    decay_=std::exp(-epsilon);oneMinusDecay_=-std::expm1(-epsilon);
    oneMinusDecaySquared_=-std::expm1(-2*epsilon);
    for(auto p:poles_)poleRotations_.push_back(std::exp(C(0,epsilon*p.real())));
}
double RationalCircleKernel::operator()(double delta) const {
    finite(delta);delta=std::remainder(delta,2*pi);double value=0;
    const C rotation=std::exp(C(0,-delta));
    for(std::size_t j=0;j<poles_.size();++j){
        // Re((alpha)(1+rho)/(1-rho)), rho=q exp(i b). Reuse fixed
        // pole rotations; no interpolation or discarded harmonics. The
        // denominator (1-q)^2+q*((1-cos b)^2+sin^2 b) avoids cancellation
        // at narrow resonances, as do expm1 evaluations of 1-q and 1-q^2.
        const C phase=rotation*poleRotations_[j];
        const double sine=phase.imag(),oneMinusCosine=1-phase.real();
        const double denominator=oneMinusDecay_*oneMinusDecay_+
            decay_*(oneMinusCosine*oneMinusCosine+sine*sine);
        value+=(residues_[j].real()*oneMinusDecaySquared_-
            2*decay_*residues_[j].imag()*sine)/denominator;
    }
    const double out=value/(2*pi);
    if(!std::isfinite(out))throw std::overflow_error("Spectral kernel overflow");return out;
}
double RationalCircleKernel::timeMultiplier(long long lag) const {
    C value{};const double k=std::abs(double(lag));
    const double decay=bandwidth_*k;
    if(decay>750)return 0; // exponential underflow, including finite-product overflow
    for(std::size_t j=0;j<poles_.size();++j)value+=residues_[j]*std::exp(C(0,decay)*poles_[j]);
    if(!std::isfinite(value.real()))throw std::overflow_error("Kernel time multiplier overflow");
    return value.real();
}
double RationalCircleKernel::symmetricTimeTailBound(std::size_t retained) const {
    double residueSum=0;for(auto a:residues_)residueSum+=std::abs(a);
    const double exponent=bandwidth_*(double(retained)+1);
    const double logarithm=std::log(residueSum)-std::log(pi)-exponent-std::log(-std::expm1(-bandwidth_));
    const double result=std::exp(logarithm);
    if(!std::isfinite(result))throw std::overflow_error("Kernel tail bound overflow");return result;
}
FourierSeries spectralPacket(const FourierSeries& series,double dt,double theta,const RationalCircleKernel& kernel){
    positive(dt,"Positive finite sample interval required");finite(theta);meanSquare(series);
    auto out=series;for(auto& [n,c]:out.coefficients)c*=kernel(theta-n*series.angularFrequency*dt);
    meanSquare(out);return out;
}
SmoothedOrbitSpectrum smoothOrbitSpectra(const std::vector<WeightedOrbitSpectrum>& orbits,
    double dt,const std::vector<double>& angles,const RationalCircleKernel& kernel){
    positive(dt,"Positive finite sample interval required");for(double a:angles)finite(a);
    if(orbits.empty()||angles.empty())throw std::invalid_argument("Orbit spectra and angles required");
    SmoothedOrbitSpectrum out;out.density.resize(angles.size());
    for(const auto& orbit:orbits){
        finite(orbit.weight);if(orbit.weight<0)throw std::invalid_argument("Nonnegative quadrature masses required");
        out.totalMass+=orbit.weight*meanSquare(orbit.series);
        for(const auto& [n,c]:orbit.series.coefficients){
            const double mass=orbit.weight*std::norm(c);
            if(n==0)out.zeroHarmonicMass+=mass;
            for(std::size_t j=0;j<angles.size();++j)out.density[j]+=mass*kernel(angles[j]-n*orbit.series.angularFrequency*dt);
        }
    }
    finite(out.totalMass);for(double d:out.density)finite(d);return out;
}
FrequencyBandPrediction predictFrequencyBands(const FourierSeries& s,double phase,
    const std::vector<double>& times,double width){
    meanSquare(s);finite(phase);finite(width);if(width<0)throw std::invalid_argument("Nonnegative frequency band width required");
    FrequencyBandPrediction out;
    for(double t:times){finite(t);C value{};double bound=0;
        for(const auto& [n,c]:s.coefficients){
            const double omega=n*s.angularFrequency;
            const double center=width>0?std::round(omega/width)*width:omega;
            finite(center);const double angle=n*phase+center*t;finite(angle);
            value+=c*std::exp(C(0,std::remainder(angle,2*pi)));
            bound+=std::abs(c)*std::min(2.,std::abs((center-omega)*t));
        }
        finite(std::abs(value));finite(bound);out.values.push_back(value);out.knownCoefficientErrorBound.push_back(bound);
    }return out;
}
}
