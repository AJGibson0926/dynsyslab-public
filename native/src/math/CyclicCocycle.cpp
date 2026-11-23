#include <dynsyslab/math/CyclicCocycle.hpp>
#include <dynsyslab/math/DiscreteCycles.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#ifdef DYNSYSLAB_MATH_NUMERICS
#include <dynsyslab/numerics/Eigensystem.hpp>
#endif
// Periodic variational dynamics: Kuznetsov, Elements of Applied Bifurcation
// Theory. Character projection: Eisner et al. (2015), chapters 14/17.
namespace dynsyslab::math { namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double number(const J& x){if(!std::holds_alternative<double>(x.storage())||!std::isfinite(x.number()))throw std::invalid_argument("Finite number required");return x.number();}
std::size_t integer(const J& x){auto q=number(x);if(q<0||q>1000000||std::floor(q)!=q)throw std::invalid_argument("Bounded integer required");return std::size_t(q);}
std::vector<double> reals(const J& x){if(!std::holds_alternative<A>(x.storage())||x.array().size()>1000000)throw std::invalid_argument("Bounded array required");std::vector<double> r;for(auto&v:x.array())r.push_back(number(v));return r;}
J encode(const std::vector<double>&x){A a;for(auto v:x){if(!std::isfinite(v))throw std::overflow_error("Nonfinite cocycle result");a.emplace_back(v);}return a;}
void validate(const std::vector<double>&x,std::size_t d){if(!d||d>256||x.empty()||x.size()%d||x.size()>1000000)throw std::invalid_argument("Invalid sample shape");for(auto v:x)if(!std::isfinite(v))throw std::invalid_argument("Finite samples required");}
}
J cyclicCocycle(const std::vector<double>&x,std::size_t d,std::size_t m,const std::vector<double>&initial,std::size_t steps,std::size_t phase){
    if(!d||d>16)throw std::invalid_argument("Cocycle dimension must be 1..16");validate(x,d*d);const auto n=x.size()/(d*d);
    if(initial.size()!=d||steps>1000000/d||n*d*d*d>20000000||steps*d*d>20000000)throw std::invalid_argument("Cocycle shape/work bound exceeded");
    validate(initial,d);if(phase>=n)throw std::invalid_argument("Initial phase outside finite cycle");
    const auto projection=projectDiscreteCycle(x,d*d,m);
    const auto projected=reals(projection.at("projected_values"));
    const auto spectrum=discreteCycleSpectrum(projected,d*d);
    // Evolve the matrix family represented by its characters. The inverse DFT
    // checks the synthesis convention and keeps numerical algorithms native.
    std::vector<double> values(projected.size());double synthesisError=0,imaginaryLeakage=0;
    const auto& coefficients=spectrum.at("coefficients").array();
    for(std::size_t c=0;c<d*d;++c){
        std::vector<std::complex<double>> signal(n);
        for(std::size_t k=0;k<n;++k){const auto& pair=coefficients[k*d*d+c].array();
            signal[k]=double(n)*std::complex<double>(pair[0].number(),pair[1].number());}
        const auto decoded=core::discreteFourierTransform(signal,true);
        for(std::size_t j=0;j<n;++j){values[j*d*d+c]=decoded[j].real();
            synthesisError=std::max(synthesisError,std::abs(decoded[j].real()-projected[j*d*d+c]));
            imaginaryLeakage=std::max(imaginaryLeakage,std::abs(decoded[j].imag()));}
    }
    std::vector<double> product(d*d),next(d*d),v=initial,trajectory=initial;
    for(std::size_t j=0;j<d;++j)product[j*d+j]=1;
    for(std::size_t k=0;k<n;++k){const auto offset=((phase+k)%n)*d*d;std::fill(next.begin(),next.end(),0);
        for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j)for(std::size_t h=0;h<d;++h)next[i*d+j]+=values[offset+i*d+h]*product[h*d+j];
        product.swap(next);
    }
    for(std::size_t k=0;k<steps;++k){std::vector<double>w(d);const auto offset=((phase+k)%n)*d*d;
        for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j)w[i]+=values[offset+i*d+j]*v[j];
        v=std::move(w);trajectory.insert(trajectory.end(),v.begin(),v.end());
    }
    A multipliers;
#ifdef DYNSYSLAB_MATH_NUMERICS
    for(auto z:numerics::eigensystem(product,d,false).eigenvalues)multipliers.emplace_back(A{z.real(),z.imag()});
#endif
    double trace=0;for(std::size_t j=0;j<d;++j)trace+=product[j*d+j];
    return O{{"dimension",double(d)},{"period",double(n)},{"quotient_period",double(m)},{"initial_phase",double(phase)},
        {"projection",projection},{"character_spectrum",spectrum},
        {"character_synthesis_maximum_error",synthesisError},{"character_synthesis_imaginary_leakage",imaginaryLeakage},
        {"monodromy",encode(product)},{"multipliers",multipliers},{"trace",trace},
        {"determinant_2d",d==2?J(product[0]*product[3]-product[1]*product[2]):J()},
        {"deviations",encode(trajectory)},{"scope","Finite periodic linear cocycle. Product spans the original period, even after projection. No nonlinear closure, global stability, or infinite-limit certificate."}};
}
J cyclicPartitionGeometry(const std::vector<double>&x,std::size_t d,std::size_t base,const std::vector<double>&images){
    validate(x,d);const auto n=x.size()/d;if(base<2||base>256)throw std::invalid_argument("Partition base must be 2..256");
    std::size_t remaining=n;while(remaining>1&&remaining%base==0)remaining/=base;
    if(remaining!=1)throw std::invalid_argument("Require a base-power cycle or explicit partition_periods");
    std::vector<std::size_t> periods;for(std::size_t m=1;;m*=base){periods.push_back(m);if(m==n)break;}
    auto result=cyclicTowerPartitionGeometry(x,d,periods,images);
    auto object=result.object();object["base"]=double(base);return object;
}
J cyclicTowerPartitionGeometry(const std::vector<double>&x,std::size_t d,const std::vector<std::size_t>&periods,const std::vector<double>&images){
    validate(x,d);const auto n=x.size()/d;
    if(periods.empty()||periods.size()>32||n>4096||n*n*d*periods.size()>30000000)throw std::invalid_argument("Bounded nonempty partition hierarchy required");
    std::size_t previous=1;for(std::size_t k=0;k<periods.size();++k){const auto m=periods[k];
        if(!m||m>n||n%m||m%previous||(k&&m<=previous))throw std::invalid_argument("Partition periods must strictly increase by divisibility and divide the sample cycle");
        previous=m;
    }
    J transition;
    if(!images.empty()){
        validate(images,d);if(images.size()!=x.size())throw std::invalid_argument("Return images must match sample shape");
        double total=0,maximum=0;
        for(std::size_t j=0;j<n;++j){double squared=0;
            for(std::size_t c=0;c<d;++c){double error=images[j*d+c]-x[((j+1)%n)*d+c];squared+=error*error;}
            total+=squared;maximum=std::max(maximum,std::sqrt(squared));}
        transition=O{{"rms_error",std::sqrt(total/double(n))},{"maximum_error",maximum},
            {"scope","Supplied point images versus the next cyclic point. No set-image enclosure or trapping certificate."}};
    }
    A levels,encodedPeriods;previous=1;
    for(std::size_t level=0;level<periods.size();++level){const auto m=periods[level];std::vector<double>diameters(m,0);double separation=std::numeric_limits<double>::infinity();
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j){double squared=0;for(std::size_t c=0;c<d;++c){double e=x[i*d+c]-x[j*d+c];squared+=e*e;}double distance=std::sqrt(squared);
            if(i%m==j%m)diameters[i%m]=std::max(diameters[i%m],distance);else separation=std::min(separation,distance);
        }
        levels.emplace_back(O{{"depth",double(level)},{"cells",double(m)},{"branching_factor",double(m/previous)},{"samples_per_cell",double(n/m)},
            {"maximum_diameter",*std::max_element(diameters.begin(),diameters.end())},{"diameters",encode(diameters)},
            {"minimum_sample_separation",std::isfinite(separation)?J(separation):J()},{"singleton_terminal_level",m==n}});
        encodedPeriods.emplace_back(double(m));previous=m;
    }
    return O{{"period",double(n)},{"partition_periods",encodedPeriods},{"levels",levels},{"sample_transition",transition},{"certified_infinite_limit",false},
        {"scope","Euclidean geometry of supplied finite orbit samples grouped by time residues. Terminal zero diameter is automatic. Does not bound invariant neighborhoods, prove trapping, or establish uniform contraction at all depths."}};
}
J evaluateCyclicCocycleQuery(const J&q){const auto d=integer(q.at("dimension"));const auto x=reals(q.at("values"));
    if(q.at("operation").string()=="cycle_partition_geometry"){
        const auto images=q.at("images").isNull()?std::vector<double>{}:reals(q.at("images"));
        if(!q.at("partition_periods").isNull()){
            if(!q.at("base").isNull())throw std::invalid_argument("Choose base or partition_periods, not both");
            if(!std::holds_alternative<A>(q.at("partition_periods").storage()))throw std::invalid_argument("partition_periods must be an array");
            std::vector<std::size_t> periods;for(const auto& p:q.at("partition_periods").array())periods.push_back(integer(p));
            return cyclicTowerPartitionGeometry(x,d,periods,images);
        }
        return cyclicPartitionGeometry(x,d,q.at("base").isNull()?2:integer(q.at("base")),images);
    }
    if(!d||d>16)throw std::invalid_argument("Cocycle dimension must be 1..16");
    const auto n=x.size()/(d*d);return cyclicCocycle(x,d,q.at("quotient_period").isNull()?n:integer(q.at("quotient_period")),
        q.at("initial").isNull()?std::vector<double>(d,0):reals(q.at("initial")),q.at("steps").isNull()?0:integer(q.at("steps")),q.at("phase").isNull()?0:integer(q.at("phase")));
}
}
