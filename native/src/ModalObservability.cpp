#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
using Complex=std::complex<double>;
Complex exponentialIntegral(Complex s,double T) {
    const auto z=s*T;
    if(std::abs(z)<1e-5)return T*(1.0+z*(0.5+z*(1.0/6.0+z/24.0)));
    return (std::exp(z)-1.0)/s;
}
Complex finiteSum(Complex q,std::size_t count) {
    // Binary geometric summation avoids cancellation near q=1.
    Complex sum{},block=1.0,power=q,offset=1.0;
    while(count) {
        if(count&1){sum+=offset*block;offset*=power;}
        count>>=1;if(count){block*=1.0+power;power*=power;}
    }
    return sum;
}
}
ModalObservabilityResult analyzeModalObservability(const KoopmanModel& model,
    const std::string& name,const ObservabilityOptions& options,const std::vector<double>& scaling) {
    const auto& map=model.reconstruction(name);const auto n=model.liftedDimension();
    auto scales=scaling;if(scales.empty())scales.assign(n,1.0);
    if(scales.size()!=n)throw std::invalid_argument("modal coordinate scales mismatch");
    for(double s:scales)if(!(s>0.0)||!std::isfinite(s))throw std::invalid_argument("invalid modal scale");
    const auto method=options.method;
    const bool discrete=method==ObservabilityMethod::discreteGramian||method==ObservabilityMethod::discreteInfiniteGramian;
    const bool infinite=method==ObservabilityMethod::continuousInfiniteGramian||method==ObservabilityMethod::discreteInfiniteGramian;
    if(method==ObservabilityMethod::classicalMatrix)throw std::invalid_argument("modal analysis requires a Gramian method");
    if(discrete&&!(options.sampleTime>0))throw std::invalid_argument("modal discrete Gramian needs sampleTime");
    if(!discrete&&!infinite&&(!(options.horizon>0)||!std::isfinite(options.horizon)))throw std::invalid_argument("invalid Gramian horizon");
    if(!(options.relativeTolerance>0)||!std::isfinite(options.relativeTolerance))throw std::invalid_argument("invalid Gramian tolerance");
    const auto lambda=discrete?model.discreteEigenvalues(options.sampleTime):model.continuousEigenvalues();
    if(infinite)for(auto l:lambda)
        if(discrete?std::abs(l)>=1.0:l.real()>=0.0)throw std::invalid_argument("infinite modal Gramian requires strict stability");
    ModalObservabilityResult result;result.modes=n;result.coordinateScales=scales;result.gramian.resize(n*n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j) {
        Complex coefficient{};
        for(std::size_t y=0;y<map.outputs;++y)
            coefficient+=std::conj(map.matrix[y*n+i])*map.matrix[y*n+j]*scales[i]*scales[j];
        Complex value{};
        if(coefficient!=Complex{}) {
            if(discrete) {
                const auto q=std::conj(lambda[i])*lambda[j];
                value=coefficient*(infinite?1.0/(1.0-q):finiteSum(q,options.blockCount?options.blockCount:n));
            } else {
                const auto s=std::conj(lambda[i])+lambda[j];
                value=coefficient*(infinite?-1.0/s:exponentialIntegral(s,options.horizon));
            }
        }
        if(!std::isfinite(value.real())||!std::isfinite(value.imag()))throw std::runtime_error("modal Gramian overflow");
        result.gramian[i*n+j]=value;result.gramian[j*n+i]=std::conj(value);
    }
    const auto eigen=numerics::hermitianEigensystem(result.gramian,n);
    result.directions.resize(n*n);
    for(std::size_t j=0;j<n;++j) {
        const auto value=std::max(0.0,eigen.eigenvalues[n-j-1]);
        result.gramianEigenvalues.push_back(value);result.singularValues.push_back(std::sqrt(value));
        for(std::size_t i=0;i<n;++i)result.directions[i*n+j]=eigen.eigenvectors[i*n+n-j-1];
    }
    result.rankTolerance=options.relativeTolerance*double(n)*result.singularValues.front();
    for(double s:result.singularValues)if(s>result.rankTolerance)++result.rank;
    return result;
}
ModalObservabilityMatrix modalObservabilityMatrix(const KoopmanModel& model,const std::string& measurement,
    double dt,std::size_t blocks,const std::vector<std::size_t>& rows) {
    const auto& map=model.reconstruction(measurement);
    const std::size_t n=model.continuousEigenvalues().size();
    if(!n||map.matrix.size()!=map.outputs*n)throw std::invalid_argument("modal observability: readout does not match the modes");
    std::vector<std::size_t> selected=rows;
    if(selected.empty())for(std::size_t y=0;y<map.outputs;++y)selected.push_back(y);
    for(auto y:selected)if(y>=map.outputs)throw std::invalid_argument("modal observability: output row out of range");
    ModalObservabilityMatrix r;
    r.modes=n;r.outputs=selected.size();r.blocks=blocks?blocks:n;
    r.multipliers=model.discreteEigenvalues(dt);
    r.sampleTime=dt>0?dt:model.spectrum().sampleInterval();
    r.matrix.assign(r.blocks*r.outputs*n,Complex{});
    std::vector<Complex> power(n,Complex(1.0,0.0));
    for(std::size_t k=0;k<r.blocks;++k) {
        for(std::size_t y=0;y<r.outputs;++y)for(std::size_t j=0;j<n;++j)
            r.matrix[(k*r.outputs+y)*n+j]=map.matrix[selected[y]*n+j]*power[j];
        for(std::size_t j=0;j<n;++j)power[j]*=r.multipliers[j];
    }
    const auto svd=numerics::complexSingularValueDecomposition(r.matrix,r.blocks*r.outputs,n);
    r.singularValues=svd.singularValues;r.directions=svd.right;
    return r;
}
} // namespace dynsyslab::core
