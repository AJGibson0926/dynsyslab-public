#include "KoopmanSpectralInternal.hpp"
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <limits>
namespace dynsyslab::core {
using namespace spectral_detail;
SpectralKROM fitKoopmanSpectral(const KoopmanSnapshotPairs& data,const SpectralDictionary& dictionary,
                               const KoopmanSpectralOptions& options) {
    if(options.method=="stochastic_resdmd")return fitStochasticKoopmanSpectral(data,dictionary,options);
    dictionary.validate();check(data,dictionary.dimension);
    if(!std::isfinite(options.sampleInterval)||options.sampleInterval<=0||
       !std::isfinite(options.relativeTolerance)||options.relativeTolerance<=0)
        throw std::invalid_argument("Positive sample interval and SVD tolerance required");
    const auto& method=options.method;
    if(method!="edmd"&&method!="resdmd"&&method!="mpedmd"&&method!="hermitian_dmd"&&
       method!="tlsdmd"&&method!="multdmd")throw std::invalid_argument("Unknown Koopman spectral estimator");
    if(method=="mpedmd"&&!options.acknowledgeMeasurePreserving)
        throw std::invalid_argument("mpEDMD requires explicit measure-preserving assumption acknowledgement");
    if(method=="hermitian_dmd"&&!options.acknowledgeSelfAdjoint)
        throw std::invalid_argument("Hermitian DMD requires explicit self-adjoint assumption acknowledgement");
    const auto m=data.samples,n=dictionary.size();
    auto w=weights(m,data.weights,options.weightPolicy);
    auto x=weighted(dictionary.evaluate(data.x),w,n),y=weighted(dictionary.evaluate(data.y),w,n);
    const auto rawX=x,rawY=y;
    if(method=="tlsdmd") {
        if(!options.rank)throw std::invalid_argument("TLS DMD requires a joint snapshot rank");
        std::vector<double> joined(m*2*n);
        for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<n;++j) {
            joined[i*2*n+j]=x[i*n+j];joined[i*2*n+n+j]=y[i*n+j];
        }
        const auto joint=numerics::singularValueDecomposition(joined,m,2*n,options.rank,options.relativeTolerance);
        if(!joint.rank)throw std::invalid_argument("TLS snapshot matrix has zero rank");
        auto qt=transpose(joint.left,m,joint.rank);
        x=product(joint.left,product(qt,x,joint.rank,m,n),m,joint.rank,n);
        y=product(joint.left,product(qt,y,joint.rank,m,n),m,joint.rank,n);
    }
    auto svd=numerics::singularValueDecomposition(x,m,n,options.rank,options.relativeTolerance);
    if(!svd.rank)throw std::invalid_argument("Source dictionary has zero numerical rank");
    SpectralKROM model;model.dictionary=dictionary;model.method=method;
    model.sampleInterval=options.sampleInterval;model.rank=svd.rank;model.samples=m;
    const auto r=model.rank;model.basis=whitening(svd);
    auto successor=product(y,model.basis,m,n,r);
    model.action=product(transpose(svd.left,m,r),successor,r,m,r);
    if(method=="mpedmd") {
        auto polar=numerics::singularValueDecomposition(model.action,r,r);
        model.action=product(polar.left,polar.rightTranspose,r,r,r);
        model.assumptions="Measure-preserving dynamics and invariant quadrature measure asserted by caller; unitary finite compression is imposed.";
    }
    if(method=="hermitian_dmd") {
        for(std::size_t i=0;i<r;++i)for(std::size_t j=i+1;j<r;++j)
            model.action[i*r+j]=model.action[j*r+i]=(model.action[i*r+j]+model.action[j*r+i])/2;
        model.assumptions="Self-adjoint Koopman action asserted by caller; symmetric finite compression is imposed.";
    }
    if(method=="multdmd") {
        if(dictionary.kind!="partition"||r!=n)throw std::invalid_argument("MultDMD requires a fully occupied, untruncated partition");
        // Multiplicative maps on a finite indicator algebra select exactly one successor cell.
        // Minimize the Gram-normalized residual objective from Bou lle/Colbrook (2024).
        auto gram=product(transpose(x,m,n),x,n,m,n),cross=product(transpose(x,m,n),y,n,m,n);
        std::vector<double> map(n*n);
        for(std::size_t i=0;i<n;++i) {
            std::size_t best=0;double cost=std::numeric_limits<double>::infinity();
            for(std::size_t j=0;j<n;++j) {
                const double value=(gram[i*n+i]-2*cross[i*n+j])/gram[j*n+j];
                if(value<cost){cost=value;best=j;}
            }
            map[i*n+best]=1;
        }
        const auto mapped=product(x,product(map,model.basis,n,n,r),m,n,r);
        model.action=product(transpose(svd.left,m,r),mapped,r,m,r);
        model.assumptions="Finite indicator algebra; one deterministic successor per cell. No measure-preservation constraint.";
    }
    if(model.assumptions.empty())model.assumptions="Finite data least-squares Koopman compression; no global convergence or nonlinear state-closure certificate.";
    if(options.weightPolicy=="smooth-birkhoff")model.assumptions+=" Smooth Birkhoff weights; accelerated asymptotics require additional ergodic regularity.";
    auto eig=numerics::eigensystem(model.action,r);
    model.multipliers=eig.eigenvalues;model.eigenfunctions.resize(n*r);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<r;++j)
        for(std::size_t k=0;k<r;++k)model.eigenfunctions[i*r+j]+=model.basis[i*r+k]*eig.rightEigenvectors[k*r+j];
    auto effective=data;effective.weights=w;
    model.spectral.residuals=koopmanResiduals(model,effective,options.residualGrid,"empirical-training");
    model.spectral.invariance=koopmanInvariance(product(rawX,model.basis,m,n,r),product(rawY,model.basis,m,n,r),m,r,options.relativeTolerance);
    model.validate();return model;
}
void fitSpectralReconstruction(SpectralKROM& model,const std::string& name,
    const std::vector<double>& points,const std::vector<double>& targets,std::size_t outputs,
    const std::vector<double>& suppliedWeights) {
    model.validate();const auto d=model.dictionary.dimension,n=model.dictionary.size(),r=model.rank;
    if(name.empty()||points.empty()||points.size()%d||!outputs||targets.size()!=points.size()/d*outputs)
        throw std::invalid_argument("Readout points/outputs shape mismatch");
    for(double v:targets)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite reconstruction target");
    const auto m=points.size()/d;auto w=weights(m,suppliedWeights);
    auto q=product(model.dictionary.evaluate(points),model.basis,m,n,r);
    auto c=numerics::solveLeastSquares(weighted(q,w,r),m,r,weighted(targets,w,outputs),outputs);
    SpectralReconstruction readout;readout.outputs=outputs;readout.coefficients=product(model.basis,c,n,r,outputs);
    auto prediction=product(q,c,m,r,outputs);double error=0,norm=0;
    for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<outputs;++j) {
        const auto k=i*outputs+j;error+=w[i]*std::pow(prediction[k]-targets[k],2);norm+=w[i]*targets[k]*targets[k];
    }
    readout.relativeTrainingError=std::sqrt(error/std::max(norm,1e-300));
    model.reconstructions[name]=std::move(readout);
}
std::vector<double> rolloutSpectralKrom(const SpectralKROM& model,const std::string& name,
    const std::vector<double>& point,std::size_t steps) {
    model.validate();if(point.size()!=model.dictionary.dimension)throw std::invalid_argument("Rollout needs one initial coordinate vector");
    return rolloutSpectralKromBatch(model,name,point,steps);
}
std::vector<double> rolloutSpectralKromBatch(const SpectralKROM& model,const std::string& name,
    const std::vector<double>& points,std::size_t steps) {
    model.validate();
    const auto dimension=model.dictionary.dimension;
    if(points.empty()||points.size()%dimension)throw std::invalid_argument("Batch rollout needs complete initial coordinate vectors");
    const auto trajectories=points.size()/dimension;
    const auto n=model.dictionary.size(),r=model.rank;const auto& readout=model.reconstructions.at(name);
    const auto all=product(model.dictionary.evaluate(points),model.basis,trajectories,n,r);
    auto c=numerics::solveLeastSquares(model.basis,n,r,readout.coefficients,readout.outputs);
    std::vector<double> out;
    if(steps==std::numeric_limits<std::size_t>::max()||
       trajectories>out.max_size()/readout.outputs/(steps+1))
        throw std::invalid_argument("Batch rollout size overflow");
    out.reserve(trajectories*(steps+1)*readout.outputs);
    // Classical Koopman finite-section propagation (same operator orientation
    // and source references as the single-trajectory API); no eigenbasis or
    // matrix logarithm is required, including for a defective compression.
    for(std::size_t trajectory=0;trajectory<trajectories;++trajectory){
        std::vector<double> q(all.begin()+trajectory*r,all.begin()+(trajectory+1)*r);
        for(std::size_t step=0;step<=steps;++step) {
            auto v=product(q,c,1,r,readout.outputs);out.insert(out.end(),v.begin(),v.end());
            if(step<steps)q=product(q,model.action,1,r,r);
            for(double value:q)if(!std::isfinite(value))throw std::runtime_error("Spectral rollout overflow");
        }
    }
    return out;
}
} // namespace dynsyslab::core
