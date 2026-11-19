#include <dynsyslab/core/EstimationSelection.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
// G = C Psi at every sample. Same product as reconstructObservable, kept local
// so the feedback module does not link the reconstruction fitters.
std::vector<std::complex<double>> evaluateReconstruction(const ReconstructionResult& r,
    const std::vector<std::complex<double>>& psi,std::size_t samples) {
    if(!samples||r.matrix.size()!=r.outputs*r.modes||psi.size()!=samples*r.modes)
        throw std::invalid_argument("observable reconstruction shape mismatch");
    std::vector<std::complex<double>> coefficients(r.modes*r.outputs);
    for(std::size_t mode=0;mode<r.modes;++mode)
        for(std::size_t output=0;output<r.outputs;++output)
            coefficients[mode*r.outputs+output]=r.matrix[output*r.modes+mode];
    return numerics::multiply(psi,samples,r.modes,coefficients,r.outputs);
}
}  // namespace
EstimationSubsetAssessment assessEstimationSubset(const KoopmanModel& model,
    const std::vector<double>& time,const std::vector<std::complex<double>>& psi,
    const std::vector<double>& x,const std::vector<double>& y,const EstimationSelectionOptions& o) {
    const auto m=model.eigenfunctionCount(),n=model.splitDimension(),N=time.size();
    const auto& rx=model.reconstruction(o.output);const auto& ry=model.reconstruction(o.measurement);
    if(N<2||psi.size()!=N*m||x.size()!=N*rx.outputs||y.size()!=N*ry.outputs
        ||o.outputScales.size()!=rx.outputs||!(o.measurementStandardDeviation>0)
        ||!(o.sampleTime>0)||!(o.horizon>0)||!std::isfinite(o.horizon/o.sampleTime)
        ||o.horizon/o.sampleTime>1e7)
        throw std::invalid_argument("invalid estimator-selection data or options");
    if(!std::isfinite(o.measurementStandardDeviation)||!std::isfinite(o.sampleTime)
        ||!std::isfinite(o.horizon)||!std::isfinite(o.dynamicsWeight)||o.dynamicsWeight<0
        ||!std::isfinite(o.measurementMismatchWeight)||o.measurementMismatchWeight<0)
        throw std::invalid_argument("selection weights and scales must be finite and nonnegative");
    for(auto z:psi)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))
        throw std::invalid_argument("selection Psi must be finite");
    for(const auto* values:{&x,&y})for(double v:*values)if(!std::isfinite(v))
        throw std::invalid_argument("selection outputs must be finite");
    for(double s:o.outputScales)if(!(s>0)||!std::isfinite(s))throw std::invalid_argument("output scales must be positive");
    for(std::size_t k=0;k<N;++k)if(!std::isfinite(time[k])||(k&&time[k]<=time[k-1]))throw std::invalid_argument("selection time must increase");
    EstimationSubsetAssessment result;result.modeScales.resize(m);
    for(std::size_t j=0;j<m;++j) {
        double energy=0;for(std::size_t k=0;k<N;++k)energy+=std::norm(psi[k*m+j]);
        result.modeScales[j]=std::max(1e-30,std::sqrt(energy/double(N)));
    }
    auto Cx=reconstructionSplitMatrix(model,rx),Cy=reconstructionSplitMatrix(model,ry);
    std::vector<double> splitScales;
    for(std::size_t j=0;j<m;++j) {
        splitScales.push_back(result.modeScales[j]);
        if(std::abs(model.continuousEigenvalues()[j].imag())>1e-12)splitScales.push_back(result.modeScales[j]);
    }
    for(std::size_t row=0;row<rx.outputs;++row)for(std::size_t col=0;col<n;++col)
        Cx[row*n+col]*=splitScales[col]/o.outputScales[row];
    for(std::size_t row=0;row<ry.outputs;++row)for(std::size_t col=0;col<n;++col)
        Cy[row*n+col]*=splitScales[col]/o.measurementStandardDeviation;
    const auto predictedX=evaluateReconstruction(rx,psi,N),predictedY=evaluateReconstruction(ry,psi,N);
    double yEnergy=0;
    for(std::size_t k=0;k<N;++k) {
        for(std::size_t row=0;row<rx.outputs;++row)
            result.stateReconstructionMse+=std::pow((predictedX[k*rx.outputs+row].real()-x[k*rx.outputs+row])/o.outputScales[row],2);
        for(std::size_t row=0;row<ry.outputs;++row) {
            result.measurementReconstructionMse+=std::pow(predictedY[k*ry.outputs+row].real()-y[k*ry.outputs+row],2);
            yEnergy+=y[k*ry.outputs+row]*y[k*ry.outputs+row];
        }
        if(k)for(std::size_t j=0;j<m;++j) {
            const auto prediction=model.spectrum().transition(j,time[k]-time[k-1])*psi[(k-1)*m+j];
            result.autonomousPropagationMse+=std::norm((psi[k*m+j]-prediction)/result.modeScales[j]);
        }
    }
    result.stateReconstructionMse/=double(N*rx.outputs);
    result.measurementReconstructionMse/=std::max(yEnergy,1e-30);
    result.autonomousPropagationMse/=double((N-1)*m);
    ObservabilityOptions g;g.method=ObservabilityMethod::discreteGramian;
    g.sampleTime=o.sampleTime;g.blockCount=std::max<std::size_t>(1,static_cast<std::size_t>(o.horizon/o.sampleTime));
    const auto A=model.generatorSplit();
    result.observability=analyzeObservability(A,Cy,n,ry.outputs,g);
    const auto transition=numerics::matrixExponential(A,n,o.horizon);
    const auto endpointMap=numerics::multiply(Cx,rx.outputs,n,transition,n);
    // W=V diag(s^2) V^T. Posterior covariance=(I+W)^-1. Propagate to
    // the final time before decoding x; weak directions retain their prior.
    for(std::size_t j=0;j<n;++j) {
        double energy=0;
        for(std::size_t row=0;row<rx.outputs;++row) {
            double value=0;for(std::size_t k=0;k<n;++k)
                value+=endpointMap[row*n+k]*result.observability.rightDirections[k*n+j];
            energy+=value*value;
        }
        result.posteriorStateNoiseMse+=energy/(1+result.observability.gramianEigenvalues[j]);
    }
    result.posteriorStateNoiseMse/=double(rx.outputs);
    result.score=result.stateReconstructionMse+result.posteriorStateNoiseMse
        +o.measurementMismatchWeight*result.measurementReconstructionMse
        +o.dynamicsWeight*result.autonomousPropagationMse;
    return result;
}
} // namespace dynsyslab::core
