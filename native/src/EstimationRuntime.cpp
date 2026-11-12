#include <dynsyslab/core/EstimationRuntime.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
void finite(const std::vector<double>& values,std::size_t size,const char* name) {
    if(values.size()!=size || std::any_of(values.begin(),values.end(),
        [](double v){return !std::isfinite(v);}))
        throw std::invalid_argument(std::string(name)+" has invalid dimensions or nonfinite values");
}
std::vector<double> decode(const std::vector<double>& matrix,
    const std::vector<double>& z,const std::vector<double>& offset) {
    auto result=offset;
    for(std::size_t row=0;row<result.size();++row)
        for(std::size_t col=0;col<z.size();++col)
            result[row]+=matrix[row*z.size()+col]*z[col];
    finite(result,result.size(),"KLQE decoded values");
    return result;
}
}

KoopmanLQEObserver::KoopmanLQEObserver(const KoopmanModel& model,
    KoopmanLQEDesign design,std::string representation,std::vector<double> initial,
    std::vector<double> offsetY,std::vector<double> offsetX)
    :model_(std::make_shared<const KoopmanModel>(model)),design_(std::move(design)),
     offsetY_(std::move(offsetY)),offsetX_(std::move(offsetX)) {
    if(representation=="exact")
        throw std::invalid_argument("Measurement-only KLQE requires a fitted control map");
    requireContinuousControlModel(*model_,representation);
    const auto n=model_->splitDimension(),m=model_->inputDimension();
    const auto p=design_.solution.measurementDimension;
    const auto q=design_.outputDimension?design_.outputDimension:design_.physicalStateDimension;
    if(!n||!m||!p||design_.solution.stateDimension!=n
        ||design_.modelIdentity!=koopmanModelIdentity(*model_))
        throw std::invalid_argument("KLQE design must belong to the same ordered KROM");
    finite(design_.generator,n*n,"KLQE generator");
    finite(design_.measurementMatrix,p*n,"KLQE measurement map");
    finite(design_.stateReconstructionMatrix,q*n,"KLQE output map");
    finite(design_.solution.gain,n*p,"KLQE gain");
    if(design_.generator!=model_->generatorSplit()
        ||design_.measurementMatrix!=reconstructionSplitMatrix(
            *model_,model_->reconstruction(design_.measurementObservable))
        ||(q&&design_.stateReconstructionMatrix!=reconstructionSplitMatrix(
            *model_,model_->reconstruction(design_.stateObservable))))
        throw std::invalid_argument("KLQE design matrices differ from deployment KROM");
    finite(initial,n,"KLQE initial lifted values");
    if(offsetY_.empty())offsetY_.assign(p,0);
    if(offsetX_.empty())offsetX_.assign(q,0);
    finite(offsetY_,p,"KLQE measurement offset");
    finite(offsetX_,q,"KLQE output offset");
    estimate_.lifted=std::move(initial);
    estimate_.measurement=decode(design_.measurementMatrix,estimate_.lifted,offsetY_);
    estimate_.output=decode(design_.stateReconstructionMatrix,estimate_.lifted,offsetX_);
    baseGenerator_=design_.generator;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
        for(std::size_t a=0;a<p;++a)
            baseGenerator_[i*n+j]-=design_.solution.gain[i*p+a]*design_.measurementMatrix[a*n+j];
    const std::vector<double> zero(n,0);
    const auto inputAt=[&](const std::vector<double>& z) {
        return model_->realControlMap(model_->controlMapFromLifted(
            model_->valuesFromSplit(z),representation));
    };
    offsetInput_=inputAt(zero);
    bilinear_=std::holds_alternative<BilinearControlMap>(model_->storedControlMap(representation));
    if(bilinear_) {
        // Exact conversion of the affine complex map (including conjugates) to
        // real operators. Basis evaluations are algebra, not a statistical fit.
        inputOperators_.resize(m*n*n);
        for(std::size_t col=0;col<n;++col) {
            auto basis=zero;basis[col]=1;
            const auto map=inputAt(basis);
            for(std::size_t a=0;a<m;++a)for(std::size_t row=0;row<n;++row)
                inputOperators_[(a*n+row)*n+col]=map[row*m+a]-offsetInput_[row*m+a];
        }
    }
}

KoopmanEstimate KoopmanLQEObserver::observe(double time,
    const std::vector<double>& y,const std::vector<double>& applied) {
    const auto n=model_->splitDimension(),m=model_->inputDimension();
    const auto p=design_.solution.measurementDimension;
    finite(y,p,"KLQE measurement");finite(applied,m,"KLQE applied input");
    if(!std::isfinite(time)||(initialized_&&time<=estimate_.time))
        throw std::invalid_argument("KLQE measurement times must strictly increase");
    auto next=estimate_;
    if(initialized_) {
        const double dt=time-estimate_.time;
        const double tolerance=64*std::numeric_limits<double>::epsilon()
            *std::max(1.0,std::abs(time));
        if(heldDt_<0||std::abs(dt-heldDt_)>tolerance
            ||(bilinear_&&heldInput_!=applied)) {
            auto generator=baseGenerator_;
            if(bilinear_)for(std::size_t a=0;a<m;++a)
                for(std::size_t k=0;k<n*n;++k)
                    generator[k]+=applied[a]*inputOperators_[a*n*n+k];
            std::vector<double> input(n*(p+m));
            for(std::size_t row=0;row<n;++row) {
                for(std::size_t a=0;a<p;++a)
                    input[row*(p+m)+a]=design_.solution.gain[row*p+a];
                for(std::size_t a=0;a<m;++a)
                    input[row*(p+m)+p+a]=offsetInput_[row*m+a];
            }
            const auto held=numerics::zeroOrderHold(generator,input,n,p+m,dt);
            heldTransition_=held.transition;heldResponse_=held.response;
            heldDt_=dt;heldInput_=applied;
        }
        next.lifted.assign(n,0);
        for(std::size_t row=0;row<n;++row) {
            for(std::size_t col=0;col<n;++col)
                next.lifted[row]+=heldTransition_[row*n+col]*estimate_.lifted[col];
            for(std::size_t a=0;a<p;++a)
                next.lifted[row]+=heldResponse_[row*(p+m)+a]
                    *(previousMeasurement_[a]-offsetY_[a]);
            for(std::size_t a=0;a<m;++a)
                next.lifted[row]+=heldResponse_[row*(p+m)+p+a]*applied[a];
        }
        finite(next.lifted,n,"KLQE propagated values");
    }
    next.time=time;
    next.measurement=decode(design_.measurementMatrix,next.lifted,offsetY_);
    next.output=decode(design_.stateReconstructionMatrix,next.lifted,offsetX_);
    next.innovation.resize(p);
    for(std::size_t a=0;a<p;++a)next.innovation[a]=y[a]-next.measurement[a];
    finite(next.innovation,p,"KLQE innovation");
    estimate_=std::move(next);previousMeasurement_=y;initialized_=true;
    return estimate_;
}
} // namespace dynsyslab::core
