#include <dynsyslab/core/Controllability.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <stdexcept>
#include <cmath>

namespace dynsyslab::core {
namespace {
std::vector<double> transpose(const std::vector<double>& x,std::size_t rows,std::size_t cols) {
    if(x.size()!=rows*cols)throw std::invalid_argument("controllability matrix dimensions mismatch");
    std::vector<double> result(x.size());
    for(std::size_t i=0;i<rows;++i)for(std::size_t j=0;j<cols;++j)result[j*rows+i]=x[i*cols+j];
    return result;
}
}
ControllabilityResult analyzeControllability(const std::vector<double>& A,
    const std::vector<double>& B,std::size_t n,std::size_t m,const ObservabilityOptions& options) {
    const auto dual=analyzeObservability(transpose(A,n,n),transpose(B,n,m),n,m,options);
    ControllabilityResult r;
    if(!dual.matrix.empty())r.matrix=transpose(dual.matrix,m*dual.blockCount,n);
    r.gramian=dual.gramian;r.singularValues=dual.singularValues;
    r.gramianEigenvalues=dual.gramianEigenvalues;
    r.leftDirections=dual.rightDirections;r.rightDirections=dual.leftDirections;
    r.method=dual.method;r.stateDimension=n;r.inputDimension=m;r.blockCount=dual.blockCount;
    r.rank=dual.rank;r.rankTolerance=dual.rankTolerance;r.horizon=dual.horizon;
    r.conditionNumber=dual.conditionNumber;r.resolvedConditionNumber=dual.resolvedConditionNumber;
    r.effectiveRank=dual.effectiveRank;return r;
}
ControllabilityResult analyzeControllability(const KoopmanModel& model,
    const std::string& name,const ObservabilityOptions& options) {
    const bool discrete=options.method==ObservabilityMethod::discreteGramian
        ||options.method==ObservabilityMethod::discreteInfiniteGramian;
    if(discrete) {
        const auto found=model.discreteInputMaps().find(name);
        if(found==model.discreteInputMaps().end())
            throw std::invalid_argument("discrete controllability requires a stored discrete input map");
        const auto& map=found->second;
        if(options.sampleTime!=0.0&&(!std::isfinite(options.sampleTime)
            ||std::abs(options.sampleTime-map.sampleInterval)>1e-12*map.sampleInterval))
            throw std::invalid_argument("controllability sampleTime differs from the stored discrete input map");
        auto local=options;local.sampleTime=0.0;
        // The input matrix is already discrete; do not exponentiate or scale it.
        const auto mu=model.discreteEigenvalues(map.sampleInterval);
        std::vector<double> A(model.splitDimension()*model.splitDimension());
        std::size_t j=0;
        for(std::size_t k=0;k<mu.size();++k) {
            A[j*model.splitDimension()+j]=mu[k].real();
            if(std::abs(model.continuousEigenvalues()[k].imag())>1e-12) {
                A[j*model.splitDimension()+j+1]=-mu[k].imag();
                A[(j+1)*model.splitDimension()+j]=mu[k].imag();
                A[(j+1)*model.splitDimension()+j+1]=mu[k].real();j+=2;
            } else ++j;
        }
        return analyzeControllability(A,model.realControlMap(map.matrix),model.splitDimension(),map.inputDimension,local);
    }
    const auto* map=std::get_if<ConstantControlMap>(&model.storedControlMap(name));
    if(!map)throw std::invalid_argument("LTI controllability requires a constant control map");
    return analyzeControllability(model.generatorSplit(),model.realControlMap(map->matrix),
        model.splitDimension(),map->inputs,options);
}
} // namespace dynsyslab::core
