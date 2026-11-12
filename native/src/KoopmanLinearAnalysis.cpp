#include <dynsyslab/core/Estimation.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
namespace { bool realMode(std::complex<double> value){return std::abs(value.imag())<=1e-12;} }
std::vector<double> reconstructionSplitMatrix(
    const KoopmanModel& model,const ReconstructionResult& reconstruction) {
    if(reconstruction.modes!=model.liftedDimension()
        ||reconstruction.matrix.size()!=reconstruction.outputs*reconstruction.modes)
        throw std::invalid_argument("reconstruction does not match KROM");
    std::vector<double> result(reconstruction.outputs*model.splitDimension());
    for(std::size_t output=0;output<reconstruction.outputs;++output) {
        std::size_t split=0;
        for(std::size_t mode=0;mode<model.liftedDimension();++mode) {
            const auto coefficient=
                reconstruction.matrix[output*reconstruction.modes+mode];
            result[output*model.splitDimension()+split++]=coefficient.real();
            if(!realMode(model.continuousEigenvalues()[mode]))
                result[output*model.splitDimension()+split++]=-coefficient.imag();
        }
    }
    return result;
}


ObservabilityResult analyzeObservability(const KoopmanModel& model,
    const std::string& name,const ObservabilityOptions& options) {
    const bool discrete=options.method==ObservabilityMethod::discreteGramian
        ||options.method==ObservabilityMethod::discreteInfiniteGramian;
    if(discrete&&!(options.sampleTime>0.0))
        throw std::invalid_argument("KROM discrete observability requires sampleTime");
    const auto& map=model.reconstruction(name);
    return analyzeObservability(model.generatorSplit(),reconstructionSplitMatrix(model,map),
        model.splitDimension(),map.outputs,options);
}
} // namespace dynsyslab::core
