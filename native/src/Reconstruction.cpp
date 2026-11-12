#include <dynsyslab/core/Reconstruction.hpp>

#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

ReconstructionResult fitObservableReconstruction(
    const std::vector<std::complex<double>>& psiValues,
    std::size_t samples, std::size_t modes,
    const std::vector<std::complex<double>>& observableValues,
    std::size_t outputs, double ridge) {
    if (!samples || !modes || !outputs || psiValues.size()!=samples*modes
        || observableValues.size()!=samples*outputs || ridge<0.0 || !std::isfinite(ridge))
        throw std::invalid_argument("invalid observable-reconstruction inputs");
    numerics::LeastSquaresOptions options;
    options.ridge=ridge;
    options.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    const auto coefficients=numerics::solveLeastSquares(
        psiValues,samples,modes,observableValues,outputs,options);
    ReconstructionResult result;
    result.outputs=outputs; result.modes=modes;
    result.matrix.resize(outputs*modes);
    for (std::size_t output=0;output<outputs;++output)
        for (std::size_t mode=0;mode<modes;++mode)
            result.matrix[output*modes+mode]=coefficients[mode*outputs+output];
    const auto predicted=numerics::multiply(
        psiValues,samples,modes,coefficients,outputs);
    double residualSquared=0.0,valueSquared=0.0;
    for(std::size_t sample=0;sample<samples;++sample)
        for(std::size_t output=0;output<outputs;++output) {
            const double residual=std::abs(
                predicted[sample*outputs+output]
                    -observableValues[sample*outputs+output]);
            result.residualMaximum=std::max(result.residualMaximum,residual);
            residualSquared+=residual*residual;
            valueSquared+=std::norm(observableValues[sample*outputs+output]);
        }
    result.residualRms=std::sqrt(
        residualSquared/static_cast<double>(samples*outputs));
    result.relativeResidual=std::sqrt(residualSquared)
        /std::max(std::sqrt(valueSquared),1.0e-30);
    return result;
}

ReconstructionResult fitObservableReconstruction(
    const KoopmanModel& model, const std::vector<double>& states,
    std::size_t samples, const Observable& observable, double ridge) {
    if (observable.stateDimension()!=model.stateDimension()
        || states.size()!=samples*model.stateDimension())
        throw std::invalid_argument("KROM observable-projection state shape mismatch");
    std::vector<std::complex<double>> psi(samples*model.liftedDimension());
    std::vector<std::complex<double>> values(samples*observable.outputDimension());
    for (std::size_t sample=0;sample<samples;++sample) {
        const double* state=states.data()+sample*model.stateDimension();
        const auto lifted=model.evaluate(state), observed=observable.evaluate(state);
        std::copy(lifted.begin(),lifted.end(),psi.begin()+sample*model.liftedDimension());
        std::copy(observed.begin(),observed.end(),
                  values.begin()+sample*observable.outputDimension());
    }
    return fitObservableReconstruction(psi,samples,model.liftedDimension(),values,
                                       observable.outputDimension(),ridge);
}

std::vector<std::complex<double>> reconstructObservable(
    const ReconstructionResult& reconstruction,
    const std::vector<std::complex<double>>& psiValues,
    std::size_t samples) {
    if (!samples || reconstruction.matrix.size()!=reconstruction.outputs*reconstruction.modes
        || psiValues.size()!=samples*reconstruction.modes)
        throw std::invalid_argument("observable reconstruction shape mismatch");
    std::vector<std::complex<double>> coefficients(
        reconstruction.modes*reconstruction.outputs);
    for(std::size_t mode=0;mode<reconstruction.modes;++mode)
        for(std::size_t output=0;output<reconstruction.outputs;++output)
            coefficients[mode*reconstruction.outputs+output]=
                reconstruction.matrix[output*reconstruction.modes+mode];
    return numerics::multiply(psiValues,samples,reconstruction.modes,
                              coefficients,reconstruction.outputs);
}

ReconstructionResult fitAndSetObservableReconstruction(KoopmanModel& model,
    const std::string& name,const std::vector<std::complex<double>>& psi,
    std::size_t samples,const std::vector<std::complex<double>>& values,
    std::size_t outputs,double ridge) {
    auto result=fitObservableReconstruction(psi,samples,model.liftedDimension(),values,outputs,ridge);
    model.setReconstruction(name,result);return result;
}

ReconstructionResult fitAndSetObservableReconstruction(KoopmanModel& model,
    const std::string& name,const std::vector<double>& states,std::size_t samples,
    const Observable& observable,double ridge) {
    auto result=fitObservableReconstruction(model,states,samples,observable,ridge);
    model.setReconstruction(name,result);return result;
}

}  // namespace dynsyslab::core
