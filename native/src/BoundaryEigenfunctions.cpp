#include <dynsyslab/core/BoundaryEigenfunctions.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::vector<std::vector<std::size_t>> groups(
    std::size_t modes,std::size_t outputs,
    std::vector<std::vector<std::size_t>> requested) {
    if(requested.empty()){
        if(outputs!=1)throw std::invalid_argument(
            "multi-output boundary fitting requires one mode group per output");
        requested.resize(1);requested[0].resize(modes);
        for(std::size_t mode=0;mode<modes;++mode)requested[0][mode]=mode;
    }
    if(requested.size()!=outputs)throw std::invalid_argument(
        "boundary mode-group count does not match observable outputs");
    std::set<std::size_t> unique;
    for(const auto& group:requested){
        if(group.empty())throw std::invalid_argument("boundary mode groups cannot be empty");
        for(const auto mode:group)
            if(mode>=modes||!unique.insert(mode).second)
                throw std::invalid_argument("boundary mode groups must partition the spectrum");
    }
    if(unique.size()!=modes)throw std::invalid_argument(
        "boundary mode groups must include every supplied eigenvalue");
    return requested;
}

void validate(const KoopmanTrajectoryEnsemble& source,const Observable& observable,
              double deltaTime,double tolerance) {
    if(source.samples<2||!source.trajectories||!source.stateDimension
        ||source.time.size()!=source.samples
        ||source.states.size()!=source.samples*source.trajectories*source.stateDimension
        ||observable.stateDimension()!=source.stateDimension
        ||!(deltaTime>0.0)||!(tolerance>=0.0))
        throw std::invalid_argument("boundary eigenfunction inputs are inconsistent");
    for(std::size_t index=1;index<source.samples;++index)
        if(!(source.time[index]>source.time[index-1])
            ||std::abs((source.time[index]-source.time[index-1])-deltaTime)
                >1.0e-9*std::max(1.0,deltaTime))
            throw std::invalid_argument(
                "boundary eigenfunction trajectories require one uniform time grid");
}

}  // namespace

BoundaryEigenfunctionResult fitBoundaryEigenfunctions(
    const std::vector<std::complex<double>>& continuousEigenvalues,
    double deltaTime,const KoopmanTrajectoryEnsemble& ensemble,
    const Observable& observable,const BoundaryEigenfunctionOptions& options) {
    validate(ensemble,observable,deltaTime,options.singularValueTolerance);
    if(continuousEigenvalues.empty())throw std::invalid_argument(
        "boundary eigenfunction fitting requires eigenvalues");
    const auto modeGroups=groups(continuousEigenvalues.size(),
        observable.outputDimension(),options.modeGroups);
    const std::size_t modes=continuousEigenvalues.size();
    const std::size_t flattened=ensemble.samples*ensemble.trajectories;
    std::vector<double> states(flattened*ensemble.stateDimension);
    std::vector<double> observed(flattened*observable.outputDimension());
    for(std::size_t trajectory=0;trajectory<ensemble.trajectories;++trajectory)
        for(std::size_t sample=0;sample<ensemble.samples;++sample){
            const std::size_t source=(trajectory*ensemble.samples+sample)
                *ensemble.stateDimension;
            const std::size_t target=(trajectory*ensemble.samples+sample)
                *ensemble.stateDimension;
            std::copy_n(ensemble.states.begin()+static_cast<std::ptrdiff_t>(source),
                ensemble.stateDimension,states.begin()+static_cast<std::ptrdiff_t>(target));
            const auto value=observable.evaluateReal(ensemble.states.data()+source);
            std::copy(value.begin(),value.end(),observed.begin()+static_cast<std::ptrdiff_t>(
                (trajectory*ensemble.samples+sample)*observable.outputDimension()));
        }
    std::vector<std::complex<double>> functionValues(flattened*modes);
    std::vector<double> residuals;
    ReconstructionResult reconstruction;
    reconstruction.outputs=observable.outputDimension();reconstruction.modes=modes;
    reconstruction.matrix.assign(reconstruction.outputs*modes,{});
    for(std::size_t output=0;output<observable.outputDimension();++output){
        const auto& assigned=modeGroups[output];
        std::vector<std::complex<double>> powers(ensemble.samples*assigned.size());
        for(std::size_t sample=0;sample<ensemble.samples;++sample)
            for(std::size_t column=0;column<assigned.size();++column)
                powers[sample*assigned.size()+column]=
                    std::exp(continuousEigenvalues[assigned[column]]
                        *deltaTime*static_cast<double>(sample));
        std::vector<std::complex<double>> measurement(
            ensemble.samples*ensemble.trajectories);
        for(std::size_t trajectory=0;trajectory<ensemble.trajectories;++trajectory)
            for(std::size_t sample=0;sample<ensemble.samples;++sample)
                measurement[sample*ensemble.trajectories+trajectory]=observed[
                        (trajectory*ensemble.samples+sample)*observable.outputDimension()+output];
        numerics::LeastSquaresOptions solveOptions;
        solveOptions.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
        solveOptions.singularValueTolerance=options.singularValueTolerance;
        const auto boundary=numerics::solveLeastSquares(
            powers,ensemble.samples,assigned.size(),measurement,
            ensemble.trajectories,solveOptions);
        const auto prediction=numerics::multiply(
            powers,ensemble.samples,assigned.size(),boundary,
            ensemble.trajectories);
        double residualSquared=0.0,measurementSquared=0.0;
        for(std::size_t sample=0;sample<ensemble.samples;++sample)
            for(std::size_t trajectory=0;trajectory<ensemble.trajectories;++trajectory) {
                const auto measured=measurement[sample*ensemble.trajectories+trajectory];
                residualSquared+=std::norm(measured
                    -prediction[sample*ensemble.trajectories+trajectory]);
                measurementSquared+=std::norm(measured);
            }
        residuals.push_back(std::sqrt(residualSquared)/std::max(
            std::sqrt(measurementSquared),std::numeric_limits<double>::min()));
        for(std::size_t local=0;local<assigned.size();++local){
            const std::size_t mode=assigned[local];
            reconstruction.matrix[output*modes+mode]=1.0;
            for(std::size_t trajectory=0;trajectory<ensemble.trajectories;++trajectory)
                for(std::size_t sample=0;sample<ensemble.samples;++sample)
                    functionValues[(trajectory*ensemble.samples+sample)*modes+mode]=
                        powers[sample*assigned.size()+local]
                        *boundary[local*ensemble.trajectories+trajectory];
        }
    }
    ScatteredEigenfunctionData functions;
    functions.continuousEigenvalues=continuousEigenvalues;
    functions.trainingStates=std::move(states);
    functions.trainingValues=std::move(functionValues);
    functions.samples=flattened;functions.stateDimension=ensemble.stateDimension;
    functions.neighbors=options.interpolationNeighbors;
    functions.finiteDifferenceRelativeStep=options.gradientRelativeStep;
    functions.provenance="fitted boundary values with scattered state interpolation";
    reconstruction.relativeResidual=*std::max_element(residuals.begin(),residuals.end());
    return {std::move(functions),std::move(reconstruction),std::move(residuals),
        modeGroups};
}

}  // namespace dynsyslab::core
