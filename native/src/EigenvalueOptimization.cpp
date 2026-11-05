#include <dynsyslab/core/EigenvalueOptimization.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Evaluation {double objective{};std::vector<double> gradient;};
struct MinimizeResult {std::vector<double> parameters;Evaluation evaluation;
    std::size_t iterations{},evaluations{};bool converged{};std::string message;};

double dot(const std::vector<double>& left,const std::vector<double>& right){
    return std::inner_product(left.begin(),left.end(),right.begin(),0.0);}
double normInfinity(const std::vector<double>& value){double result=0.0;
    for(const double item:value)result=std::max(result,std::abs(item));return result;}
void clamp(std::vector<double>& value,double lower,double upper){
    for(auto& item:value)item=std::clamp(item,lower,upper);}

Evaluation objective(const std::vector<double>& parameters,
                     const std::vector<std::complex<double>>& measurements,
                     std::size_t samples,std::size_t trajectories){
    const std::size_t modes=parameters.size()/2;
    std::vector<std::complex<double>> eigenvalues(modes);
    for(std::size_t mode=0;mode<modes;++mode)eigenvalues[mode]=
        {parameters[mode],parameters[modes+mode]};
    std::vector<std::complex<double>> powers(samples*modes);
    for(std::size_t sample=0;sample<samples;++sample)
        for(std::size_t mode=0;mode<modes;++mode)
            powers[sample*modes+mode]=
                std::pow(eigenvalues[mode],static_cast<double>(sample));
    numerics::LeastSquaresOptions options;
    options.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    const auto boundary=numerics::solveLeastSquares(
        powers,samples,modes,measurements,trajectories,options);
    const auto prediction=numerics::multiply(
        powers,samples,modes,boundary,trajectories);
    std::vector<std::complex<double>> residual(measurements.size());
    double normalizer=0.0,residualSquared=0.0;
    for(std::size_t sample=0;sample<samples;++sample)
        for(std::size_t trajectory=0;trajectory<trajectories;++trajectory) {
            auto& value=residual[sample*trajectories+trajectory];
            value=measurements[sample*trajectories+trajectory]
                -prediction[sample*trajectories+trajectory];
            residualSquared+=std::norm(value);
            normalizer+=std::norm(measurements[sample*trajectories+trajectory]);
        }
    normalizer=std::max(normalizer,std::numeric_limits<double>::min());
    Evaluation result;result.objective=residualSquared/normalizer;
    result.gradient.assign(parameters.size(),0.0);
    for(std::size_t mode=0;mode<modes;++mode){
        std::complex<double> projection{};
        for(std::size_t sample=1;sample<samples;++sample){
            const auto derivative=static_cast<double>(sample)*std::pow(
                eigenvalues[mode],static_cast<double>(sample-1));
            for(std::size_t trajectory=0;trajectory<trajectories;++trajectory)
                projection+=std::conj(residual[sample*trajectories+trajectory])
                    *derivative*boundary[mode*trajectories+trajectory];
        }
        result.gradient[mode]=-2.0*projection.real()/normalizer;
        result.gradient[modes+mode]=2.0*projection.imag()/normalizer;
    }
    return result;
}

std::vector<double> projectedGradient(const std::vector<double>& parameters,
    const std::vector<double>& gradient,double lower,double upper){
    auto result=gradient;const double edge=1.0e-13*std::max(1.0,std::abs(upper-lower));
    for(std::size_t index=0;index<result.size();++index){
        if((parameters[index]<=lower+edge&&result[index]>0.0)
            ||(parameters[index]>=upper-edge&&result[index]<0.0))result[index]=0.0;}
    return result;
}

MinimizeResult minimize(std::vector<double> parameters,
                        const std::vector<std::complex<double>>& measurements,
                        std::size_t samples,std::size_t trajectories,
                        const EigenvalueOptimizationOptions& options){
    clamp(parameters,options.lowerBound,options.upperBound);
    Evaluation current=objective(parameters,measurements,samples,trajectories);
    std::size_t evaluations=1;
    std::vector<std::vector<double>> historyS,historyY;std::vector<double> rho;
    for(std::size_t iteration=0;iteration<options.maximumIterations;++iteration){
        const auto projected=projectedGradient(parameters,current.gradient,
            options.lowerBound,options.upperBound);
        if(normInfinity(projected)<=options.tolerance)
            return {parameters,current,iteration,evaluations,true,"projected gradient tolerance"};
        std::vector<double> direction=projected,alpha(historyS.size());
        for(std::size_t reverse=historyS.size();reverse>0;--reverse){
            const std::size_t index=reverse-1;alpha[index]=rho[index]*dot(historyS[index],direction);
            for(std::size_t component=0;component<direction.size();++component)
                direction[component]-=alpha[index]*historyY[index][component];}
        if(!historyS.empty()){
            const double yy=dot(historyY.back(),historyY.back());
            const double scale=yy>0.0?dot(historyS.back(),historyY.back())/yy:1.0;
            for(auto& item:direction)item*=scale;}
        for(std::size_t index=0;index<historyS.size();++index){
            const double beta=rho[index]*dot(historyY[index],direction);
            for(std::size_t component=0;component<direction.size();++component)
                direction[component]+=historyS[index][component]*(alpha[index]-beta);}
        for(auto& item:direction)item=-item;
        if(dot(direction,projected)>=0.0)
            for(std::size_t index=0;index<direction.size();++index)direction[index]=-projected[index];
        double step=1.0;std::vector<double> candidate(parameters.size());Evaluation trial;
        bool accepted=false;const double previous=current.objective;
        for(std::size_t line=0;line<40;++line){
            for(std::size_t index=0;index<candidate.size();++index)
                candidate[index]=std::clamp(parameters[index]+step*direction[index],
                    options.lowerBound,options.upperBound);
            std::vector<double> displacement(candidate.size());
            for(std::size_t index=0;index<candidate.size();++index)
                displacement[index]=candidate[index]-parameters[index];
            if(normInfinity(displacement)<=std::numeric_limits<double>::epsilon()){
                step*=0.5;continue;}
            trial=objective(candidate,measurements,samples,trajectories);++evaluations;
            if(trial.objective<=current.objective+1.0e-4*dot(current.gradient,displacement)){
                accepted=true;break;}step*=0.5;
        }
        if(!accepted)return {parameters,current,iteration,evaluations,false,"line search failed"};
        std::vector<double> s(parameters.size()),y(parameters.size());
        for(std::size_t index=0;index<s.size();++index){
            s[index]=candidate[index]-parameters[index];
            y[index]=trial.gradient[index]-current.gradient[index];}
        const double curvature=dot(s,y);
        if(curvature>1.0e-14*std::max(1.0,dot(s,s))){
            if(historyS.size()==options.memory){historyS.erase(historyS.begin());
                historyY.erase(historyY.begin());rho.erase(rho.begin());}
            historyS.push_back(s);historyY.push_back(y);rho.push_back(1.0/curvature);}
        parameters=std::move(candidate);current=std::move(trial);
        if(std::abs(previous-current.objective)<=options.tolerance
            *std::max(1.0,std::abs(previous)))
            return {parameters,current,iteration+1,evaluations,true,"objective tolerance"};
    }
    return {parameters,current,options.maximumIterations,evaluations,false,
        "maximum iterations reached"};
}

void symmetrize(std::vector<std::complex<double>>& values,
                const std::vector<std::complex<double>>& initial){
    std::vector<bool> used(values.size());
    for(std::size_t index=0;index<initial.size();++index){
        if(used[index]||std::abs(initial[index].imag())<=1.0e-12)continue;
        for(std::size_t partner=index+1;partner<initial.size();++partner)
            if(!used[partner]&&std::abs(initial[partner]-std::conj(initial[index]))
                <=1.0e-8*std::max(1.0,std::abs(initial[index]))){
                const auto representative=0.5*(values[index]+std::conj(values[partner]));
                values[index]=representative;values[partner]=std::conj(representative);
                used[index]=used[partner]=true;break;}}
}

}  // namespace

EigenvalueOptimizationResult optimizeTrajectoryEigenvalues(
    const KoopmanTrajectoryEnsemble& ensemble,const Observable& observable,
    const EigenvalueOptimizationOptions& options){
    if(ensemble.samples<3||!ensemble.trajectories||!ensemble.stateDimension
        ||ensemble.time.size()!=ensemble.samples
        ||ensemble.states.size()!=ensemble.samples*ensemble.trajectories*ensemble.stateDimension
        ||observable.stateDimension()!=ensemble.stateDimension
        ||options.initialDiscreteGroups.size()!=observable.outputDimension()
        ||!(options.lowerBound<options.upperBound)||!options.maximumIterations
        ||!options.memory||!(options.tolerance>0.0))
        throw std::invalid_argument("eigenvalue optimization inputs are inconsistent");
    const double dt=ensemble.time[1]-ensemble.time[0];
    if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument(
        "eigenvalue optimization requires increasing finite time");
    for(std::size_t sample=1;sample<ensemble.samples;++sample)
        if(std::abs((ensemble.time[sample]-ensemble.time[sample-1])-dt)
            >1.0e-9*std::max(1.0,dt))throw std::invalid_argument(
                "eigenvalue optimization requires a uniform shared time grid");
    EigenvalueOptimizationResult result;result.deltaTime=dt;std::size_t offset=0;
    for(std::size_t output=0;output<observable.outputDimension();++output){
        const auto& initial=options.initialDiscreteGroups[output];
        if(initial.empty())throw std::invalid_argument("initial eigenvalue groups cannot be empty");
        std::vector<std::complex<double>> measurements(
            ensemble.samples*ensemble.trajectories);
        for(std::size_t trajectory=0;trajectory<ensemble.trajectories;++trajectory)
            for(std::size_t sample=0;sample<ensemble.samples;++sample){
                const auto value=observable.evaluate(ensemble.states.data()+
                    (trajectory*ensemble.samples+sample)*ensemble.stateDimension);
                measurements[sample*ensemble.trajectories+trajectory]=value[output];}
        std::vector<double> parameters(initial.size()*2);
        for(std::size_t mode=0;mode<initial.size();++mode){parameters[mode]=initial[mode].real();
            parameters[initial.size()+mode]=initial[mode].imag();}
        auto optimized=minimize(std::move(parameters),measurements,
            ensemble.samples,ensemble.trajectories,options);
        std::vector<std::complex<double>> values(initial.size());
        for(std::size_t mode=0;mode<values.size();++mode)values[mode]=
            {optimized.parameters[mode],optimized.parameters[values.size()+mode]};
        if(options.preserveConjugatePairs)symmetrize(values,initial);
        std::vector<std::size_t> group;
        for(const auto value:values){
            if(!std::isfinite(value.real())||!std::isfinite(value.imag())
                ||std::abs(value)<=std::numeric_limits<double>::min())
                throw std::runtime_error(
                    "optimized discrete eigenvalue cannot be converted to generator form");
            result.discreteEigenvalues.push_back(value);
            result.continuousEigenvalues.push_back(std::log(value)/dt);
            group.push_back(offset++);}
        const double relative=std::sqrt(std::max(0.0,optimized.evaluation.objective));
        result.modeGroups.push_back(std::move(group));
        result.diagnostics.push_back({output,values.size(),optimized.iterations,
            optimized.evaluations,optimized.evaluation.objective,relative,
            optimized.converged,std::move(optimized.message)});
    }
    return result;
}

}  // namespace dynsyslab::core
