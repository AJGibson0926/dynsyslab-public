#include <dynsyslab/core/CoordinateClosure.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct PairCandidate {
    std::vector<double> left;
    std::vector<double> right;
    double score{};
};

double normalizedSeparation(const std::vector<double>& left,
    const std::vector<double>& right, const ObservableClosureOptions& options) {
    double squared{};
    for (std::size_t axis=0;axis<left.size();++axis) {
        const double range=options.stateMaximum[axis]-options.stateMinimum[axis];
        const double value=(left[axis]-right[axis])/range;
        squared+=value*value;
    }
    return std::sqrt(squared/double(left.size()));
}

std::vector<double> normalizedState(const std::vector<double>& state,
    const ObservableClosureOptions& options) {
    std::vector<double> result(state.size());
    for(std::size_t axis=0;axis<state.size();++axis)
        result[axis]=(state[axis]-options.stateMinimum[axis])
            /(options.stateMaximum[axis]-options.stateMinimum[axis]);
    return result;
}

std::vector<double> physicalState(const double* normalized,std::size_t dimension,
    const ObservableClosureOptions& options) {
    std::vector<double> result(dimension);
    for(std::size_t axis=0;axis<dimension;++axis)
        result[axis]=options.stateMinimum[axis]+normalized[axis]
            *(options.stateMaximum[axis]-options.stateMinimum[axis]);
    return result;
}

double outputResidual(const std::vector<double>& left,
    const std::vector<double>& right,const std::vector<double>& scale) {
    double result{};
    for(std::size_t index=0;index<left.size();++index)
        result=std::max(result,std::abs(left[index]-right[index])/scale[index]);
    return result;
}

std::vector<double> coordinateScale(const Observable& jet,
    const std::vector<std::vector<double>>& samples) {
    const std::size_t outputs=jet.outputDimension();
    std::vector<double> low(outputs,std::numeric_limits<double>::infinity());
    std::vector<double> high(outputs,-std::numeric_limits<double>::infinity());
    for(const auto& state:samples) {
        const auto value=jet.evaluateReal(state.data());
        for(std::size_t output=0;output<outputs;++output) {
            low[output]=std::min(low[output],value[output]);
            high[output]=std::max(high[output],value[output]);
        }
    }
    std::vector<double> scale(outputs);
    for(std::size_t output=0;output<outputs;++output)
        scale[output]=std::max({1.0e-12,high[output]-low[output],
            std::abs(low[output]),std::abs(high[output])});
    return scale;
}

double halton(std::uint64_t index,std::size_t base) {
    double fraction=1.0,value=0.0;
    while(index) { fraction/=double(base);value+=fraction*double(index%base);index/=base; }
    return value;
}

std::vector<std::size_t> primes(std::size_t count) {
    std::vector<std::size_t> result;
    for(std::size_t candidate=2;result.size()<count;++candidate) {
        bool prime=true;
        for(const auto divisor:result) {
            if(divisor*divisor>candidate)break;
            if(candidate%divisor==0){prime=false;break;}
        }
        if(prime)result.push_back(candidate);
    }
    return result;
}

bool inDomain(const std::vector<double>& state,const ObservableClosureOptions& options) {
    for(std::size_t axis=0;axis<state.size();++axis)
        if(state[axis]<options.stateMinimum[axis]||state[axis]>options.stateMaximum[axis])return false;
    return true;
}

std::vector<std::vector<double>> sampleStates(const ODESystem& system,
    const Observable& maximumJet,const ObservableClosureOptions& options) {
    const auto dimension=system.descriptor().stateDimension;
    std::vector<std::vector<double>> samples;
    auto append=[&](std::vector<double> state) {
        if(!inDomain(state,options)||system.invalid(0.0,state.data()))return;
        try { (void)maximumJet.evaluateReal(state.data()); }
        catch(const std::exception&) { return; }
        for(const auto& existing:samples)if(existing==state)return;
        samples.push_back(std::move(state));
    };
    for(std::size_t offset=0;offset<options.pairSeeds.size();offset+=dimension) {
        append({options.pairSeeds.begin()+static_cast<std::ptrdiff_t>(offset),
            options.pairSeeds.begin()+static_cast<std::ptrdiff_t>(offset+dimension)});
    }
    std::vector<double> center(dimension);
    for(std::size_t axis=0;axis<dimension;++axis)
        center[axis]=0.5*(options.stateMinimum[axis]+options.stateMaximum[axis]);
    append(center);
    const auto bases=primes(dimension);
    const std::uint64_t offset=1+(options.seed%1000003ULL);
    for(std::uint64_t index=0;samples.size()<options.sampleCount&&index<50*options.sampleCount;++index) {
        std::vector<double> state(dimension);
        for(std::size_t axis=0;axis<dimension;++axis) {
            const double unit=halton(index+offset,bases[axis]);
            state[axis]=options.stateMinimum[axis]
                +unit*(options.stateMaximum[axis]-options.stateMinimum[axis]);
        }
        append(std::move(state));
    }
    if(samples.size()<std::min<std::size_t>(options.sampleCount,dimension+1))
        throw std::runtime_error("closure search could not sample enough valid states");
    return samples;
}

std::vector<PairCandidate> candidatePairs(const Observable& jet,
    const std::vector<std::vector<double>>& samples,const std::vector<double>& scale,
    const ObservableClosureOptions& options) {
    std::vector<std::vector<double>> values;
    values.reserve(samples.size());
    for(const auto& state:samples)values.push_back(jet.evaluateReal(state.data()));
    std::vector<PairCandidate> candidates;
    const auto limit=options.pairCandidateCount;
    auto consider=[&](const std::vector<double>& left,const std::vector<double>& right,
        const std::vector<double>& leftValue,const std::vector<double>& rightValue) {
        const double separation=normalizedSeparation(left,right,options);
        if(separation<options.minimumNormalizedStateSeparation)return;
        const double score=outputResidual(leftValue,rightValue,scale)/separation;
        candidates.push_back({left,right,score});
    };
    const auto dimension=jet.stateDimension();
    for(std::size_t offset=0;offset<options.pairSeeds.size();offset+=2*dimension) {
        std::vector<double> left(options.pairSeeds.begin()+static_cast<std::ptrdiff_t>(offset),
            options.pairSeeds.begin()+static_cast<std::ptrdiff_t>(offset+dimension));
        std::vector<double> right(options.pairSeeds.begin()+static_cast<std::ptrdiff_t>(offset+dimension),
            options.pairSeeds.begin()+static_cast<std::ptrdiff_t>(offset+2*dimension));
        consider(left,right,jet.evaluateReal(left.data()),jet.evaluateReal(right.data()));
    }
    for(std::size_t left=0;left<samples.size();++left)
        for(std::size_t right=left+1;right<samples.size();++right)
            consider(samples[left],samples[right],values[left],values[right]);
    std::stable_sort(candidates.begin(),candidates.end(),
        [](const PairCandidate& left,const PairCandidate& right){return left.score<right.score;});
    if(candidates.size()>limit)candidates.resize(limit);
    return candidates;
}

struct RefinementEvaluation {
    std::vector<double> residual;
    std::vector<double> jacobian;
    std::vector<double> leftState,rightState,leftValue,rightValue;
    double objective{},separation{},outputError{};
};

RefinementEvaluation evaluatePair(const ODESystem& system,const Observable& jet,
    const std::vector<double>& normalized,
    const std::vector<double>& scale,const ObservableClosureOptions& options,bool withJacobian) {
    const auto dimension=jet.stateDimension(),outputs=jet.outputDimension();
    RefinementEvaluation result;
    result.leftState=physicalState(normalized.data(),dimension,options);
    result.rightState=physicalState(normalized.data()+dimension,dimension,options);
    if(system.invalid(0.0,result.leftState.data())||system.invalid(0.0,result.rightState.data()))
        throw std::domain_error("closure refinement entered an invalid system state");
    result.leftValue=jet.evaluateReal(result.leftState.data());
    result.rightValue=jet.evaluateReal(result.rightState.data());
    result.residual.resize(outputs+1);
    for(std::size_t output=0;output<outputs;++output)
        result.residual[output]=(result.leftValue[output]-result.rightValue[output])/scale[output];
    double separationSquared{};
    for(std::size_t axis=0;axis<dimension;++axis) {
        const double difference=normalized[axis]-normalized[dimension+axis];
        separationSquared+=difference*difference/double(dimension);
    }
    result.separation=std::sqrt(separationSquared);
    constexpr double barrierWeight=10.0;
    result.residual[outputs]=barrierWeight*std::max(0.0,
        options.minimumNormalizedStateSeparation-result.separation);
    result.objective=std::inner_product(result.residual.begin(),result.residual.end(),
        result.residual.begin(),0.0);
    result.outputError=outputResidual(result.leftValue,result.rightValue,scale);
    if(!withJacobian)return result;
    result.jacobian.assign((outputs+1)*2*dimension,0.0);
    const auto leftJacobian=jet.jacobian(result.leftState.data(),ObservableJacobianBackend::analytic);
    const auto rightJacobian=jet.jacobian(result.rightState.data(),ObservableJacobianBackend::analytic);
    for(std::size_t output=0;output<outputs;++output)
        for(std::size_t axis=0;axis<dimension;++axis) {
            const double range=options.stateMaximum[axis]-options.stateMinimum[axis];
            result.jacobian[output*2*dimension+axis]=leftJacobian[output*dimension+axis].real()*range/scale[output];
            result.jacobian[output*2*dimension+dimension+axis]=-rightJacobian[output*dimension+axis].real()*range/scale[output];
        }
    if(result.separation<options.minimumNormalizedStateSeparation&&result.separation>1.0e-14)
        for(std::size_t axis=0;axis<dimension;++axis) {
            const double derivative=-barrierWeight*(normalized[axis]-normalized[dimension+axis])
                /(double(dimension)*result.separation);
            result.jacobian[outputs*2*dimension+axis]=derivative;
            result.jacobian[outputs*2*dimension+dimension+axis]=-derivative;
        }
    return result;
}

ObservableCollisionWitness refinePair(const ODESystem& system,const Observable& jet,
    const PairCandidate& candidate,const std::vector<double>& scale,
    const ObservableClosureOptions& options) {
    const auto dimension=jet.stateDimension();
    auto normalized=normalizedState(candidate.left,options);
    const auto right=normalizedState(candidate.right,options);
    normalized.insert(normalized.end(),right.begin(),right.end());
    auto current=evaluatePair(system,jet,normalized,scale,options,true);
    std::size_t completed{};
    double ridge=1.0e-8;
    for(;completed<options.maximumRefinementIterations;++completed) {
        if(current.outputError<=options.collisionTolerance
            &&current.separation>=options.minimumNormalizedStateSeparation)break;
        std::vector<double> rhs(current.residual.size());
        std::transform(current.residual.begin(),current.residual.end(),rhs.begin(),
            [](double value){return -value;});
        numerics::LeastSquaresOptions leastSquares;
        leastSquares.ridge=ridge;
        const auto step=numerics::solveLeastSquares(current.jacobian,current.residual.size(),
            2*dimension,rhs,1,leastSquares);
        bool accepted=false;
        for(double fraction=1.0;fraction>=1.0/1024.0;fraction*=0.5) {
            auto trial=normalized;
            for(std::size_t index=0;index<trial.size();++index)
                trial[index]=std::clamp(trial[index]+fraction*step[index],0.0,1.0);
            try {
                auto evaluated=evaluatePair(system,jet,trial,scale,options,true);
                if(evaluated.objective<current.objective) {
                    normalized=std::move(trial);current=std::move(evaluated);
                    ridge=std::max(1.0e-14,ridge*0.25);accepted=true;break;
                }
            } catch(const std::exception&) { }
        }
        if(!accepted) { ridge*=100.0;if(ridge>1.0e8)break; }
    }
    ObservableCollisionWitness witness;
    witness.found=current.outputError<=options.collisionTolerance
        &&current.separation>=options.minimumNormalizedStateSeparation;
    if(!witness.found)return witness;
    witness.leftState=current.leftState;witness.rightState=current.rightState;
    witness.leftCoordinates=current.leftValue;witness.rightCoordinates=current.rightValue;
    witness.normalizedResidual=current.outputError;
    witness.normalizedStateSeparation=current.separation;
    witness.refinementIterations=completed;
    return witness;
}

ObservableClosureOrderReport analyzeOrder(const ODESystem& system,const Observable& jet,
    const Observable& observable,
    const std::vector<AnalyticExpression>& field,
    const std::vector<std::vector<double>>& samples,const ObservableClosureOptions& options,
    std::size_t order) {
    ObservableClosureOrderReport report;
    report.highestDerivativeOrder=order;report.coordinateCount=jet.outputDimension();
    report.minimumSampledRank=jet.stateDimension();
    report.minimumSampledSingularValue=std::numeric_limits<double>::infinity();
    const auto scale=coordinateScale(jet,samples);
    for(const auto& state:samples) {
        const auto jacobian=jet.jacobian(state.data(),ObservableJacobianBackend::analytic);
        std::vector<double> normalized(jacobian.size());
        for(std::size_t output=0;output<jet.outputDimension();++output)
            for(std::size_t axis=0;axis<jet.stateDimension();++axis)
                normalized[output*jet.stateDimension()+axis]=jacobian[output*jet.stateDimension()+axis].real()
                    *(options.stateMaximum[axis]-options.stateMinimum[axis])/scale[output];
        const auto svd=numerics::singularValueDecomposition(normalized,jet.outputDimension(),jet.stateDimension());
        const double largest=svd.singularValues.empty()?0.0:svd.singularValues.front();
        const double cutoff=options.rankTolerance*largest;
        const std::size_t rank=static_cast<std::size_t>(std::count_if(svd.singularValues.begin(),
            svd.singularValues.end(),[&](double value){return value>cutoff;}));
        const double smallest=svd.singularValues.size()<jet.stateDimension()?0.0:svd.singularValues.back();
        if(rank<report.minimumSampledRank||
            (rank==report.minimumSampledRank&&smallest<report.minimumSampledSingularValue)) {
            report.minimumSampledRank=rank;report.minimumSampledSingularValue=smallest;
            report.minimumRankState=state;
        }
    }
    const auto successor=observable.flowDerivative(field,order+1);
    const auto successorScale=coordinateScale(successor,samples);
    for(const auto& candidate:candidatePairs(jet,samples,scale,options)) {
        auto witness=refinePair(system,jet,candidate,scale,options);
        if(witness.found) {
            witness.successorDerivativeLeft=successor.evaluateReal(witness.leftState.data());
            witness.successorDerivativeRight=successor.evaluateReal(witness.rightState.data());
            witness.normalizedSuccessorDifference=outputResidual(
                witness.successorDerivativeLeft,witness.successorDerivativeRight,successorScale);
            if(!report.coordinateCollision.found)report.coordinateCollision=witness;
            if(witness.normalizedSuccessorDifference>options.successorDifferenceTolerance) {
                report.closureObstruction=std::move(witness);break;
            }
        }
    }
    return report;
}

void validate(const ODESystem& system,const Observable& observable,
    const ObservableClosureOptions& options) {
    const auto dimension=system.descriptor().stateDimension;
    if(observable.stateDimension()!=dimension||observable.expressions().empty())
        throw std::invalid_argument("closure checks require an analytic system observable");
    if(options.stateMinimum.size()!=dimension||options.stateMaximum.size()!=dimension)
        throw std::invalid_argument("closure domain must bound every state component");
    for(std::size_t axis=0;axis<dimension;++axis)
        if(!std::isfinite(options.stateMinimum[axis])||!std::isfinite(options.stateMaximum[axis])
            ||options.stateMinimum[axis]>=options.stateMaximum[axis])
            throw std::invalid_argument("closure domain bounds must be finite and increasing");
    if(options.highestDerivativeOrder>8||options.sampleCount<2||options.sampleCount>4096
        ||!options.pairCandidateCount||options.pairCandidateCount>4096
        ||!options.maximumRefinementIterations||options.maximumRefinementIterations>1000)
        throw std::invalid_argument("closure search limits are outside supported bounds");
    if(!(options.minimumNormalizedStateSeparation>0.0)
        ||!(options.minimumNormalizedStateSeparation<1.0)
        ||!std::isfinite(options.minimumNormalizedStateSeparation)
        ||!(options.collisionTolerance>0.0)||!std::isfinite(options.collisionTolerance)
        ||!(options.successorDifferenceTolerance>0.0)
        ||!std::isfinite(options.successorDifferenceTolerance)
        ||!(options.rankTolerance>0.0)||!std::isfinite(options.rankTolerance))
        throw std::invalid_argument("closure tolerances must be positive and finite");
    if(options.pairSeeds.size()%(2*dimension)!=0)
        throw std::invalid_argument("closure pair seeds must contain complete state pairs");
    if(system.driftExpressions().size()!=dimension)
        throw std::invalid_argument("closure checks require an analytic autonomous drift field");
}

}  // namespace

Observable observableJet(const Observable& observable,
    const std::vector<AnalyticExpression>& autonomousField,
    std::size_t highestDerivativeOrder) {
    if(observable.expressions().empty())
        throw std::invalid_argument("observable jet requires analytic expressions");
    std::vector<Observable> members{observable};
    for(std::size_t order=1;order<=highestDerivativeOrder;++order)
        members.push_back(observable.flowDerivative(autonomousField,order));
    return Observable::stack(observable.descriptor().name+"_jet_through_"
        +std::to_string(highestDerivativeOrder),"Observable jet through derivative order "
        +std::to_string(highestDerivativeOrder),members);
}

ObservableClosureReport analyzeObservableClosure(const ODESystem& system,
    const Observable& observable,const ObservableClosureOptions& options) {
    validate(system,observable,options);
    const auto field=system.driftExpressions();
    std::vector<Observable> jets;
    for(std::size_t order=0;order<=options.highestDerivativeOrder;++order)
        jets.push_back(observableJet(observable,field,order));
    const auto samples=sampleStates(system,jets.back(),options);
    ObservableClosureReport report;
    report.stateDimension=system.descriptor().stateDimension;
    report.observableDimension=observable.outputDimension();
    report.requestedSamples=options.sampleCount;report.acceptedSamples=samples.size();
    for(std::size_t order=0;order<=options.highestDerivativeOrder;++order)
        report.orders.push_back(analyzeOrder(
            system,jets[order],observable,field,samples,options,order));
    return report;
}

}  // namespace dynsyslab::core
