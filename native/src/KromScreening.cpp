#include <dynsyslab/core/KromScreening.hpp>

#include <dynsyslab/core/Spectral.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

double norm(const std::vector<double>& values) {
    double sum=0.0;for(const double value:values)sum+=value*value;
    return std::sqrt(sum);
}

std::vector<double> multiply(const std::vector<double>& matrix,
    const std::vector<double>& right,std::size_t rows,std::size_t columns) {
    const std::size_t outputs=right.size()/columns;
    std::vector<double> result(rows*outputs);
    for(std::size_t row=0;row<rows;++row)for(std::size_t output=0;output<outputs;++output)
        for(std::size_t column=0;column<columns;++column)
            result[row*outputs+output]+=matrix[row*columns+column]
                *right[column*outputs+output];
    return result;
}

KromCandidateScore score(const KoopmanModel& parent,
    const std::vector<std::size_t>& indices,const std::vector<double>& initial,
    const std::vector<double>& target) {
    const auto model=parent.extract(indices);const std::size_t n=model.splitDimension();
    const auto a=model.generatorSplit();const std::size_t inputs=model.inputDimension();
    std::vector<double> middle(initial.size());for(std::size_t i=0;i<middle.size();++i)
        middle[i]=0.5*(initial[i]+target[i]);
    const std::vector<std::vector<double>> probes{initial,middle,target};
    double worst=std::numeric_limits<double>::infinity();
    double authority=std::numeric_limits<double>::infinity();std::size_t rank=n;
    for(const auto& probe:probes){const auto b=model.controlMapSplit(probe.data(),"exact");
        authority=std::min(authority,norm(b));std::vector<double> matrix(n*n*inputs);
        auto block=b;
        for(std::size_t power=0;power<n;++power){for(std::size_t row=0;row<n;++row)
            for(std::size_t input=0;input<inputs;++input)
                matrix[row*(n*inputs)+power*inputs+input]=block[row*inputs+input];
            block=multiply(a,block,n,n);}
        const auto svd=singularValueDecomposition(matrix,n,n*inputs);
        const double largest=svd.singularValues.empty()?0.0:svd.singularValues.front();
        const double threshold=std::max(1.0,largest)*1.0e-10;std::size_t localRank=0;
        double smallest=largest;
        for(const double value:svd.singularValues)if(value>threshold){++localRank;smallest=value;}
        rank=std::min(rank,localRank);worst=std::min(worst,
            largest>std::numeric_limits<double>::min()?smallest/largest:0.0);}
    const auto psi0=model.evaluateSplit(initial.data()),psi1=model.evaluateSplit(target.data());
    std::vector<double> delta(psi0.size());for(std::size_t i=0;i<delta.size();++i)
        delta[i]=psi0[i]-psi1[i];
    double low=std::numeric_limits<double>::infinity(),high=0.0;
    for(const auto value:model.continuousEigenvalues()){
        const double frequency=std::abs(value.imag());low=std::min(low,frequency);
        high=std::max(high,frequency);}
    const double spread=model.liftedDimension()>1?high-low:0.0;
    const double displacement=norm(delta);const double scoreValue=displacement
        *std::pow(std::max(worst,std::numeric_limits<double>::epsilon()),0.25)
        *std::log1p(authority)*(1.0+0.05*spread)
        /std::pow(static_cast<double>(indices.size()),0.15);
    return {indices,scoreValue,displacement,worst,authority,spread,rank,n};
}

void combinations(const std::vector<std::size_t>& seeds,std::size_t count,
    std::size_t offset,std::vector<std::size_t>& active,
    std::set<std::vector<std::size_t>>& output,std::size_t limit) {
    if(output.size()>=limit)return;if(active.size()==count){output.insert(active);return;}
    for(std::size_t i=offset;i<seeds.size()&&output.size()<limit;++i){active.push_back(seeds[i]);
        combinations(seeds,count,i+1,active,output,limit);active.pop_back();}
}

}  // namespace

std::vector<KromCandidateScore> rankKromControlCandidates(
    const KoopmanModel& parent,const std::vector<double>& initial,
    const std::vector<double>& target,const KromScreeningOptions& options) {
    if(initial.size()!=parent.stateDimension()||target.size()!=parent.stateDimension())
        throw std::invalid_argument("KROM screening state dimension mismatch");
    if(!options.maximumModeCount||options.maximumModeCount>10
        ||!options.maximumVirtualCandidates||!options.retainedCandidates)
        throw std::invalid_argument("invalid KROM screening limits");
    const std::size_t modes=parent.liftedDimension();if(!modes)
        throw std::invalid_argument("KROM screening requires modes");
    std::vector<KromCandidateScore> singles;singles.reserve(modes);
    for(std::size_t mode=0;mode<modes;++mode)singles.push_back(
        score(parent,{mode},initial,target));
    std::sort(singles.begin(),singles.end(),[](const auto& left,const auto& right){
        return left.selectionScore>right.selectionScore;});
    std::vector<std::size_t> seeds;for(std::size_t i=0;i<std::min({modes,
        options.seedModeCount,options.maximumVirtualCandidates});++i)
        seeds.push_back(singles[i].modeIndices.front());
    std::sort(seeds.begin(),seeds.end());std::set<std::vector<std::size_t>> selections;
    for(std::size_t mode=0;mode<modes&&selections.size()<options.maximumVirtualCandidates;++mode)
        selections.insert({mode});
    std::vector<std::size_t> active;
    for(std::size_t size=2;size<=std::min(options.maximumModeCount,seeds.size());++size)
        combinations(seeds,size,0,active,selections,options.maximumVirtualCandidates);
    std::vector<std::size_t> frequencyOrder(modes);std::iota(frequencyOrder.begin(),frequencyOrder.end(),0);
    std::sort(frequencyOrder.begin(),frequencyOrder.end(),[&](std::size_t left,std::size_t right){
        return std::abs(parent.continuousEigenvalues()[left].imag())
            <std::abs(parent.continuousEigenvalues()[right].imag());});
    for(std::size_t size=2;size<=std::min(options.maximumModeCount,modes);++size)
        for(std::size_t begin=0;begin+size<=modes&&selections.size()<options.maximumVirtualCandidates;++begin){
            std::vector<std::size_t> value(frequencyOrder.begin()+begin,frequencyOrder.begin()+begin+size);
            std::sort(value.begin(),value.end());selections.insert(std::move(value));}
    std::vector<KromCandidateScore> result;result.reserve(selections.size());
    for(const auto& selection:selections)result.push_back(score(parent,selection,initial,target));
    std::sort(result.begin(),result.end(),[](const auto& left,const auto& right){
        if(left.selectionScore!=right.selectionScore)return left.selectionScore>right.selectionScore;
        if(left.modeIndices.size()!=right.modeIndices.size())return left.modeIndices.size()<right.modeIndices.size();
        return left.modeIndices<right.modeIndices;});
    if(result.size()>options.retainedCandidates)result.resize(options.retainedCandidates);
    return result;
}

}  // namespace dynsyslab::core
