#include <dynsyslab/core/PCA.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace dynsyslab::core {
namespace {

void validate(const std::vector<double>& values,std::size_t samples,
              std::size_t features,const char* method) {
    if (!samples || !features || values.size()!=samples*features)
        throw std::invalid_argument(std::string(method)+" requires a nonempty samples-by-features matrix");
    if (!std::all_of(values.begin(),values.end(),[](double value) { return std::isfinite(value); }))
        throw std::invalid_argument(std::string(method)+" requires finite values");
}

double frobenius(const std::vector<double>& values) {
    double sum=0.0;
    for (const double value:values) sum+=value*value;
    return std::sqrt(sum);
}

std::vector<double> columnMean(const std::vector<double>& values,
                               std::size_t samples,std::size_t features) {
    std::vector<double> result(features);
    for (std::size_t row=0;row<samples;++row)
        for (std::size_t column=0;column<features;++column)
            result[column]+=values[row*features+column];
    for (double& value:result) value/=static_cast<double>(samples);
    return result;
}

std::vector<double> subtractMean(const std::vector<double>& values,
                                 const std::vector<double>& mean,
                                 std::size_t samples,std::size_t features) {
    std::vector<double> result(values);
    for (std::size_t row=0;row<samples;++row)
        for (std::size_t column=0;column<features;++column)
            result[row*features+column]-=mean[column];
    return result;
}

std::vector<double> singularValueThreshold(const std::vector<double>& values,
    std::size_t samples,std::size_t features,double threshold) {
    const auto decomposition=singularValueDecomposition(values,samples,features);
    std::size_t keep=0;
    while (keep<decomposition.rank
        && decomposition.singularValues[keep]>threshold) ++keep;
    std::vector<double> result(samples*features);
    for (std::size_t row=0;row<samples;++row)
        for (std::size_t column=0;column<features;++column)
            for (std::size_t mode=0;mode<keep;++mode)
                result[row*features+column]+=decomposition.left[row*decomposition.rank+mode]
                    *(decomposition.singularValues[mode]-threshold)
                    *decomposition.right[mode*features+column];
    return result;
}

}  // namespace

PCAResult principalComponentAnalysis(const std::vector<double>& values,
    std::size_t samples,std::size_t features,const PCAOptions& options) {
    validate(values,samples,features,"PCA");
    PCAResult result; result.samples=samples; result.features=features;
    result.mean=options.center?columnMean(values,samples,features)
                              :std::vector<double>(features);
    const auto centered=subtractMean(values,result.mean,samples,features);
    const auto full=singularValueDecomposition(centered,samples,features);
    result.rank=options.rank?std::min(options.rank,full.rank):full.rank;
    result.singularValues.assign(full.singularValues.begin(),
        full.singularValues.begin()+static_cast<std::ptrdiff_t>(result.rank));
    result.components.resize(result.rank*features);
    std::copy_n(full.right.begin(),result.components.size(),result.components.begin());
    result.scores.resize(samples*result.rank);
    for (std::size_t row=0;row<samples;++row)
        for (std::size_t mode=0;mode<result.rank;++mode)
            result.scores[row*result.rank+mode]=full.left[row*full.rank+mode]
                *full.singularValues[mode];
    double total=0.0;
    for (const double value:full.singularValues) total+=value*value;
    result.explainedVarianceRatio.resize(result.rank);
    if (total>0.0) for (std::size_t mode=0;mode<result.rank;++mode)
        result.explainedVarianceRatio[mode]=result.singularValues[mode]
            *result.singularValues[mode]/total;
    return result;
}

BaggedPCAResult baggedPrincipalComponentAnalysis(
    const std::vector<double>& values,std::size_t samples,std::size_t features,
    const BaggedPCAOptions& options) {
    validate(values,samples,features,"bagged PCA");
    if (!options.bags) throw std::invalid_argument("bagged PCA requires at least one bag");
    const std::uint64_t seed=options.seed.value_or([] {
        std::random_device source;
        return (static_cast<std::uint64_t>(source())<<32U)^source();
    }());
    PCAOptions baseOptions{options.rank,options.center};
    const auto reference=principalComponentAnalysis(
        values,samples,features,baseOptions);
    std::vector<double> components(reference.components.size());
    std::vector<double> ratios(reference.rank);
    std::vector<double> sample(values.size());
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<std::size_t> select(0,samples-1);
    for (std::size_t bag=0;bag<options.bags;++bag) {
        for (std::size_t row=0;row<samples;++row) {
            const std::size_t origin=select(generator);
            std::copy_n(values.begin()+static_cast<std::ptrdiff_t>(origin*features),
                        features,sample.begin()+static_cast<std::ptrdiff_t>(row*features));
        }
        auto trial=principalComponentAnalysis(sample,samples,features,
                                              {reference.rank,options.center});
        for (std::size_t mode=0;mode<reference.rank;++mode) {
            double alignment=0.0;
            for (std::size_t column=0;column<features;++column)
                alignment+=trial.components[mode*features+column]
                    *reference.components[mode*features+column];
            const double sign=alignment<0.0?-1.0:1.0;
            for (std::size_t column=0;column<features;++column)
                components[mode*features+column]+=sign
                    *trial.components[mode*features+column];
            ratios[mode]+=trial.explainedVarianceRatio[mode];
        }
    }
    for (std::size_t mode=0;mode<reference.rank;++mode) {
        double norm=0.0;
        for (std::size_t column=0;column<features;++column)
            norm+=components[mode*features+column]*components[mode*features+column];
        norm=std::sqrt(norm);
        if (!(norm>std::numeric_limits<double>::min()))
            throw std::runtime_error("bagged PCA produced a degenerate mean component");
        for (std::size_t column=0;column<features;++column)
            components[mode*features+column]/=norm;
        ratios[mode]/=static_cast<double>(options.bags);
    }
    PCAResult aggregate=reference; aggregate.components=std::move(components);
    aggregate.explainedVarianceRatio=std::move(ratios);
    const auto centered=subtractMean(values,aggregate.mean,samples,features);
    for (std::size_t row=0;row<samples;++row)
        for (std::size_t mode=0;mode<aggregate.rank;++mode) {
            aggregate.scores[row*aggregate.rank+mode]=0.0;
            for (std::size_t column=0;column<features;++column)
                aggregate.scores[row*aggregate.rank+mode]+=centered[row*features+column]
                    *aggregate.components[mode*features+column];
        }
    return {std::move(aggregate),options.bags,seed};
}

RobustPCAResult robustPrincipalComponentAnalysis(
    const std::vector<double>& values,std::size_t samples,std::size_t features,
    const RobustPCAOptions& options) {
    validate(values,samples,features,"robust PCA");
    if (options.lambda<0.0 || !std::isfinite(options.lambda)
        || !(options.tolerance>0.0) || !std::isfinite(options.tolerance)
        || !options.maximumIterations)
        throw std::invalid_argument("invalid robust PCA options");
    RobustPCAResult result; result.samples=samples; result.features=features;
    result.lambda=options.lambda>0.0?options.lambda
        :1.0/std::sqrt(static_cast<double>(std::max(samples,features)));
    result.lowRank.assign(values.size(),0.0); result.sparse.assign(values.size(),0.0);
    std::vector<double> dual(values.size()),argument(values.size()),residual(values.size());
    const double rawNorm=frobenius(values);
    const double norm=rawNorm>std::numeric_limits<double>::min()?rawNorm:1.0;
    const auto initial=singularValueDecomposition(values,samples,features,1);
    const double spectral=initial.singularValues.empty()?0.0:initial.singularValues.front();
    const double mu=spectral>std::numeric_limits<double>::min()
        ?1.25/spectral:1.25;
    for (std::size_t iteration=1;iteration<=options.maximumIterations;++iteration) {
        for (std::size_t index=0;index<values.size();++index)
            argument[index]=values[index]-result.sparse[index]+dual[index]/mu;
        result.lowRank=singularValueThreshold(argument,samples,features,1.0/mu);
        for (std::size_t index=0;index<values.size();++index) {
            const double value=values[index]-result.lowRank[index]+dual[index]/mu;
            result.sparse[index]=std::copysign(
                std::max(std::abs(value)-result.lambda/mu,0.0),value);
            residual[index]=values[index]-result.lowRank[index]-result.sparse[index];
            dual[index]+=mu*residual[index];
        }
        result.iterations=iteration; result.relativeResidual=frobenius(residual)/norm;
        if (result.relativeResidual<options.tolerance) {
            result.converged=true; break;
        }
    }
    return result;
}

HankelPODResult hankelProperOrthogonalDecomposition(
    const std::vector<double>& values,std::size_t samples,std::size_t features,
    const HankelPODOptions& options) {
    validate(values,samples,features,"Hankel POD");
    if (options.matrixSize<2 || !options.delayStride)
        throw std::invalid_argument("Hankel POD requires at least two delays and positive stride");
    const std::size_t last=(options.matrixSize-1)*options.delayStride;
    if (samples<=last+2)
        throw std::invalid_argument("Hankel POD delay settings leave fewer than three samples");
    HankelPODResult result; result.embeddedSamples=samples-last;
    result.embeddedFeatures=options.matrixSize*features;
    result.matrixSize=options.matrixSize; result.delayStride=options.delayStride;
    std::vector<double> hankel(result.embeddedSamples*result.embeddedFeatures);
    for (std::size_t row=0;row<result.embeddedSamples;++row)
        for (std::size_t delay=0;delay<options.matrixSize;++delay)
            for (std::size_t column=0;column<features;++column)
                hankel[row*result.embeddedFeatures+delay*features+column]
                    =values[(row+delay*options.delayStride)*features+column];
    result.mean=options.center?columnMean(
        hankel,result.embeddedSamples,result.embeddedFeatures)
        :std::vector<double>(result.embeddedFeatures);
    const auto centered=subtractMean(hankel,result.mean,
                                     result.embeddedSamples,result.embeddedFeatures);
    result.decomposition=singularValueDecomposition(
        centered,result.embeddedSamples,result.embeddedFeatures,options.rank);
    return result;
}

}  // namespace dynsyslab::core
