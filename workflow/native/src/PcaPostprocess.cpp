#include "PostprocessInternal.hpp"

#include <dynsyslab/core/PCA.hpp>

#include <cstdint>
#include <cmath>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

void addPcaArrays(dynsyslab::core::NpzRecord& archive,
                  const dynsyslab::core::PCAResult& result) {
    archive.arrays["mean"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.features},result.mean);
    archive.arrays["components"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.rank,result.features},result.components);
    archive.arrays["singular_values"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.rank},result.singularValues);
    archive.arrays["scores"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.samples,result.rank},result.scores);
    archive.arrays["explained_variance_ratio"]=
        dynsyslab::core::NpyArray::fromDoubles(
            {result.rank},result.explainedVarianceRatio);
}

std::string numeric(double value) {
    std::ostringstream out; out<<std::setprecision(17)<<value; return out.str();
}

std::size_t count(const json::Value& value,std::size_t fallback,
                  const char* name,bool allowZero=false) {
    const double number=value.number(static_cast<double>(fallback));
    if (number<0.0 || (!allowZero && number==0.0)
        || number>9007199254740991.0 || std::floor(number)!=number)
        throw std::invalid_argument(std::string(name)+" must be a valid integer");
    return static_cast<std::size_t>(number);
}

}  // namespace

std::filesystem::path runPcaPostprocess(
    const std::string& method,const SourceData& source,
    const json::Value& settings,const std::filesystem::path& root,
    const std::string& label) {
    dynsyslab::core::NpzRecord archive; std::string details;
    if (method=="svd" || method=="pca" || method=="pod") {
        const bool centered=method=="pca" || method=="pod"
            || settings.at("center").boolean(false);
        const auto result=dynsyslab::core::principalComponentAnalysis(
            source.values,source.samples,source.dimension,
            {count(settings.at("rank"),0,"rank",true),centered});
        addPcaArrays(archive,result);
        archive.arrays["left_singular_vectors"]=
            dynsyslab::core::NpyArray::fromDoubles(
                {result.samples,result.rank},[&] {
                    std::vector<double> left(result.scores.size());
                    for (std::size_t row=0;row<result.samples;++row)
                        for (std::size_t mode=0;mode<result.rank;++mode) {
                            const double singular=result.singularValues[mode];
                            left[row*result.rank+mode]=singular>0.0
                                ?result.scores[row*result.rank+mode]/singular:0.0;
                        }
                    return left;
                }());
        archive.arrays["right_singular_vectors"]=archive.arrays["components"];
        archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples},source.time);
        details="\"centered\":"+std::string(centered?"true":"false")
            +",\"rank\":"+std::to_string(result.rank);
    } else if (method=="bagged_pca") {
        dynsyslab::core::BaggedPCAOptions options;
        options.rank=count(settings.at("rank"),0,"rank",true);
        options.center=settings.at("bagged_center").boolean(
            settings.at("center").boolean(true));
        options.bags=count(settings.at("bags"),32,"bags");
        if (!settings.at("seed").isNull())
            options.seed=static_cast<std::uint64_t>(
                count(settings.at("seed"),0,"seed",true));
        const auto result=dynsyslab::core::baggedPrincipalComponentAnalysis(
            source.values,source.samples,source.dimension,options);
        addPcaArrays(archive,result.decomposition);
        archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples},source.time);
        details="\"centered\":"+std::string(options.center?"true":"false")
            +",\"rank\":"+std::to_string(result.decomposition.rank)
            +",\"bags\":"+std::to_string(result.bags)
            +",\"seed\":"+std::to_string(result.seed)
            +",\"rng\":\"std::mt19937_64\"";
    } else if (method=="robust_pca") {
        dynsyslab::core::RobustPCAOptions options;
        options.lambda=settings.at("lambda").number(
            settings.at("lambda_value").number(0.0));
        options.tolerance=settings.at("tolerance").number(1.0e-7);
        options.maximumIterations=count(
            settings.at("max_iterations"),1000,"max_iterations");
        const auto result=dynsyslab::core::robustPrincipalComponentAnalysis(
            source.values,source.samples,source.dimension,options);
        archive.arrays["low_rank"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples,source.dimension},result.lowRank);
        archive.arrays["sparse"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples,source.dimension},result.sparse);
        archive.arrays["iterations"]=dynsyslab::core::NpyArray::fromDoubles(
            {1},{static_cast<double>(result.iterations)});
        archive.arrays["relative_residual"]=
            dynsyslab::core::NpyArray::fromDoubles({1},{result.relativeResidual});
        archive.arrays["lambda_value"]=dynsyslab::core::NpyArray::fromDoubles(
            {1},{result.lambda});
        archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples},source.time);
        details="\"iterations\":"+std::to_string(result.iterations)
            +",\"relative_residual\":"+numeric(result.relativeResidual)
            +",\"lambda\":"+numeric(result.lambda)
            +",\"converged\":"+std::string(result.converged?"true":"false");
    } else if (method=="hankel_pod") {
        dynsyslab::core::HankelPODOptions options;
        options.matrixSize=count(settings.at("matrix_size"),10,"matrix_size");
        options.delayStride=count(settings.at("delay_stride"),1,"delay_stride");
        options.rank=count(settings.at("rank"),0,"rank",true);
        options.center=settings.at("center").boolean(false);
        const auto result=dynsyslab::core::hankelProperOrthogonalDecomposition(
            source.values,source.samples,source.dimension,options);
        const auto& decomposition=result.decomposition;
        archive.arrays["left_singular_vectors"]=
            dynsyslab::core::NpyArray::fromDoubles(
                {decomposition.rows,decomposition.rank},decomposition.left);
        archive.arrays["singular_values"]=dynsyslab::core::NpyArray::fromDoubles(
            {decomposition.rank},decomposition.singularValues);
        archive.arrays["right_singular_vectors"]=
            dynsyslab::core::NpyArray::fromDoubles(
                {decomposition.rank,decomposition.columns},decomposition.right);
        archive.arrays["mean"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.embeddedFeatures},result.mean);
        archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.embeddedSamples},std::vector<double>(
                source.time.begin(),source.time.begin()
                    +static_cast<std::ptrdiff_t>(result.embeddedSamples)));
        details="\"centered\":"+std::string(options.center?"true":"false")
            +",\"rank\":"+std::to_string(decomposition.rank)
            +",\"matrix_size\":"+std::to_string(result.matrixSize)
            +",\"delay_stride_samples\":"+std::to_string(result.delayStride)
            +",\"state_dimension\":"+std::to_string(source.dimension)
            +",\"embedded_feature_dimension\":"
            +std::to_string(result.embeddedFeatures)
            +",\"embedded_samples\":"+std::to_string(result.embeddedSamples);
    } else throw std::invalid_argument("unknown PCA/POD stage");
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        metadata(source,method,label,details));
    return publish(root,source,method,label,std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
