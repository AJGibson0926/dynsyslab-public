#include "HankelAnalysis.hpp"

#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/PCA.hpp>

#include <complex>
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

// POD modes as functions of the state: the modal amplitudes a_i(x_k) = sigma_i U_ki at the sample
// states, interpolated over the states by inverse distance (8 neighbours, weights 1/d^2; Shepard
// 1968, Proc. 23rd ACM National Conference; core FunctionInterpolation), give phi_i(x) and the
// reconstruction mean + sum_i phi_i(x) v_i at any state (Holmes, Lumley & Berkooz, Turbulence,
// Coherent Structures, Dynamical Systems and Symmetry, 1996). "amplitude_queries": [[x...], ...].
void addAmplitudeQueries(dynsyslab::core::NpzRecord& archive, const json::Value& queries,
                         const std::vector<double>& states, const dynsyslab::core::PCAResult& result) {
    if (queries.isNull()) return;
    const auto d = result.features, r = result.rank;
    std::vector<std::complex<double>> amplitudes(result.scores.begin(), result.scores.end());
    const auto interpolant = dynsyslab::core::interpolateSamples(states, result.samples,
        dynsyslab::core::FunctionDomain::physicalState(d), std::move(amplitudes), {}, 8, "pod-amplitude");
    std::vector<double> phi, reconstruction;
    for (const auto& query : queries.array()) {
        const auto& point = query.array();
        if (point.size() != d) throw std::invalid_argument("amplitude_queries must be states of the source dimension");
        std::vector<double> x;
        for (const auto& v : point) x.push_back(v.number());
        const auto a = interpolant->evaluate(x.data());
        for (std::size_t i = 0; i < r; ++i) phi.push_back(a[i].real());
        for (std::size_t j = 0; j < d; ++j) {
            double value = result.mean.empty() ? 0.0 : result.mean[j];
            for (std::size_t i = 0; i < r; ++i) value += a[i].real() * result.components[i * d + j];
            reconstruction.push_back(value);
        }
    }
    const auto n = queries.array().size();
    archive.arrays["pod_amplitude_queries"] = dynsyslab::core::NpyArray::fromDoubles({n, r}, phi);
    archive.arrays["pod_reconstruction_queries"] = dynsyslab::core::NpyArray::fromDoubles({n, d}, reconstruction);
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
            source.values(),source.samples(),source.dimension(),
            {count(settings.at("rank"),0,"rank",true),centered});
        addPcaArrays(archive,result);
        addAmplitudeQueries(archive,settings.at("amplitude_queries"),source.values(),result);
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
            {source.samples()},source.time());
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
            source.values(),source.samples(),source.dimension(),options);
        addPcaArrays(archive,result.decomposition);
        archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples()},source.time());
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
            source.values(),source.samples(),source.dimension(),options);
        archive.arrays["low_rank"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples(),source.dimension()},result.lowRank);
        archive.arrays["sparse"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples(),source.dimension()},result.sparse);
        archive.arrays["iterations"]=dynsyslab::core::NpyArray::fromDoubles(
            {1},{static_cast<double>(result.iterations)});
        archive.arrays["relative_residual"]=
            dynsyslab::core::NpyArray::fromDoubles({1},{result.relativeResidual});
        archive.arrays["lambda_value"]=dynsyslab::core::NpyArray::fromDoubles(
            {1},{result.lambda});
        archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {source.samples()},source.time());
        details="\"iterations\":"+std::to_string(result.iterations)
            +",\"relative_residual\":"+numeric(result.relativeResidual)
            +",\"lambda\":"+numeric(result.lambda)
            +",\"converged\":"+std::string(result.converged?"true":"false");
    } else throw std::invalid_argument("unknown PCA/POD stage");
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        metadata(source,method,label,details));
    return publishObject(root,source,dynsyslab::core::object_kind::coordinateMap,method,std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
