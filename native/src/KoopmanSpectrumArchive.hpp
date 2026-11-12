#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core::detail {
inline void readPolynomialSpectrumMetadata(const NpzRecord& archive,KoopmanSINDyResult& fit) {
    if(archive.contains("eigenvalue_sample_time")) {
        fit.sampleTime=archive.at("eigenvalue_sample_time").doubles().at(0);
        if(!std::isfinite(fit.sampleTime)||fit.sampleTime<=0)
            throw std::invalid_argument("invalid polynomial eigenvalue sampling interval");
    } else if(!archive.metadataJson().empty()) {
        const auto metadata=json::parse(archive.metadataJson());
        // Older native SINDy archives already record the sampling interval in
        // metadata rather than a dedicated array. Never infer it from physics.
        const auto& dt=metadata.at("delta_t");
        if(!dt.isNull()) {
            fit.sampleTime=dt.number(-1);
            if(!std::isfinite(fit.sampleTime)||fit.sampleTime<=0)
                throw std::invalid_argument("invalid archived polynomial sampling interval");
        }
    }
    if(archive.contains("eigenvalue_origin")) fit.eigenvalueOrigin=archive.at("eigenvalue_origin").string();
    if(archive.contains("eigenvalue_log_branches")) {
        for(const double branch:archive.at("eigenvalue_log_branches").doubles()) {
            if(!std::isfinite(branch)||branch!=std::floor(branch)||std::abs(branch)>2147483647)
                throw std::invalid_argument("invalid eigenvalue logarithm branch");
            fit.logarithmBranches.push_back(static_cast<int>(branch));
        }
        if(fit.logarithmBranches.size()!=fit.continuousEigenvalues.size())
            throw std::invalid_argument("eigenvalue logarithm branches are not mode-aligned");
    }
}
inline void writePolynomialSpectrumMetadata(NpzRecord& archive,const KoopmanSINDyResult& fit) {
    if(!std::isfinite(fit.sampleTime)||fit.sampleTime<0)
        throw std::invalid_argument("invalid polynomial eigenvalue sampling interval");
    if(fit.sampleTime>0) {
        archive.arrays["eigenvalue_sample_time"]=NpyArray::fromDoubles({1},{fit.sampleTime});
        archive.arrays["discrete_eigenvalues"]=NpyArray::fromComplex(
            {fit.continuousEigenvalues.size()},fit.discreteEigenvalues());
    }
    archive.arrays["eigenvalue_origin"]=NpyArray::fromString(fit.eigenvalueOrigin);
    if(!fit.logarithmBranches.empty()) {
        if(fit.logarithmBranches.size()!=fit.continuousEigenvalues.size())
            throw std::invalid_argument("eigenvalue logarithm branches are not mode-aligned");
        archive.arrays["eigenvalue_log_branches"]=NpyArray::fromInt64s(
            {fit.logarithmBranches.size()},{fit.logarithmBranches.begin(),fit.logarithmBranches.end()});
    }
}
} // namespace dynsyslab::core::detail
