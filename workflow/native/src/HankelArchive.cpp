#include "HankelAnalysis.hpp"
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c=dynsyslab::core;
void archiveHankelLayout(c::NpzRecord& archive,const HankelAnalysis& input,bool saveMatrix) {
    const auto& h=input.matrix;
    const auto ticks=[&](const char* key,const std::vector<std::size_t>& values) {
        archive.arrays[key]=c::NpyArray::fromInt64s({values.size()},{values.begin(),values.end()});
    };
    ticks("hankel_start_samples",h.starts);
    ticks("hankel_end_samples",h.ends);
    if (saveMatrix) {
        archive.arrays["hankel_matrix"]=c::NpyArray::fromDoubles({h.rows,h.columns},h.x);
        if (!h.y.empty()) archive.arrays["shifted_hankel_matrix"]=c::NpyArray::fromDoubles({h.rows,h.columns},h.y);
    }
}

c::NpzRecord hankelMatrixObject(const SourceData& source,const HankelAnalysis& input,const std::string& label) {
    c::NpzRecord archive;
    archiveHankelLayout(archive,input,true);
    const auto& h=input.matrix;
    const auto& ticks=input.alignEnd?h.ends:h.starts;
    std::vector<double> time,successors;
    time.reserve(ticks.size());
    for (auto tick:ticks) {
        time.push_back(source.time()[tick]);
        if (!h.y.empty()) successors.push_back(source.time()[tick+h.window.sampleStride]);
    }
    archive.arrays["time"]=c::NpyArray::fromDoubles({time.size()},time);
    if (!successors.empty()) archive.arrays["successor_time"]=c::NpyArray::fromDoubles({successors.size()},successors);
    auto meta=json::parse(metadata(source,"hankel_matrix",label,""));
    meta.at("hankel_window")=input.provenance;
    meta.at("source_sha256")=c::sha256File(source.path);
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return archive;
}

std::filesystem::path runHankelMatrix(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    if (settings.at("window").isNull() || !settings.at("hankel_archive").isNull())
        throw std::invalid_argument("hankel_matrix requires a window on the supplied source");
    const auto input=prepareHankel(source,settings,settings.at("shifted_pair").boolean(true));
    return publishObject(root, source, dynsyslab::core::object_kind::matrix, "hankel_matrix",hankelMatrixObject(source,input,label));
}

HankelAnalysis readHankel(const SourceData& source,const std::filesystem::path& path,bool requirePair) {
    const auto archive=c::readNpz(path);
    const auto meta=json::parse(archive.metadataJson());
    const auto& layout=meta.at("hankel_window");
    // A Hankel matrix object, whichever step wrote it (hankel_matrix, or H-DMD with save_matrix).
    if (c::dataObjectKind(archive.metadataJson())!=c::object_kind::matrix || layout.at("schema").string()!="dynsyslab.hankel-window/v1")
        throw std::invalid_argument("expected a published Hankel matrix");
    if (meta.at("source_sha256").string()!=c::sha256File(source.path))
        throw std::invalid_argument("Hankel matrix and source trajectory have different content hashes");
    const auto count=[](const json::Value& v,bool zero=false) {
        const double x=v.number(-1.);
        if (!std::isfinite(x) || std::floor(x)!=x || x>9007199254740991. || x<(zero?0.:1.))
            throw std::invalid_argument("invalid archived Hankel dimension or sample offset");
        return static_cast<std::size_t>(x);
    };
    HankelAnalysis result;
    result.provenance=layout;
    result.provenance.at("matrix_archive")=std::filesystem::absolute(path).lexically_normal().generic_u8string();
    result.alignEnd=layout.at("alignment").string()=="end";
    auto& h=result.matrix;
    h.sourceSamples=source.samples();
    h.rows=count(layout.at("rows")); h.columns=count(layout.at("columns"));
    h.window={count(layout.at("delays")),count(layout.at("sample_stride")),count(layout.at("delay_stride")),
        count(layout.at("snapshot_stride")),h.columns};
    for (const auto& item:layout.at("channels").array()) {
        result.names.push_back(item.at("name").string());
        const double scale=item.at("scale").number();
        if (!(scale>0) || !std::isfinite(scale)) throw std::invalid_argument("invalid archived Hankel scale");
        result.scales.push_back(scale);
        h.channels.push_back({count(item.at("measurement_column"),true),count(item.at("lag_samples"),true),scale});
    }
    if (h.channels.empty() || h.rows%h.channels.size() || h.rows/h.channels.size()!=h.window.delays)
        throw std::invalid_argument("inconsistent archived Hankel row layout");
    auto readMatrix=[&](const char* key) {
        const auto& a=archive.at(key);
        if (a.shape!=std::vector<std::size_t>{h.rows,h.columns})
            throw std::invalid_argument("archived Hankel matrix has inconsistent dimensions");
        auto values=a.doubles();
        if (!std::all_of(values.begin(),values.end(),[](double x){return std::isfinite(x);}))
            throw std::invalid_argument("archived Hankel matrix contains nonfinite measurements");
        return values;
    };
    h.x=readMatrix("hankel_matrix");
    if (requirePair) h.y=readMatrix("shifted_hankel_matrix");
    auto offsets=[&](const char* key) {
        const auto& a=archive.at(key);
        if (a.shape!=std::vector<std::size_t>{h.columns}) throw std::invalid_argument("inconsistent Hankel clock length");
        std::vector<std::size_t> result;
        for (double value:a.doubles()) {
            const auto i=count(value,true);
            if (i>=source.samples() || (!result.empty() && i<=result.back())) throw std::invalid_argument("invalid Hankel clock");
            result.push_back(i);
        }
        return result;
    };
    h.starts=offsets("hankel_start_samples");h.ends=offsets("hankel_end_samples");
    if (requirePair && h.window.sampleStride>=source.samples()-h.ends.back())
        throw std::invalid_argument("archived Hankel successor exceeds source clock");
    return result;
}
} // namespace dynsyslab::workflow::postprocess_detail
