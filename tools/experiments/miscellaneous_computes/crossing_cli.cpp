#include "crossing_cli.hpp"
#include "crossing_kernel.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace core=dynsyslab::core;
void crossingTable(const std::filesystem::path& input,const std::filesystem::path& output) {
    auto record=core::readNpz(input);
    const auto t=record.at("time").doubles(),xy=record.at("values").doubles();
    if(record.at("time").shape!=std::vector<std::size_t>{t.size()} ||
       record.at("values").shape!=std::vector<std::size_t>{t.size(),2})
        throw std::invalid_argument("Expected aligned time and two-dimensional plane");
    for(std::size_t i=0;i<t.size();++i)if(!std::isfinite(t[i])||(i&&t[i]<=t[i-1]))
        throw std::invalid_argument("Expected an increasing clock");
    std::size_t count=0;
    if(mc_crossings_v1(xy.data(),t.size(),nullptr,0,&count))throw std::runtime_error("Invalid crossing input");
    std::vector<mc_crossing_v1> hits(count);
    if(mc_crossings_v1(xy.data(),t.size(),hits.data(),count,&count))throw std::runtime_error("Crossing enumeration failed");
    std::vector<double> times,points,segments;
    for(auto h:hits) {
        times.push_back(t[h.first]+h.u*(t[h.first+1]-t[h.first]));
        times.push_back(t[h.second]+h.v*(t[h.second+1]-t[h.second]));
        for(int j=0;j<2;++j)points.push_back(xy[2*h.first+j]+h.u*(xy[2*h.first+2+j]-xy[2*h.first+j]));
        segments.push_back(static_cast<double>(h.first));segments.push_back(static_cast<double>(h.second));
    }
    core::NpzRecord result;
    result.arrays["times"]=core::NpyArray::fromDoubles({count,2},times);
    result.arrays["points"]=core::NpyArray::fromDoubles({count,2},points);
    result.arrays["segments"]=core::NpyArray::fromDoubles({count,2},segments);
    result.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(core::json::Value::Object{
        {"schema","dynsyslab.miscellaneous-computes/planar-crossings/v1"},
        {"source",std::filesystem::absolute(input).generic_u8string()}}));
    core::writeNpzAtomicNew(output,result);
}
