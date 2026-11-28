#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Colbrook, Drysdale & Horning, Rigged DMD, arXiv:2405.00782.
// This exports a causal linear feature map, not an autonomous packet law.
namespace dynsyslab::core {
PeriodicDelayPacketFeatures periodicDelayPacketFeatures(const KoopmanModel& krom,
    std::size_t channel,const std::vector<double>& centers,double bandwidth,std::size_t order){
    if(!krom.compression())throw std::invalid_argument("Packet features need a periodic-delay krom with its compressed operator");
    const auto& frozen=*krom.compression();
    frozen.validate();const auto contract=json::parse(frozen.coordinateContract);
    if(contract.at("schema").string()!="dynsyslab.periodic-delay-coordinates/v1"||
       contract.at("sample_interval").number()!=frozen.sampleInterval||frozen.dictionary.kind!="linear"||
       frozen.dictionary.constant||!frozen.dictionary.center.empty()||!frozen.dictionary.scale.empty()||
       centers.empty()||centers.size()>1024)throw std::invalid_argument("Verified linear periodic-delay model and bounded packet centers required");
    PeriodicDelayPacketFeatures out;out.centers=centers;out.outputs=2*centers.size();
    for(const auto& v:contract.at("channels").array()){double c=v.number(-1);
        if(!std::isfinite(c)||c<0||c>1000000||c!=std::floor(c))throw std::invalid_argument("Invalid packet channel contract");out.channels.push_back(std::size_t(c));}
    for(const auto& v:contract.at("time_offsets").array()){
        if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite packet offsets required");out.timeOffsets.push_back(v.number());}
    const auto target=std::find(out.channels.begin(),out.channels.end(),channel);
    const auto n=frozen.dictionary.size();
    if(out.timeOffsets.empty()||out.channels.empty()||target==out.channels.end()||
       n!=out.timeOffsets.size()*out.channels.size()||frozen.dictionary.dimension!=n||n>50000000/out.outputs)
        throw std::invalid_argument("Packet history dimensions/channel mismatch");
    const auto latest=std::max_element(out.timeOffsets.begin(),out.timeOffsets.end());
    const auto column=std::size_t(latest-out.timeOffsets.begin())*out.channels.size()+std::size_t(target-out.channels.begin());
    const double cutoff=*latest;for(auto& t:out.timeOffsets)t-=cutoff;
    const auto& readout=frozen.readouts.at("coordinates");
    if(readout.outputs!=n)throw std::invalid_argument("Packet map requires coordinate readout");
    std::vector<double> g(n);for(std::size_t i=0;i<n;++i)g[i]=readout.coefficients[i*n+column];
    auto model=krom;riggedSpectralMeasure(model,g,"newest-coordinate",0,"packet-centers",centers,bandwidth,order);
    const auto& packets=model.spectrum().measure("packet-centers","newest-coordinate")->packets->coefficients();
    out.coefficients.resize(n*out.outputs);
    for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<centers.size();++k){
        const auto value=packets[k*n+i];
        if(!std::isfinite(value.real())||!std::isfinite(value.imag()))throw std::runtime_error("Nonfinite packet map");
        out.coefficients[i*out.outputs+2*k]=value.real();out.coefficients[i*out.outputs+2*k+1]=value.imag();}
    return out;
}
}
