#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <dynsyslab/numerics/StreamingQr.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Fourier translation/Parseval: Katznelson, Introduction to Harmonic Analysis,
// Ch. I. EDMD: Williams/Kevrekidis/Rowley, https://arxiv.org/abs/1408.4408.
// Ch. I. mpEDMD and Rigged DMD: Colbrook (2023); Colbrook/Drysdale/Horning,
// https://arxiv.org/abs/2405.00782. This module composes existing algorithms.
namespace dynsyslab::core {
namespace {
std::vector<math::PeriodicSnapshotOrbit> select(
    const std::vector<math::PeriodicSnapshotOrbit>&bank,const std::vector<std::size_t>&channels){
    if(channels.empty())throw std::invalid_argument("Select at least one channel");
    auto sorted=channels;std::sort(sorted.begin(),sorted.end());
    if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::invalid_argument("Duplicate selected channel");
    std::vector<math::PeriodicSnapshotOrbit> out;
    for(const auto&o:bank){math::PeriodicSnapshotOrbit row;row.weight=o.weight;
        for(auto c:channels){if(c>=o.channels.size())throw std::invalid_argument("Selected channel out of range");row.channels.push_back(o.channels[c]);}
        out.push_back(std::move(row));}return out;
}
KoopmanSnapshotPairs compressedPairs(const std::vector<math::PeriodicSnapshotOrbit>&bank,
    std::size_t phases,const std::vector<double>&offsets,double dt,std::size_t&jointFeatures){
    const auto channels=bank.front().channels.size();
    // Reuse only EXACTLY identical offsets. No near-equality threshold changes
    // the sampled dictionary. Regular delay grids often share successor taps.
    std::vector<double> all;std::vector<std::size_t> source,successor;
    const auto insert=[&](double t){if(!std::isfinite(t))throw std::invalid_argument("Finite delay offsets required");
        const auto it=std::find(all.begin(),all.end(),t);if(it!=all.end())return std::size_t(it-all.begin());
        all.push_back(t);return all.size()-1;};
    for(double t:offsets)source.push_back(insert(t));
    for(double t:offsets)successor.push_back(insert(t+dt));
    jointFeatures=all.size()*channels;numerics::WeightedQrAccumulator qr(jointFeatures);
    for(const auto&o:bank){const auto block=math::periodicDelaySamples(o.channels,phases,all);
        qr.append(block,phases,std::vector<double>(phases,o.weight/double(phases)));}
    KoopmanSnapshotPairs out;out.samples=qr.rows();out.dimension=offsets.size()*channels;
    // The downstream fitter normalizes supplied row weights to one. These
    // synthetic feature rows therefore need sqrt(compressed_rows/total_mass).
    const double scale=std::sqrt(double(out.samples)/qr.totalWeight());
    out.x.resize(out.samples*out.dimension);out.y.resize(out.x.size());
    for(std::size_t i=0;i<out.samples;++i)for(std::size_t tap=0;tap<offsets.size();++tap)
        for(std::size_t c=0;c<channels;++c){const auto j=i*out.dimension+tap*channels+c;
            out.x[j]=scale*qr.factor()[i*jointFeatures+source[tap]*channels+c];
            out.y[j]=scale*qr.factor()[i*jointFeatures+successor[tap]*channels+c];}
    return out;
}
}
SpectralKROM fitPeriodicDelaySpectrum(const std::vector<math::PeriodicSnapshotOrbit>&bank,
    const std::vector<std::size_t>&channels,const std::vector<double>&offsets,
    std::size_t phases,const KoopmanSpectralOptions&options,bool requireResolved,const std::string&storage){
    if(options.method!="mpedmd")throw std::invalid_argument("Periodic spectral comparison requires mpEDMD");
    return fitPeriodicDelayModel(bank,channels,offsets,phases,options,requireResolved,storage);
}
SpectralKROM fitPeriodicDelayModel(const std::vector<math::PeriodicSnapshotOrbit>&bank,
    const std::vector<std::size_t>&channels,const std::vector<double>&offsets,
    std::size_t phases,const KoopmanSpectralOptions&options,bool requireResolved,const std::string&storage){
    if(options.method!="edmd"&&options.method!="resdmd"&&options.method!="mpedmd")
        throw std::invalid_argument("Periodic delay fitting supports EDMD, ResDMD and mpEDMD");
    if(options.weightPolicy!="uniform")throw std::invalid_argument("Use the explicitly supplied product weights without an additional weight policy");
    const auto selected=select(bank,channels);const auto audit=math::auditPeriodicQuadrature(selected,{});
    if(phases<3||phases>1000000||selected.size()>1000000/phases||offsets.empty()||offsets.size()>4096
        ||channels.size()>4096/offsets.size()||!std::isfinite(options.sampleInterval)||options.sampleInterval<=0)
        throw std::invalid_argument("Periodic fit requires at most one million samples and 4096 features, with a positive clock");
    if(storage!="auto"&&storage!="dense"&&storage!="streaming_qr")throw std::invalid_argument("Storage must be auto, dense, or streaming_qr");
    if(requireResolved&&phases<audit.aliasFreePhaseCount)
        throw std::invalid_argument("Phase quadrature under-resolves supplied Fourier products; increase phase_count or explicitly set require_resolved_quadrature=false");
    const auto originalSamples=phases*selected.size(),dimension=offsets.size()*channels.size();
    const bool streaming=storage=="streaming_qr"||(storage=="auto"&&originalSamples>2000000/dimension);
    std::size_t jointFeatures=2*dimension;KoopmanSnapshotPairs data;
    if(streaming)data=compressedPairs(selected,phases,offsets,options.sampleInterval,jointFeatures);
    else{auto p=math::periodicSnapshotPairs(selected,phases,offsets,options.sampleInterval);
        data.x=std::move(p.x);data.y=std::move(p.y);data.weights=std::move(p.weights);data.samples=p.samples;data.dimension=p.dimension;}
    SpectralDictionary dictionary;dictionary.dimension=data.dimension;dictionary.constant=false;
    auto model=fitKoopmanSpectral(data,dictionary,options);
    fitSpectralReconstruction(model,"coordinates",data.x,data.x,data.dimension,data.weights);
    model.samples=originalSamples;
    if(streaming)model.assumptions+=" Weighted row-block QR preserves the joint feature products; training residuals use equivalent compressed feature rows, not physical states.";
    json::Value::Array cs,ts;for(auto c:channels)cs.emplace_back(double(c));for(auto t:offsets)ts.emplace_back(t);
    model.coordinateContract=json::dump(json::Value::Object{{"schema","dynsyslab.periodic-delay-coordinates/v1"},
        {"channels",cs},{"time_offsets",ts},{"phase_count",double(phases)},{"sample_interval",options.sampleInterval},
        {"resolved_finite_products",phases>=audit.aliasFreePhaseCount},
        {"storage",streaming?"streaming_qr":"dense"},{"original_samples",double(originalSamples)},
        {"fit_rows",double(data.samples)},{"joint_features",double(jointFeatures)}});
    return model;
}
PeriodicSpectralReference periodicSpectralReference(const std::vector<math::PeriodicSnapshotOrbit>&training,
    const std::vector<math::PeriodicSnapshotOrbit>&validation,std::size_t channel,double dt,
    const std::vector<double>&angles,const std::vector<double>&centers,double bandwidth,std::size_t order){
    const auto audit=math::auditPeriodicQuadrature(training,{});math::auditPeriodicQuadrature(validation,{});
    if(channel>=audit.channels||centers.empty())throw std::invalid_argument("Reference channel and packet centers required");
    const math::RationalCircleKernel kernel(bandwidth,order);PeriodicSpectralReference ref;
    ref.channel=channel;ref.validationOrbits=validation.size();ref.sampleInterval=dt;ref.bandwidth=bandwidth;
    ref.order=order;ref.angles=angles;ref.centers=centers;
    std::vector<math::WeightedOrbitSpectrum> orbits;
    for(const auto&o:training)orbits.push_back({o.channels.at(channel),o.weight/audit.suppliedWeight});
    ref.spectrum=math::smoothOrbitSpectra(orbits,dt,angles,kernel);
    for(double center:centers)for(const auto&o:validation)
        ref.packetRms.push_back(std::sqrt(math::meanSquare(math::spectralPacket(o.channels.at(channel),dt,center,kernel))));
    return ref;
}
PeriodicSpectralComparison comparePeriodicDelaySpectrum(const SpectralKROM&frozen,
    const std::vector<std::size_t>&channels,const std::vector<double>&offsets,
    const std::vector<math::PeriodicSnapshotOrbit>&validation,const PeriodicSpectralReference&ref,bool legacy){
    frozen.validate();const auto selected=select(validation,channels);const auto audit=math::auditPeriodicQuadrature(validation,{});
    json::Value contract;try{contract=json::parse(frozen.coordinateContract);}catch(const std::exception&){if(!legacy)throw std::invalid_argument("Legacy coordinate contract cannot be verified; explicit opt-in required");}
    if(contract.at("schema").string()=="dynsyslab.periodic-delay-coordinates/v1"){
        if(contract.at("sample_interval").number()!=frozen.sampleInterval)throw std::invalid_argument("Frozen coordinate/sample-clock mismatch");
        const auto&cs=contract.at("channels").array();const auto&ts=contract.at("time_offsets").array();
        if(cs.size()!=channels.size()||ts.size()!=offsets.size())throw std::invalid_argument("Frozen delay contract shape mismatch");
        for(std::size_t i=0;i<cs.size();++i)if(cs[i].number()!=double(channels[i]))throw std::invalid_argument("Frozen delay channel order mismatch");
        for(std::size_t i=0;i<ts.size();++i)if(ts[i].number()!=offsets[i])throw std::invalid_argument("Frozen time-offset contract mismatch");
    }else if(!legacy)throw std::invalid_argument("Missing verified periodic-delay coordinate contract");
    if(ref.spectrum.density.size()!=ref.angles.size()||ref.packetRms.size()!=ref.centers.size()*validation.size())
        throw std::invalid_argument("Malformed periodic reference dimensions");
    const auto zero=std::find(offsets.begin(),offsets.end(),0.);const auto target=std::find(channels.begin(),channels.end(),ref.channel);
    const auto dimension=offsets.size()*channels.size();
    if(zero==offsets.end()||target==channels.end())throw std::invalid_argument("A current-time reference channel must be present in the delay dictionary");
    if(frozen.method!="mpedmd"||frozen.dictionary.kind!="linear"||frozen.dictionary.constant||
       !frozen.dictionary.center.empty()||!frozen.dictionary.scale.empty()||frozen.dictionary.dimension!=dimension||
       frozen.sampleInterval!=ref.sampleInterval||validation.size()!=ref.validationOrbits)
        throw std::invalid_argument("Frozen periodic model/dictionary/clock mismatch");
    const auto&readout=frozen.reconstructions.at("coordinates");const auto column=(zero-offsets.begin())*channels.size()+(target-channels.begin());
    if(readout.outputs!=dimension)throw std::invalid_argument("Expected coordinate reconstruction");
    std::vector<double> g(dimension);for(std::size_t i=0;i<dimension;++i)g[i]=readout.coefficients[i*dimension+column];
    auto model=frozen;auto grid=ref.angles;grid.insert(grid.end(),ref.centers.begin(),ref.centers.end());
    riggedSpectralMeasure(model,g,"periodic-comparison",grid,ref.bandwidth,ref.order);
    const auto&m=model.spectral.measures.back();const auto&p=model.spectral.wavePackets.back();
    PeriodicSpectralComparison out;out.rank=model.rank;out.features=dimension;out.fittedMass=m.totalMass;
    out.density.assign(m.density.begin(),m.density.begin()+ref.angles.size());
    for(std::size_t j=0;j<out.density.size();++j){out.densityErrorNorm+=std::pow(out.density[j]-ref.spectrum.density[j],2);out.densityReferenceNorm+=std::pow(ref.spectrum.density[j],2);}
    out.densityErrorNorm=std::sqrt(out.densityErrorNorm);out.densityReferenceNorm=std::sqrt(out.densityReferenceNorm);
    const math::RationalCircleKernel kernel(ref.bandwidth,ref.order);
    for(std::size_t k=0;k<ref.centers.size();++k){
        const auto start=(ref.angles.size()+k)*dimension;
        const std::vector<std::complex<double>> coefficients(p.coefficients.begin()+start,p.coefficients.begin()+start+dimension);
        double amplitudeError=0,complexError=0,referenceNorm=0;
        for(std::size_t j=0;j<validation.size();++j){
            const auto predicted=math::periodicDelayReadout(selected[j].channels,offsets,coefficients);
            const auto truth=math::spectralPacket(validation[j].channels.at(ref.channel),ref.sampleInterval,ref.centers[k],kernel);
            auto difference=predicted;for(const auto&[n,a]:truth.coefficients)difference.coefficients[n]-=a;
            const auto norm=std::sqrt(math::meanSquare(predicted)),error=std::sqrt(math::meanSquare(difference));
            const auto exact=ref.packetRms[k*validation.size()+j],weight=validation[j].weight/audit.suppliedWeight;
            out.packetRms.push_back(norm);out.packetErrorRms.push_back(error);
            amplitudeError+=weight*std::pow(norm-exact,2);complexError+=weight*error*error;referenceNorm+=weight*exact*exact;
        }
        out.amplitudeErrorNorm.push_back(std::sqrt(amplitudeError));out.complexErrorNorm.push_back(std::sqrt(complexError));out.packetReferenceNorm.push_back(std::sqrt(referenceNorm));
    }return out;
}
}
