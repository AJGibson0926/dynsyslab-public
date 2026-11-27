#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double number(const J&v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");return v.number();}
std::size_t integer(const J&v,std::size_t max=1000000){double n=number(v);if(n<0||n>max||n!=std::floor(n))throw std::invalid_argument("Integer out of range");return std::size_t(n);}
const A&array(const J&v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>1000000)throw std::invalid_argument("Bounded array required");return v.array();}
std::vector<double> reals(const J&v){std::vector<double>r;for(const auto&x:array(v))r.push_back(number(x));return r;}
std::vector<std::size_t> indices(const J&v){std::vector<std::size_t>r;for(const auto&x:array(v))r.push_back(integer(x));return r;}
template<class T>J values(const std::vector<T>&v){A a;for(auto x:v)a.emplace_back(double(x));return a;}
J relative(double error,double reference){return reference>0?J(error/reference):J{};}
J relatives(const std::vector<double>&errors,const std::vector<double>&reference){A a;for(std::size_t i=0;i<errors.size();++i)a.push_back(relative(errors[i],reference[i]));return a;}
std::vector<math::PeriodicSnapshotOrbit> bank(const J&v){std::vector<math::PeriodicSnapshotOrbit>out;
    for(const auto&o:array(v)){math::PeriodicSnapshotOrbit orbit;orbit.weight=number(o.at("weight"));
        for(const auto&channel:array(o.at("channels"))){math::FourierSeries s{number(channel.at("angular_frequency")),{}};
            for(const auto&row:array(channel.at("coefficients"))){const auto&a=array(row);if(a.size()!=3)throw std::invalid_argument("Fourier rows need [n,re,im]");
                double n=number(a[0]);if(std::abs(n)>1000000||n!=std::floor(n))throw std::invalid_argument("Invalid harmonic");
                if(!s.coefficients.emplace(int(n),std::complex<double>(number(a[1]),number(a[2]))).second)throw std::invalid_argument("Duplicate harmonic");}
            orbit.channels.push_back(std::move(s));}out.push_back(std::move(orbit));}return out;}
J auditJson(const math::PeriodicQuadratureAudit&a){A relativeErrors;
    for(std::size_t k=0;k<a.phaseCounts.size();++k)for(std::size_t c=0;c<a.channels;++c)
        relativeErrors.push_back(relative(a.sampledMass[k*a.channels+c]-a.exactMass[c],a.exactMass[c]));
    return O{{"channels",double(a.channels)},{"maximum_harmonic",double(a.maximumHarmonic)},
        {"alias_free_phase_count",double(a.aliasFreePhaseCount)},{"supplied_weight",a.suppliedWeight},
        {"phase_counts",values(a.phaseCounts)},{"exact_mass",values(a.exactMass)},{"sampled_mass",values(a.sampledMass)},
        {"signed_relative_error",relativeErrors},{"layout","phase-count/channel"},
        {"scope","Normalized supplied orbit weights; finite Fourier products only, not omitted tails or transverse quadrature"}};}
std::string identifier(const J&v){const auto s=v.string();if(s.empty()||s.size()>120||
    !std::all_of(s.begin(),s.end(),[](unsigned char c){return std::isalnum(c)||c=='_'||c=='-';}))throw std::invalid_argument("Group name needs letters, digits, hyphens or underscores");return s;}
struct Model {std::string label,archive;std::vector<std::size_t>channels;std::vector<double>offsets;core::SpectralKROM model;};
}
core::json::Value evaluatePeriodicSpectralQuery(const J&q){
    if(q.at("operation").string()=="packet_features"){
        const auto model=q.at("model").isNull()?core::loadSpectralKrom(std::filesystem::u8path(q.at("archive").string())):core::spectralKromFromJson(q.at("model"));
        const auto r=core::periodicDelayPacketFeatures(model,integer(q.at("channel")),reals(q.at("packet_centers")),number(q.at("bandwidth")),integer(q.at("order"),8));
        return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-packet-features"},
            {"model_identity",core::spectralKromIdentity(model)},{"coordinate_contract",model.coordinateContract},
            {"sample_interval",model.sampleInterval},{"packet_centers",values(r.centers)},{"bandwidth",q.at("bandwidth")},{"order",q.at("order")},
            {"features",O{{"channels",values(r.channels)},{"time_offsets",values(r.timeOffsets)},
                {"outputs",double(r.outputs)},{"coefficients",values(r.coefficients)}}},
            {"layout","coefficients input x output; input offset/channel, output center/(real,imaginary)"},
            {"scope","Frozen newest-observable Rigged-DMD packet map. History shifted to end at zero; no future observations, refitting or source mutation."}};
    }
    if(q.at("operation").string()=="delay_forecast_audit"){
        const auto validation=bank(q.at("validation"));std::vector<std::size_t>steps;
        std::vector<std::vector<std::size_t>> selected;
        if(!q.at("evaluation_steps").isNull()){
            if(!q.at("steps_per_orbit").isNull())throw std::invalid_argument("Use evaluation_steps OR steps_per_orbit");
            for(const auto&v:array(q.at("evaluation_steps"))){selected.push_back(indices(v));
                if(selected.back().empty())throw std::invalid_argument("Evaluation-step lists cannot be empty");steps.push_back(selected.back().back());}
        }else steps=indices(q.at("steps_per_orbit"));
        const auto model=q.at("model").isNull()?core::loadSpectralKrom(std::filesystem::u8path(q.at("archive").string())):core::spectralKromFromJson(q.at("model"));
        const auto result=core::auditPeriodicDelayForecast(model,validation,steps,selected);A rows;
        for(std::size_t j=0;j<result.size();++j){const auto&r=result[j];rows.emplace_back(O{{"steps",double(r.steps)},
            {"parameter",q.at("validation").array()[j].at("parameter")},{"channels",double(r.channels)},{"step_indices",values(r.stepIndices)},
            {"reference_mean_square",values(r.referenceMeanSquare)},{"prediction_mean_square",values(r.predictionMeanSquare)},
            {"error_mean_square",values(r.errorMeanSquare)},{"persistence_error_mean_square",values(r.persistenceErrorMeanSquare)}});}
        return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-delay-forecast-audit"},
            {"model_identity",core::spectralKromIdentity(model)},{"sample_interval",model.sampleInterval},{"orbits",rows},
            {"validation_sha256",core::sha256Text(core::json::dump(q.at("validation")))},
            {"layout","one row per orbit; forecast arrays time/channel, including current time"},
            {"scope","Newest-coordinate, past-only free forecasts. Exact uniform-phase Parseval integrals for supplied finite Fourier coefficients, up to roundoff. No temporal interpolation, phase quadrature, omitted-tail bound or continuum-transverse certificate."}};
    }
    const auto training=bank(q.at("training"));std::vector<std::size_t>counts;
    if(!q.at("audit_phase_counts").isNull())counts=indices(q.at("audit_phase_counts"));
    const auto audit=math::auditPeriodicQuadrature(training,counts);
    if(q.at("operation").string()=="quadrature_audit")return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-quadrature-audit"},{"audit",auditJson(audit)}};
    if(q.at("operation").string()=="delay_fit"){
        core::KoopmanSpectralOptions options;options.method=q.at("method").string("edmd");
        options.sampleInterval=number(q.at("sample_interval"));
        options.relativeTolerance=q.at("svd_tolerance").number(1e-10);
        options.rank=q.at("rank").isNull()?0:integer(q.at("rank"));
        options.acknowledgeMeasurePreserving=q.at("assume_measure_preserving").boolean();
        const auto output=q.at("output_archive").string();
        if(!output.empty()&&std::filesystem::exists(std::filesystem::u8path(output)))
            throw std::invalid_argument("Delay-fit archive already exists; use a new destination");
        auto model=core::fitPeriodicDelayModel(training,indices(q.at("channels")),reals(q.at("time_offsets")),
            integer(q.at("phase_count")),options,q.at("require_resolved_quadrature").boolean(true),q.at("storage").string("auto"));
        const O provenance{{"training_sha256",core::sha256Text(core::json::dump(q.at("training")))},
            {"source",q.at("source")},{"coordinate_contract",model.coordinateContract},
            {"method",options.method},{"svd_tolerance",options.relativeTolerance},{"requested_rank",double(options.rank)}};
        if(!output.empty()){
            const auto path=std::filesystem::u8path(output);
            if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
            auto archive=core::makeSpectralKromArchive(model);
            archive.arrays["periodic_fit_provenance_json"]=core::NpyArray::fromString(core::json::dump(provenance));
            core::writeNpzAtomicNew(path,archive);
        }
        O result{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-delay-fit"},
            {"fitting_method",model.method},{"model_identity",core::spectralKromIdentity(model)},
            {"rank",double(model.rank)},{"samples",double(model.samples)},{"features",double(model.dictionary.size())},
            {"coordinate_contract",model.coordinateContract},{"provenance",provenance},{"quadrature_audit",auditJson(audit)}};
        if(!output.empty())result["archive"]=output;
        if(q.at("include_model").boolean(output.empty()))result["model"]=core::spectralKromJson(model);
        return result;
    }
    const auto validation=bank(q.at("validation"));math::auditPeriodicQuadrature(validation,{});
    const double dt=number(q.at("sample_interval")),width=number(q.at("bandwidth"));const auto order=integer(q.at("order"),8);
    const auto angles=reals(q.at("angles")),centers=reals(q.at("packet_centers"));
    if(angles.empty()||angles.size()>16384||centers.empty()||centers.size()>64)throw std::invalid_argument("Spectral grid/packet center limits exceeded");
    std::vector<Model>models;
    if(q.at("operation").string()=="spectral_sweep"){
        const auto delays=indices(q.at("delay_counts"));const auto phases=integer(q.at("phase_count"));
        if(delays.empty()||delays.size()>32||array(q.at("groups")).empty()||array(q.at("groups")).size()>32)throw std::invalid_argument("Bounded nonempty sweep required");
        core::KoopmanSpectralOptions options;options.method="mpedmd";options.sampleInterval=dt;
        options.acknowledgeMeasurePreserving=q.at("assume_measure_preserving").boolean();
        options.relativeTolerance=q.at("svd_tolerance").number(1e-10);options.rank=q.at("rank").isNull()?0:integer(q.at("rank"));
        if(!options.acknowledgeMeasurePreserving)throw std::invalid_argument("Assert invariant-measure sampling explicitly");
        const auto output=q.at("output_directory").string();
        const auto requestArchive=std::filesystem::u8path(output)/"periodic-spectral-request.npz";
        if(!output.empty()&&std::filesystem::exists(requestArchive))throw std::invalid_argument("Sweep input archive already exists; use a new output directory or compare frozen models");
        // Preflight every immutable destination before fitting/writing anything.
        for(const auto&g:array(q.at("groups")))for(auto n:delays){if(n<1||n>4096)throw std::invalid_argument("Delay count must be 1..4096");
            Model m;m.label=identifier(g.at("name"))+"-"+std::to_string(n);m.channels=indices(g.at("channels"));
            for(std::size_t k=0;k<n;++k)m.offsets.push_back((double(k)-double(n/2))*dt);
            if(!output.empty())m.archive=(std::filesystem::u8path(output)/(m.label+".npz")).u8string();
            if(!m.archive.empty()&&std::filesystem::exists(std::filesystem::u8path(m.archive)))throw std::invalid_argument("Spectral sweep output already exists; use spectral_compare on frozen archives");
            for(const auto&previous:models)if(previous.label==m.label)throw std::invalid_argument("Duplicate sweep label");
            models.push_back(std::move(m));}
        if(!output.empty()){
            std::filesystem::create_directories(std::filesystem::u8path(output));core::NpzRecord input;
            input.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(q));core::writeNpzAtomicNew(requestArchive,input);
        }
        for(auto&m:models){m.model=core::fitPeriodicDelaySpectrum(training,m.channels,m.offsets,phases,options,q.at("require_resolved_quadrature").boolean(true),q.at("storage").string("auto"));
            if(!m.archive.empty()){std::filesystem::create_directories(std::filesystem::u8path(output));core::writeNpzAtomicNew(std::filesystem::u8path(m.archive),core::makeSpectralKromArchive(m.model));}}
    }else if(q.at("operation").string()=="spectral_compare"){
        if(array(q.at("models")).empty()||array(q.at("models")).size()>128)throw std::invalid_argument("Need 1..128 frozen models");
        for(const auto&v:array(q.at("models"))){Model m;m.label=v.at("label").string();m.channels=indices(v.at("channels"));m.offsets=reals(v.at("time_offsets"));m.archive=v.at("archive").string();
            m.model=v.at("model").isNull()?core::loadSpectralKrom(std::filesystem::u8path(m.archive)):core::spectralKromFromJson(v.at("model"));models.push_back(std::move(m));}
    }else throw std::invalid_argument("Unknown periodic spectral operation");
    if(models.empty()||models.size()>128)throw std::invalid_argument("Need 1..128 models");
    std::vector<std::size_t>channels;for(const auto&m:models)for(auto c:m.channels)if(std::find(channels.begin(),channels.end(),c)==channels.end())channels.push_back(c);
    A references,results,modelRecords,frequencies,parameters,families;
    for(std::size_t i=0;i<validation.size();++i){frequencies.emplace_back(validation[i].channels.front().angularFrequency);parameters.push_back(q.at("validation").array()[i].at("parameter"));families.emplace_back(q.at("validation").array()[i].at("family").string("all"));}
    for(auto c:channels){const auto ref=core::periodicSpectralReference(training,validation,c,dt,angles,centers,width,order);
        const auto&names=q.at("channel_names").array();const std::string name=c<names.size()?names[c].string():"channel_"+std::to_string(c);
        references.emplace_back(O{{"channel",double(c)},{"name",name},{"density",values(ref.spectrum.density)},
            {"total_mass",ref.spectrum.totalMass},{"zero_harmonic_mass",ref.spectrum.zeroHarmonicMass},{"packet_rms",values(ref.packetRms)}});
        for(const auto&m:models){if(std::find(m.channels.begin(),m.channels.end(),c)==m.channels.end())continue;
            const auto r=core::comparePeriodicDelaySpectrum(m.model,m.channels,m.offsets,validation,ref);
            results.emplace_back(O{{"channel",double(c)},{"observable",name},{"model",m.label},{"rank",double(r.rank)},{"features",double(r.features)},
                {"fitted_mass",r.fittedMass},{"density",values(r.density)},{"packet_rms",values(r.packetRms)},{"packet_error_rms",values(r.packetErrorRms)},
                {"density_absolute_l2",r.densityErrorNorm},{"density_relative_l2",relative(r.densityErrorNorm,r.densityReferenceNorm)},
                {"packet_amplitude_error",values(r.amplitudeErrorNorm)},{"packet_complex_error",values(r.complexErrorNorm)},
                {"packet_reference_norm",values(r.packetReferenceNorm)},{"packet_amplitude_relative_l2",relatives(r.amplitudeErrorNorm,r.packetReferenceNorm)},
                {"packet_complex_relative_l2",relatives(r.complexErrorNorm,r.packetReferenceNorm)}});
        }
    }
    for(const auto&m:models){O r{{"label",m.label},{"archive",m.archive},{"model_identity",core::spectralKromIdentity(m.model)},
        {"channels",values(m.channels)},{"time_offsets",values(m.offsets)},{"rank",double(m.model.rank)},
        {"coordinate_contract",m.model.coordinateContract}};
        if(q.at("include_models").boolean())r["model"]=core::spectralKromJson(m.model);modelRecords.emplace_back(r);}
    return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-spectral-comparison"},
        {"angles",values(angles)},{"packet_centers",values(centers)},{"sample_interval",dt},{"bandwidth",width},{"order",double(order)},
        {"validation_frequencies",frequencies},{"validation_parameters",parameters},{"validation_families",families},{"references",references},{"comparisons",results},
        {"models",modelRecords},{"quadrature_audit",auditJson(audit)},
        {"training_sha256",core::sha256Text(core::json::dump(q.at("training")))},
        {"validation_sha256",core::sha256Text(core::json::dump(q.at("validation")))},{"source",q.at("source")},
        {"layout","packet arrays center-major/orbit; errors use normalized supplied validation weights; density error is discrete grid L2"},
        {"assumptions","Finite Fourier reference; physical-time phase; invariant product measure asserted by fitter. No unit-mass normalization of observables, clipping, or infinite-dimensional convergence claim."}};
}
}
