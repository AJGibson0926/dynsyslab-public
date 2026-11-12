#include <dynsyslab/core/SpectralKROM.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
using V=json::Value; using O=V::Object; using A=V::Array;
V numbers(const std::vector<double>& values){A a;for(double v:values)a.emplace_back(v);return a;}
V complexes(const std::vector<std::complex<double>>& values){A a;for(auto v:values)a.emplace_back(A{V(v.real()),V(v.imag())});return a;}
std::vector<double> numbers(const V& value){std::vector<double> out;for(const auto& v:value.array()){double x=v.number(std::numeric_limits<double>::quiet_NaN());if(!std::isfinite(x))throw std::invalid_argument("Expected finite numeric array");out.push_back(x);}return out;}
std::vector<std::complex<double>> complexes(const V& value){std::vector<std::complex<double>> out;for(const auto& v:value.array()){auto a=numbers(v);if(a.size()!=2)throw std::invalid_argument("Complex entry needs [real,imag]");out.emplace_back(a[0],a[1]);}return out;}
std::size_t size(const V& v){double n=v.number(-1);if(!std::isfinite(n)||n<0||std::floor(n)!=n||n>double(std::numeric_limits<std::size_t>::max()/2))throw std::invalid_argument("Expected nonnegative integer");return std::size_t(n);}
V residualJson(const SpectralResiduals& r){return O{{"scope",r.scope},{"definition",r.definition},{"eigenpair",numbers(r.eigenpair)},{"minimum",numbers(r.minimum)},{"grid",complexes(r.grid)},{"minimizing_coefficients",complexes(r.minimizingCoefficients)},{"features",double(r.features)}};}
}
json::Value spectralDictionaryJson(const SpectralDictionary& d) {
    A exponents;for(auto v:d.exponents)exponents.emplace_back(double(v));
    return O{{"kind",d.kind},{"dimension",double(d.dimension)},{"constant",d.constant},
        {"center",numbers(d.center)},{"scale",numbers(d.scale)},{"nodes",numbers(d.nodes)},
        {"exponents",exponents},{"bandwidth",d.bandwidth}};
}
SpectralDictionary spectralDictionaryFromJson(const V& v) {
    SpectralDictionary d;d.kind=v.at("kind").string("linear");d.dimension=size(v.at("dimension"));
    d.constant=v.at("constant").boolean(true);d.center=numbers(v.at("center"));d.scale=numbers(v.at("scale"));
    d.nodes=numbers(v.at("nodes"));for(const auto& e:v.at("exponents").array())d.exponents.push_back(size(e));
    d.bandwidth=v.at("bandwidth").number(1);d.validate();return d;
}
json::Value spectralDataJson(const SpectralData& s) {
    A measures,packets,rkhs;
    for(const auto& k:s.rkhsAdjoints)rkhs.emplace_back(O{{"dimension",double(k.dimension)},{"nodes_count",double(k.nodesCount)},
        {"rank",double(k.rank)},{"nodes",numbers(k.nodes)},{"successors",numbers(k.successors)},{"bandwidth",k.bandwidth},
        {"eigenvalues",complexes(k.eigenvalues)},{"section_coefficients",complexes(k.sectionCoefficients)},
        {"residuals",residualJson(k.residuals)},{"assumptions",k.assumptions}});
    for(const auto& m:s.measures)measures.emplace_back(O{{"name",m.name},{"domain",m.domain},{"method",m.method},
        {"assumptions",m.assumptions},{"grid",numbers(m.grid)},{"density",numbers(m.density)},
        {"atom_locations",numbers(m.atomLocations)},{"atom_masses",numbers(m.atomMasses)},
        {"moments",complexes(m.moments)},{"total_mass",m.totalMass},{"bandwidth",m.bandwidth},{"order",double(m.order)}});
    for(const auto& p:s.wavePackets)packets.emplace_back(O{{"observable",p.observable},{"dictionary_identity",p.dictionaryIdentity},
        {"assumptions",p.assumptions},{"grid",numbers(p.grid)},{"coefficients",complexes(p.coefficients)},
        {"features",double(p.features)},{"order",double(p.order)},{"bandwidth",p.bandwidth}});
    return O{{"schema","dynsyslab.spectral-data/v1"},{"measures",measures},{"wave_packets",packets},{"rkhs_adjoints",rkhs},
        {"stochastic",O{{"sample_mean",numbers(s.stochastic.sampleMean)},{"innovation",numbers(s.stochastic.innovation)},
        {"total",numbers(s.stochastic.total)},{"replicates",double(s.stochastic.replicates)},{"interpretation",s.stochastic.interpretation}}},
        {"residuals",residualJson(s.residuals)},{"invariance",O{{"principal_angles",numbers(s.invariance.principalAngles)},
        {"projection_defect",s.invariance.projectionDefect},{"source_rank",double(s.invariance.sourceRank)},
        {"successor_rank",double(s.invariance.successorRank)}}}};
}
SpectralData spectralDataFromJson(const V& v) {
    if(v.at("schema").string()!="dynsyslab.spectral-data/v1")throw std::invalid_argument("Unknown spectral-data schema");
    SpectralData s;
    const auto& stochastic=v.at("stochastic");
    if(!stochastic.isNull()) {
        s.stochastic.sampleMean=numbers(stochastic.at("sample_mean"));s.stochastic.innovation=numbers(stochastic.at("innovation"));
        s.stochastic.total=numbers(stochastic.at("total"));s.stochastic.replicates=size(stochastic.at("replicates"));
        s.stochastic.interpretation=stochastic.at("interpretation").string();
        if(s.stochastic.sampleMean.size()!=s.stochastic.innovation.size()||s.stochastic.total.size()!=s.stochastic.sampleMean.size())
            throw std::invalid_argument("Stochastic residual shape mismatch");
    }
    for(const auto& a:v.at("measures").array()) {
        SpectralMeasure m;m.name=a.at("name").string();m.domain=a.at("domain").string();m.method=a.at("method").string();
        m.assumptions=a.at("assumptions").string();m.grid=numbers(a.at("grid"));m.density=numbers(a.at("density"));
        m.atomLocations=numbers(a.at("atom_locations"));m.atomMasses=numbers(a.at("atom_masses"));m.moments=complexes(a.at("moments"));
        m.totalMass=a.at("total_mass").number();m.bandwidth=a.at("bandwidth").number();m.order=size(a.at("order"));
        if(m.grid.size()!=m.density.size()||m.atomLocations.size()!=m.atomMasses.size())throw std::invalid_argument("Spectral measure shape mismatch");
        s.measures.push_back(std::move(m));
    }
    for(const auto& a:v.at("rkhs_adjoints").array()) {
        RkhsAdjointSpectrum k;k.dimension=size(a.at("dimension"));k.nodesCount=size(a.at("nodes_count"));k.rank=size(a.at("rank"));
        k.nodes=numbers(a.at("nodes"));k.successors=numbers(a.at("successors"));k.bandwidth=a.at("bandwidth").number();
        k.eigenvalues=complexes(a.at("eigenvalues"));k.sectionCoefficients=complexes(a.at("section_coefficients"));k.assumptions=a.at("assumptions").string();
        const auto& r=a.at("residuals");k.residuals.scope=r.at("scope").string();k.residuals.definition=r.at("definition").string();
        k.residuals.eigenpair=numbers(r.at("eigenpair"));k.residuals.minimum=numbers(r.at("minimum"));k.residuals.grid=complexes(r.at("grid"));
        k.residuals.minimizingCoefficients=complexes(r.at("minimizing_coefficients"));k.residuals.features=size(r.at("features"));
        if(!k.dimension||!k.rank||k.rank>k.nodesCount||k.nodes.size()!=k.nodesCount*k.dimension||k.successors.size()!=k.nodes.size()||
            k.eigenvalues.size()!=k.rank||k.sectionCoefficients.size()!=k.nodesCount*k.rank)throw std::invalid_argument("RKHS section spectrum shape mismatch");
        s.rkhsAdjoints.push_back(std::move(k));
    }
    for(const auto& a:v.at("wave_packets").array()) {
        SpectralWavePackets p;p.observable=a.at("observable").string();p.dictionaryIdentity=a.at("dictionary_identity").string();p.assumptions=a.at("assumptions").string();
        p.grid=numbers(a.at("grid"));p.coefficients=complexes(a.at("coefficients"));p.features=size(a.at("features"));p.order=size(a.at("order"));p.bandwidth=a.at("bandwidth").number();
        if(p.coefficients.size()!=p.grid.size()*p.features)throw std::invalid_argument("Wave packet shape mismatch");
        s.wavePackets.push_back(std::move(p));
    }
    const auto& r=v.at("residuals");s.residuals.scope=r.at("scope").string();s.residuals.definition=r.at("definition").string();
    s.residuals.eigenpair=numbers(r.at("eigenpair"));s.residuals.minimum=numbers(r.at("minimum"));s.residuals.grid=complexes(r.at("grid"));
    s.residuals.minimizingCoefficients=complexes(r.at("minimizing_coefficients"));s.residuals.features=size(r.at("features"));
    if(s.residuals.grid.size()!=s.residuals.minimum.size()||s.residuals.minimizingCoefficients.size()!=s.residuals.grid.size()*s.residuals.features)
        throw std::invalid_argument("Residual grid shape mismatch");
    const auto& i=v.at("invariance");s.invariance.principalAngles=numbers(i.at("principal_angles"));
    s.invariance.projectionDefect=i.at("projection_defect").number();s.invariance.sourceRank=size(i.at("source_rank"));s.invariance.successorRank=size(i.at("successor_rank"));
    return s;
}
json::Value spectralKromJson(const SpectralKROM& m) {
    m.validate();O readouts;for(const auto& r:m.reconstructions)readouts[r.first]=O{{"outputs",double(r.second.outputs)},
        {"coefficients",numbers(r.second.coefficients)},{"relative_training_error",r.second.relativeTrainingError}};
    return O{{"schema","dynsyslab.spectral-krom/v1"},{"dictionary",spectralDictionaryJson(m.dictionary)},
        {"method",m.method},{"assumptions",m.assumptions},{"coordinate_contract",m.coordinateContract},
        {"sample_interval",m.sampleInterval},{"rank",double(m.rank)},{"samples",double(m.samples)},
        {"basis",numbers(m.basis)},{"action",numbers(m.action)},{"multipliers",complexes(m.multipliers)},
        {"eigenfunctions",complexes(m.eigenfunctions)},{"spectral",spectralDataJson(m.spectral)},{"reconstructions",readouts}};
}
SpectralKROM spectralKromFromJson(const V& v) {
    if(v.at("schema").string()!="dynsyslab.spectral-krom/v1")throw std::invalid_argument("Unknown spectral KROM schema");
    SpectralKROM m;m.dictionary=spectralDictionaryFromJson(v.at("dictionary"));m.method=v.at("method").string();
    m.assumptions=v.at("assumptions").string();m.coordinateContract=v.at("coordinate_contract").string();
    m.sampleInterval=v.at("sample_interval").number();m.rank=size(v.at("rank"));m.samples=size(v.at("samples"));
    m.basis=numbers(v.at("basis"));m.action=numbers(v.at("action"));m.multipliers=complexes(v.at("multipliers"));
    m.eigenfunctions=complexes(v.at("eigenfunctions"));m.spectral=spectralDataFromJson(v.at("spectral"));
    for(const auto& r:v.at("reconstructions").object())m.reconstructions[r.first]={size(r.second.at("outputs")),
        numbers(r.second.at("coefficients")),r.second.at("relative_training_error").number()};
    m.validate();return m;
}
} // namespace dynsyslab::core
