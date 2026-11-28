// JSON of a krom's spectral parts: the frozen dictionary, spectral bands and measures, the fit's diagnostics and
// the compressed operator (KoopmanCompression.hpp). Bulk arrays of a krom are written as NPZ arrays by
// KromArchiveParts.cpp; this JSON is for identities and query results.
#include <dynsyslab/core/KoopmanCompression.hpp>
#include <dynsyslab/core/FileDigest.hpp>
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
}
json::Value spectralBandJson(const SpectralBand& c) {
    O out{{"name",c.name},{"support",c.support},{"method",c.method},{"assumptions",c.assumptions},
        {"grid",numbers(c.grid)},{"bandwidth",c.bandwidth},{"order",double(c.order)}};
    if(!c.continuity.empty())out["continuity"]=c.continuity;
    return out;
}
SpectralBand spectralBandFromJson(const V& a) {
    SpectralBand c;c.name=a.at("name").string();c.continuity=a.at("continuity").string();
    c.support=a.at("support").string("unit-circle");c.method=a.at("method").string();c.assumptions=a.at("assumptions").string();
    c.grid=numbers(a.at("grid"));c.bandwidth=a.at("bandwidth").number(0);c.order=size(a.at("order"));
    return c;
}
json::Value spectralMeasureJson(const SpectralMeasure& m) {
    return O{{"observable",m.observable},{"output",double(m.output)},{"band",m.band},
        {"density",numbers(m.density)},{"moments",complexes(m.moments)},{"total_mass",m.totalMass},{"assumptions",m.assumptions}};
}
SpectralMeasure spectralMeasureFromJson(const V& a) {
    SpectralMeasure m;m.observable=a.at("observable").string();m.output=size(a.at("output"));m.band=a.at("band").string();
    m.density=numbers(a.at("density"));m.moments=complexes(a.at("moments"));m.totalMass=a.at("total_mass").number();
    m.assumptions=a.at("assumptions").string();
    return m;
}
json::Value spectralDictionaryJson(const SpectralDictionary& d) {
    A exponents;for(auto v:d.exponents)exponents.emplace_back(double(v));
    O result{{"kind",d.kind},{"dimension",double(d.dimension)},{"constant",d.constant},
        {"center",numbers(d.center)},{"scale",numbers(d.scale)},{"nodes",numbers(d.nodes)},
        {"exponents",exponents},{"bandwidth",d.bandwidth}};
    // Preserve byte-for-byte identities of existing v1 dictionaries.
    if(d.kind=="matern")result["smoothness"]=d.smoothness;
    if(d.amplitude!=1)result["amplitude"]=d.amplitude;
    return result;
}
SpectralDictionary spectralDictionaryFromJson(const V& v) {
    SpectralDictionary d;d.kind=v.at("kind").string("linear");d.dimension=size(v.at("dimension"));
    d.constant=v.at("constant").boolean(true);d.center=numbers(v.at("center"));d.scale=numbers(v.at("scale"));
    d.nodes=numbers(v.at("nodes"));for(const auto& e:v.at("exponents").array())d.exponents.push_back(size(e));
    d.bandwidth=v.at("bandwidth").number(1);d.smoothness=v.at("smoothness").number(2);
    d.amplitude=v.at("amplitude").number(1);d.validate();return d;
}
json::Value spectralDiagnosticsJson(const SpectralDiagnostics& s) {
    A rkhs;
    for(const auto& k:s.rkhsAdjoints)rkhs.emplace_back(O{{"dimension",double(k.dimension)},{"nodes_count",double(k.nodesCount)},
        {"rank",double(k.rank)},{"nodes",numbers(k.nodes)},{"successors",numbers(k.successors)},{"bandwidth",k.bandwidth},
        {"eigenvalues",complexes(k.eigenvalues)},{"section_coefficients",complexes(k.sectionCoefficients)},
        {"residuals",numbers(k.residuals)},{"pseudospectrum_grid",complexes(k.pseudospectrum.grid)},
        {"pseudospectrum_minimum",numbers(k.pseudospectrum.minimum)},
        {"pseudospectrum_minimizers",complexes(k.pseudospectrumMinimizers)},{"assumptions",k.assumptions}});
    return O{{"schema","dynsyslab.spectral-diagnostics/v1"},{"rkhs_adjoints",rkhs},
        {"stochastic",O{{"sample_mean",numbers(s.stochastic.sampleMean)},{"innovation",numbers(s.stochastic.innovation)},
        {"total",numbers(s.stochastic.total)},{"replicates",double(s.stochastic.replicates)},{"interpretation",s.stochastic.interpretation}}},
        {"invariance",O{{"principal_angles",numbers(s.invariance.principalAngles)},
        {"projection_defect",s.invariance.projectionDefect},{"source_rank",double(s.invariance.sourceRank)},
        {"successor_rank",double(s.invariance.successorRank)}}}};
}
SpectralDiagnostics spectralDiagnosticsFromJson(const V& v) {
    SpectralDiagnostics s;
    const auto& stochastic=v.at("stochastic");
    if(!stochastic.isNull()) {
        s.stochastic.sampleMean=numbers(stochastic.at("sample_mean"));s.stochastic.innovation=numbers(stochastic.at("innovation"));
        s.stochastic.total=numbers(stochastic.at("total"));s.stochastic.replicates=size(stochastic.at("replicates"));
        s.stochastic.interpretation=stochastic.at("interpretation").string();
        if(s.stochastic.sampleMean.size()!=s.stochastic.innovation.size()||s.stochastic.total.size()!=s.stochastic.sampleMean.size())
            throw std::invalid_argument("Stochastic residual shape mismatch");
    }
    for(const auto& a:v.at("rkhs_adjoints").array()) {
        RkhsAdjointSpectrum k;k.dimension=size(a.at("dimension"));k.nodesCount=size(a.at("nodes_count"));k.rank=size(a.at("rank"));
        k.nodes=numbers(a.at("nodes"));k.successors=numbers(a.at("successors"));k.bandwidth=a.at("bandwidth").number();
        k.eigenvalues=complexes(a.at("eigenvalues"));k.sectionCoefficients=complexes(a.at("section_coefficients"));k.assumptions=a.at("assumptions").string();
        k.residuals=numbers(a.at("residuals"));k.pseudospectrum.grid=complexes(a.at("pseudospectrum_grid"));
        k.pseudospectrum.minimum=numbers(a.at("pseudospectrum_minimum"));k.pseudospectrumMinimizers=complexes(a.at("pseudospectrum_minimizers"));
        k.pseudospectrum.validate();
        if(!k.dimension||!k.rank||k.rank>k.nodesCount||k.nodes.size()!=k.nodesCount*k.dimension||k.successors.size()!=k.nodes.size()||
            k.eigenvalues.size()!=k.rank||k.sectionCoefficients.size()!=k.nodesCount*k.rank||
            k.pseudospectrumMinimizers.size()!=k.pseudospectrum.grid.size()*k.nodesCount)throw std::invalid_argument("RKHS section spectrum shape mismatch");
        s.rkhsAdjoints.push_back(std::move(k));
    }
    const auto& i=v.at("invariance");s.invariance.principalAngles=numbers(i.at("principal_angles"));
    s.invariance.projectionDefect=i.at("projection_defect").number();s.invariance.sourceRank=size(i.at("source_rank"));s.invariance.successorRank=size(i.at("successor_rank"));
    return s;
}
json::Value koopmanCompressionJson(const KoopmanCompression& m) {
    m.validate();O readouts;for(const auto& r:m.readouts)readouts[r.first]=O{{"outputs",double(r.second.outputs)},
        {"coefficients",numbers(r.second.coefficients)},{"relative_training_error",r.second.relativeTrainingError}};
    return O{{"schema","dynsyslab.koopman-compression/v1"},{"dictionary",spectralDictionaryJson(m.dictionary)},
        {"method",m.method},{"assumptions",m.assumptions},{"coordinate_contract",m.coordinateContract},
        {"sample_interval",m.sampleInterval},{"rank",double(m.rank)},{"samples",double(m.samples)},
        {"basis",numbers(m.basis)},{"action",numbers(m.action)},{"readouts",readouts}};
}
KoopmanCompression koopmanCompressionFromJson(const V& v) {
    if(v.at("schema").string()!="dynsyslab.koopman-compression/v1")throw std::invalid_argument("Unknown Koopman compression schema");
    KoopmanCompression m;m.dictionary=spectralDictionaryFromJson(v.at("dictionary"));m.method=v.at("method").string();
    m.assumptions=v.at("assumptions").string();m.coordinateContract=v.at("coordinate_contract").string();
    m.sampleInterval=v.at("sample_interval").number();m.rank=size(v.at("rank"));m.samples=size(v.at("samples"));
    m.basis=numbers(v.at("basis"));m.action=numbers(v.at("action"));
    for(const auto& r:v.at("readouts").object())m.readouts[r.first]={size(r.second.at("outputs")),
        numbers(r.second.at("coefficients")),r.second.at("relative_training_error").number()};
    m.validate();return m;
}
std::string koopmanCompressionIdentity(const KoopmanCompression& model) {
    // Dictionary, operator and readouts: changing a referenced decomposition changes the identity.
    return sha256Text(json::dump(koopmanCompressionJson(model)));
}
} // namespace dynsyslab::core
