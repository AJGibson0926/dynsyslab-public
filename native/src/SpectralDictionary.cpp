#include <dynsyslab/core/SpectralKROM.hpp>
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <algorithm>
#include <set>
#include <dynsyslab/core/BasisFunctions.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
void SpectralDictionary::validate() const {
    if(!dimension)throw std::invalid_argument("Dictionary dimension must be positive");
    const bool product=basis::isProduct(kind),radial=basis::isRadial(kind);
    if(kind!="linear"&&!product&&!radial&&kind!="partition")throw std::invalid_argument("Unknown spectral dictionary");
    if((!center.empty()&&center.size()!=dimension)||(!scale.empty()&&scale.size()!=dimension))
        throw std::invalid_argument("Dictionary normalization shape mismatch");
    for(double v:center)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary center");
    for(double v:scale)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Dictionary scales must be positive");
    for(double v:nodes)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary node");
    if((radial||kind=="partition")&&(nodes.empty()||nodes.size()%dimension))
        throw std::invalid_argument("Dictionary nodes must be centers x dimension");
    if(radial)basis::radialKernel(kind,0,bandwidth,smoothness,amplitude);
    if(product) {
        if(exponents.empty()||exponents.size()%dimension)throw std::invalid_argument("Basis indices must be features x dimension");
        std::set<std::vector<std::size_t>> seen;
        for(std::size_t j=0;j<exponents.size();j+=dimension) {
            if(!seen.emplace(exponents.begin()+j,exponents.begin()+j+dimension).second)
                throw std::invalid_argument("Duplicate dictionary multi-index");
            for(std::size_t k=0;k<dimension;++k)
                if(exponents[j+k]>(kind=="fourier_hermite"&&k==0?1024:512))
                    throw std::invalid_argument("Dictionary basis order exceeds limit");
        }
    }
    if(kind=="fourier_hermite"&&(dimension!=2||(!scale.empty()&&scale[0]!=1)))
        throw std::invalid_argument("Fourier-Hermite requires (angle,velocity), with unscaled radian angle");
}
std::size_t SpectralDictionary::size() const {
    if(kind=="linear")return dimension+std::size_t(constant);
    if(basis::isProduct(kind))return dimension?exponents.size()/dimension:0;
    return dimension?nodes.size()/dimension:0;
}
std::vector<double> SpectralDictionary::evaluate(const std::vector<double>& points) const {
    validate();
    if(points.empty()||points.size()%dimension)throw std::invalid_argument("Dictionary points shape mismatch");
    const auto m=points.size()/dimension,n=size();
    std::vector<double> out(m*n),q(dimension);
    const double pi=3.14159265358979323846;
    for(std::size_t i=0;i<m;++i) {
        for(std::size_t d=0;d<dimension;++d) {
            if(!std::isfinite(points[i*dimension+d]))throw std::invalid_argument("Nonfinite dictionary input");
            q[d]=(points[i*dimension+d]-(center.empty()?0:center[d]))/(scale.empty()?1:scale[d]);
        }
        if(kind=="linear") {
            if(constant)out[i*n]=1;
            for(std::size_t d=0;d<dimension;++d)out[i*n+d+std::size_t(constant)]=q[d];
        } else if(basis::isProduct(kind)) {
            for(std::size_t j=0;j<n;++j) {
                double v=1;
                for(std::size_t d=0;d<dimension;++d) {
                    const auto order=exponents[j*dimension+d];
                    if(kind=="polynomial")v*=std::pow(q[d],double(order));
                    else if(kind=="chebyshev")v*=basis::chebyshev(order,q[d]);
                    else if(kind=="laguerre")v*=basis::laguerre(order,q[d]);
                    else if(kind=="hermite"||d==1)v*=basis::hermiteFunction(order,q[d]);
                    else if(order==0)v/=std::sqrt(2*pi);
                    else {
                        const auto k=(order+1)/2;
                        v*=(order%2?std::cos(k*q[d]):std::sin(k*q[d]))/std::sqrt(pi);
                    }
                }
                out[i*n+j]=v;
            }
        } else {
            double nearest=std::numeric_limits<double>::infinity();std::size_t best=0;
            for(std::size_t j=0;j<n;++j) {
                double distance=0;
                for(std::size_t d=0;d<dimension;++d) {
                    const double node=(nodes[j*dimension+d]-(center.empty()?0:center[d]))/(scale.empty()?1:scale[d]);
                    distance=std::hypot(distance,q[d]-node);
                }
                if(basis::isRadial(kind))out[i*n+j]=basis::radialKernel(kind,distance,bandwidth,smoothness,amplitude);
                if(distance<nearest){nearest=distance;best=j;}
            }
            if(kind=="partition")out[i*n+best]=1;
        }
    }
    for(double v:out)if(!std::isfinite(v))throw std::invalid_argument("Dictionary evaluation overflow");
    return out;
}
void SpectralKROM::validate() const {
    dictionary.validate(); const auto n=dictionary.size();
    if(!rank||rank>n||!samples||!std::isfinite(sampleInterval)||sampleInterval<=0||
        basis.size()!=n*rank||action.size()!=rank*rank||multipliers.size()!=rank||
        eigenfunctions.size()!=n*rank)throw std::invalid_argument("Invalid spectral KROM dimensions or clock");
    for(const auto& r:reconstructions)if(!r.second.outputs||r.second.coefficients.size()!=n*r.second.outputs)
        throw std::invalid_argument("Invalid spectral reconstruction shape");
    for(const auto* values:{&action,&basis})for(double v:*values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite spectral operator/basis");
    for(const auto* values:{&multipliers,&eigenfunctions})for(auto v:*values)
        if(!std::isfinite(v.real())||!std::isfinite(v.imag()))throw std::invalid_argument("Nonfinite spectral eigenpair");
    for(const auto& item:reconstructions)for(double v:item.second.coefficients)
        if(!std::isfinite(v))throw std::invalid_argument("Nonfinite spectral readout");
    if(!spectral.residuals.eigenpair.empty()&&spectral.residuals.eigenpair.size()!=rank)
        throw std::invalid_argument("Spectral residual count differs from retained modes");
    for(const auto& p:spectral.wavePackets) {
        if(p.features!=n||p.coefficients.size()!=p.grid.size()*n||
           p.dictionaryIdentity!=sha256Text(json::dump(spectralDictionaryJson(dictionary))))
            throw std::invalid_argument("Wave packet dictionary identity or dimensions differ");
    }
}
std::vector<std::complex<double>> SpectralKROM::evaluateEigenfunctions(const std::vector<double>& points) const {
    validate();const auto p=dictionary.evaluate(points);const auto n=dictionary.size(),m=points.size()/dictionary.dimension;
    std::vector<std::complex<double>> out(m*rank);
    for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<rank;++j)
        for(std::size_t k=0;k<n;++k)out[i*rank+j]+=p[i*n+k]*eigenfunctions[k*rank+j];
    return out;
}
std::vector<double> SpectralKROM::reconstruct(const std::string& name,const std::vector<double>& points) const {
    validate();const auto& r=reconstructions.at(name);const auto p=dictionary.evaluate(points);
    const auto n=dictionary.size(),m=points.size()/dictionary.dimension;std::vector<double> out(m*r.outputs);
    for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<r.outputs;++j)
        for(std::size_t k=0;k<n;++k)out[i*r.outputs+j]+=p[i*n+k]*r.coefficients[k*r.outputs+j];
    return out;
}
// The spectral dictionary as dictionary terms (kept here, with SpectralDictionary).
ExpansionDictionary ExpansionDictionary::fromSpectral(const SpectralDictionary& s) {
    s.validate();
    const auto d = s.dimension;
    std::vector<DictionaryTerm> terms;
    auto base = [&](const std::string& kind) { DictionaryTerm t; t.kind = kind; t.center = s.center; t.scale = s.scale; return t; };
    if (s.kind == "linear") {
        if (s.constant) { auto t = base("monomial"); t.orders.assign(d, 0); terms.push_back(t); }
        for (std::size_t i = 0; i < d; ++i) { auto t = base("monomial"); t.orders.assign(d, 0); t.orders[i] = 1; terms.push_back(t); }
    } else if (basis::isProduct(s.kind)) {
        for (std::size_t j = 0; j < s.exponents.size() / d; ++j) {
            auto t = base(s.kind == "polynomial" ? "monomial" : s.kind);
            t.orders.assign(s.exponents.begin() + std::ptrdiff_t(j * d), s.exponents.begin() + std::ptrdiff_t((j + 1) * d));
            terms.push_back(std::move(t));
        }
    } else {
        for (std::size_t j = 0; j < s.nodes.size() / d; ++j) {
            auto t = base(s.kind);
            if (s.kind == "partition") { t.nodes = s.nodes; t.index = j; }
            else {
                t.vector.assign(s.nodes.begin() + std::ptrdiff_t(j * d), s.nodes.begin() + std::ptrdiff_t((j + 1) * d));
                t.width = s.bandwidth; t.smoothness = s.smoothness; t.amplitude = s.amplitude;
            }
            terms.push_back(std::move(t));
        }
    }
    return ExpansionDictionary(d, std::move(terms));
}

} // namespace dynsyslab::core
