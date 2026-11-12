#include <dynsyslab/core/SpectralKROM.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
void SpectralDictionary::validate() const {
    if(!dimension)throw std::invalid_argument("Dictionary dimension must be positive");
    if(kind!="linear"&&kind!="polynomial"&&kind!="gaussian"&&kind!="partition")
        throw std::invalid_argument("Unknown spectral dictionary");
    if((!center.empty()&&center.size()!=dimension)||(!scale.empty()&&scale.size()!=dimension))
        throw std::invalid_argument("Dictionary normalization shape mismatch");
    for(double v:center)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary center");
    for(double v:scale)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Dictionary scales must be positive");
    for(double v:nodes)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary node");
    if((kind=="gaussian"||kind=="partition")&&(nodes.empty()||nodes.size()%dimension))
        throw std::invalid_argument("Dictionary nodes must be centers x dimension");
    if(kind=="gaussian"&&(!std::isfinite(bandwidth)||bandwidth<=0))
        throw std::invalid_argument("Gaussian bandwidth must be positive");
    if(kind=="polynomial"&&(exponents.empty()||exponents.size()%dimension))
        throw std::invalid_argument("Polynomial exponents must be features x dimension");
}
std::size_t SpectralDictionary::size() const {
    if(kind=="linear")return dimension+std::size_t(constant);
    if(kind=="polynomial")return dimension?exponents.size()/dimension:0;
    return dimension?nodes.size()/dimension:0;
}
std::vector<double> SpectralDictionary::evaluate(const std::vector<double>& points) const {
    validate();
    if(points.empty()||points.size()%dimension)throw std::invalid_argument("Dictionary points shape mismatch");
    const auto m=points.size()/dimension,n=size();
    std::vector<double> out(m*n),q(dimension);
    for(std::size_t i=0;i<m;++i) {
        for(std::size_t d=0;d<dimension;++d) {
            if(!std::isfinite(points[i*dimension+d]))throw std::invalid_argument("Nonfinite dictionary input");
            q[d]=(points[i*dimension+d]-(center.empty()?0:center[d]))/(scale.empty()?1:scale[d]);
        }
        if(kind=="linear") {
            if(constant)out[i*n]=1;
            for(std::size_t d=0;d<dimension;++d)out[i*n+d+std::size_t(constant)]=q[d];
        } else if(kind=="polynomial") {
            for(std::size_t j=0;j<n;++j) {
                double v=1;
                for(std::size_t d=0;d<dimension;++d)v*=std::pow(q[d],double(exponents[j*dimension+d]));
                out[i*n+j]=v;
            }
        } else {
            double nearest=std::numeric_limits<double>::infinity();std::size_t best=0;
            for(std::size_t j=0;j<n;++j) {
                double r2=0;
                for(std::size_t d=0;d<dimension;++d) {
                    const double node=(nodes[j*dimension+d]-(center.empty()?0:center[d]))/(scale.empty()?1:scale[d]);
                    r2+=(q[d]-node)*(q[d]-node);
                }
                if(kind=="gaussian")out[i*n+j]=std::exp(-r2/(2*bandwidth*bandwidth));
                if(r2<nearest){nearest=r2;best=j;}
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
} // namespace dynsyslab::core
