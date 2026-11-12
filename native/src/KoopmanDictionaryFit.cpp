#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
SpectralDictionary buildSpectralDictionary(const std::vector<double>& points,std::size_t dimension,
                                           const SpectralDictionaryOptions& options) {
    if(!dimension||points.empty()||points.size()%dimension||!options.maximumFeatures)
        throw std::invalid_argument("Invalid dictionary training shape or feature limit");
    for(double v:points)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary training data");
    const auto m=points.size()/dimension;
    SpectralDictionary d;d.dimension=dimension;d.kind=options.kind;d.constant=options.constant;d.bandwidth=options.bandwidth;
    if(options.normalize) {
        d.center.resize(dimension);d.scale.resize(dimension);
        for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<dimension;++j)d.center[j]+=points[i*dimension+j]/m;
        for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<dimension;++j)d.scale[j]+=std::pow(points[i*dimension+j]-d.center[j],2)/m;
        for(auto& s:d.scale)s=std::sqrt(s)>1e-14?std::sqrt(s):1.;
    }
    if(d.kind=="polynomial") {
        if(options.degree>32)throw std::invalid_argument("Polynomial degree exceeds supported limit 32");
        std::vector<std::size_t> exponent(dimension);
        std::function<void(std::size_t,std::size_t)> enumerate=[&](std::size_t axis,std::size_t left) {
            if(axis+1==dimension) {
                exponent[axis]=left;
                if(d.exponents.size()/dimension>=options.maximumFeatures)throw std::invalid_argument("Polynomial dictionary exceeds maximum feature count");
                d.exponents.insert(d.exponents.end(),exponent.begin(),exponent.end());return;
            }
            for(std::size_t power=0;power<=left;++power){exponent[axis]=power;enumerate(axis+1,left-power);}
        };
        for(std::size_t degree=options.constant?0:1;degree<=options.degree;++degree)enumerate(0,degree);
    } else if(d.kind=="gaussian"||d.kind=="partition") {
        if(!options.centers||options.centers>options.maximumFeatures)throw std::invalid_argument("Invalid dictionary center count");
        // Deterministic farthest-point centers are selected using training coordinates only.
        std::vector<double> nearest(m,std::numeric_limits<double>::infinity());std::size_t index=0;
        for(std::size_t j=0;j<std::min(options.centers,m);++j) {
            d.nodes.insert(d.nodes.end(),points.begin()+index*dimension,points.begin()+(index+1)*dimension);
            for(std::size_t i=0;i<m;++i) {
                double r2=0;for(std::size_t k=0;k<dimension;++k)r2+=std::pow((points[i*dimension+k]-points[index*dimension+k])/(d.scale.empty()?1:d.scale[k]),2);
                nearest[i]=std::min(nearest[i],r2);
            }
            index=std::size_t(std::max_element(nearest.begin(),nearest.end())-nearest.begin());
            if(nearest[index]<=1e-24)break;
        }
    }
    d.validate();if(d.size()>options.maximumFeatures)throw std::invalid_argument("Dictionary exceeds maximum feature count");
    return d;
}
}
