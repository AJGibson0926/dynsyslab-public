#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/BasisFunctions.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
void productIndices(SpectralDictionary& d,const SpectralDictionaryOptions& o) {
    std::vector<std::size_t> index(d.dimension);
    auto append=[&]() {
        if(!o.constant&&std::all_of(index.begin(),index.end(),[](auto x){return x==0;}))return;
        if(d.exponents.size()/d.dimension>=o.maximumFeatures)
            throw std::invalid_argument("Dictionary exceeds maximum feature count");
        d.exponents.insert(d.exponents.end(),index.begin(),index.end());
    };
    if(d.kind=="fourier_hermite") {
        if(d.dimension!=2||o.angularOrder>512||o.hermiteOrder>512)
            throw std::invalid_argument("Fourier-Hermite requires dimension 2 and orders <=512");
        for(std::size_t k=0;k<=2*o.angularOrder;++k)for(std::size_t n=0;n<=o.hermiteOrder;++n) {
            index={k,n};append();
        }
    } else {
        if(o.degree>32)throw std::invalid_argument("Product degree exceeds supported limit 32");
        if(o.indexSet=="total_degree") {
            std::function<void(std::size_t,std::size_t)> enumerate=[&](std::size_t axis,std::size_t left) {
                if(axis+1==d.dimension){index[axis]=left;append();return;}
                for(std::size_t p=0;p<=left;++p){index[axis]=p;enumerate(axis+1,left-p);}
            };
            for(std::size_t degree=0;degree<=o.degree;++degree)enumerate(0,degree);
        } else if(o.indexSet=="tensor_product") {
            std::function<void(std::size_t)> enumerate=[&](std::size_t axis) {
                if(axis==d.dimension){append();return;}
                for(std::size_t p=0;p<=o.degree;++p){index[axis]=p;enumerate(axis+1);}
            };
            enumerate(0);
        } else throw std::invalid_argument("index_set must be total_degree or tensor_product");
    }
}
void centers(SpectralDictionary& d,const std::vector<double>& points,const SpectralDictionaryOptions& o) {
    if(!o.centers||o.centers>o.maximumFeatures)throw std::invalid_argument("Invalid dictionary center count");
    if(o.centerSelection!="farthest"&&o.centerSelection!="kmeans")
        throw std::invalid_argument("center_selection must be farthest or kmeans");
    if(o.centerSelection=="kmeans"&&(!o.centerIterations||o.centerIterations>1000))
        throw std::invalid_argument("center_iterations must be in [1,1000]");
    const auto m=points.size()/d.dimension;
    auto distance=[&](const double* a,const double* b) {
        double result=0;
        for(std::size_t k=0;k<d.dimension;++k)
            result=std::hypot(result,(a[k]-b[k])/(d.scale.empty()?1:d.scale[k]));
        return result;
    };
    // Deterministic initialization; use training X only, never held-out points.
    std::vector<double> nearest(m,std::numeric_limits<double>::infinity());std::size_t index=0;
    for(std::size_t j=0;j<std::min(o.centers,m);++j) {
        d.nodes.insert(d.nodes.end(),points.begin()+index*d.dimension,points.begin()+(index+1)*d.dimension);
        for(std::size_t i=0;i<m;++i)nearest[i]=std::min(nearest[i],distance(&points[i*d.dimension],&points[index*d.dimension]));
        index=std::size_t(std::max_element(nearest.begin(),nearest.end())-nearest.begin());
        if(nearest[index]<=1e-12)break;
    }
    if(o.centerSelection!="kmeans")return;
    const auto n=d.nodes.size()/d.dimension;
    for(std::size_t iteration=0;iteration<o.centerIterations;++iteration) {
        std::vector<double> means(d.nodes.size());std::vector<std::size_t> counts(n);
        for(std::size_t i=0;i<m;++i) {
            std::size_t best=0;double shortest=std::numeric_limits<double>::infinity();
            for(std::size_t j=0;j<n;++j) {
                const auto r=distance(&points[i*d.dimension],&d.nodes[j*d.dimension]);
                if(r<shortest){shortest=r;best=j;}
            }
            ++counts[best];
            for(std::size_t k=0;k<d.dimension;++k) {
                auto& mean=means[best*d.dimension+k];mean+=(points[i*d.dimension+k]-mean)/counts[best];
            }
        }
        double movement=0;
        for(std::size_t j=0;j<n;++j) {
            if(!counts[j])std::copy_n(&d.nodes[j*d.dimension],d.dimension,&means[j*d.dimension]);
            movement=std::max(movement,distance(&means[j*d.dimension],&d.nodes[j*d.dimension]));
        }
        d.nodes=std::move(means);if(movement<=1e-12)break;
    }
}
}
SpectralDictionary buildSpectralDictionary(const std::vector<double>& points,std::size_t dimension,
                                           const SpectralDictionaryOptions& o) {
    if(!dimension||points.empty()||points.size()%dimension||!o.maximumFeatures)
        throw std::invalid_argument("Invalid dictionary training shape or feature limit");
    for(double v:points)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary training data");
    if((!o.center.empty()&&o.center.size()!=dimension)||(!o.scale.empty()&&o.scale.size()!=dimension))
        throw std::invalid_argument("Dictionary normalization shape mismatch");
    if(o.normalize&&(!o.center.empty()||!o.scale.empty()))
        throw std::invalid_argument("Choose training normalization or explicit center/scale");
    if(o.kind=="fourier_hermite"&&o.normalize)
        throw std::invalid_argument("Do not normalize the periodic angle; set explicit velocity scale instead");
    SpectralDictionary d;d.dimension=dimension;d.kind=o.kind;d.constant=o.constant;
    d.bandwidth=o.bandwidth;d.smoothness=o.smoothness;d.amplitude=o.amplitude;d.center=o.center;d.scale=o.scale;
    for(double v:d.center)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite dictionary center");
    for(double v:d.scale)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Dictionary scales must be positive");
    const auto m=points.size()/dimension;
    if(o.normalize) {
        d.center.resize(dimension);d.scale.resize(dimension);
        for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<dimension;++j)d.center[j]+=points[i*dimension+j]/m;
        for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<dimension;++j)d.scale[j]+=std::pow(points[i*dimension+j]-d.center[j],2)/m;
        for(auto& s:d.scale)s=std::sqrt(s)>1e-14?std::sqrt(s):1.;
    }
    if(basis::isProduct(d.kind))productIndices(d,o);
    else if(basis::isRadial(d.kind)||d.kind=="partition")centers(d,points,o);
    d.validate();if(d.size()>o.maximumFeatures)throw std::invalid_argument("Dictionary exceeds maximum feature count");
    return d;
}
}
