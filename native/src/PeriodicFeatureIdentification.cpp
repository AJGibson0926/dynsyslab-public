#include <dynsyslab/core/PeriodicFeatureIdentification.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
// SVD whitening / generalized least squares: Golub & Van Loan, Matrix
// Computations. Bounded nonlinear fitting: Nocedal & Wright, Numerical
// Optimization, Ch. 10. This is atlas identification, not a linear KROM.
namespace dynsyslab::core {
PeriodicFeatureIdentificationBatch identifyPeriodicHistoryFeatures(
    const math::PeriodicAtlas& atlas,const PeriodicHistoryFeatureMap& map,
    const std::vector<double>& histories,const std::vector<double>& times,
    const PeriodicFeatureIdentificationOptions& options){
    const auto nc=atlas.knots().front().channels.size(),nd=map.timeOffsets.size(),nm=map.channels.size();
    if(!nd||!nm||nd>4096/nm||!map.outputs||map.outputs>4096)throw std::invalid_argument("Bounded nonempty history/feature dimensions required");
    const auto n=nd*nm,m=map.outputs;
    if(map.coefficients.size()!=n*m||map.noiseStandardDeviations.size()!=n||histories.empty()||histories.size()%n||
       histories.size()>10000000||atlas.knots().size()<2)throw std::invalid_argument("History/feature/noise/atlas dimensions mismatch");
    if(!std::isfinite(map.relativeTolerance)||map.relativeTolerance<=0||map.relativeTolerance>=1)throw std::invalid_argument("Feature tolerance must lie in (0,1)");
    for(auto c:map.channels)if(c>=nc)throw std::invalid_argument("Feature channel outside atlas");
    for(double t:map.timeOffsets)if(!std::isfinite(t)||t>0)throw std::invalid_argument("Feature histories must be past-only");
    if(*std::max_element(map.timeOffsets.begin(),map.timeOffsets.end())!=0)throw std::invalid_argument("History must end at zero");
    for(double v:map.noiseStandardDeviations)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Positive independent input noise scales required");
    for(const auto* a:{&map.coefficients,&histories,&times})for(double v:*a)if(!std::isfinite(v))throw std::invalid_argument("Finite feature inputs required");
    for(double t:times)if(t<0)throw std::invalid_argument("Forecast times must be nonnegative");
    if(options.phaseSeeds<4||options.phaseSeeds>16384||!options.starts||options.starts>64||
       !options.parameterSubdivisions||options.parameterSubdivisions>64||
       atlas.knots().size()>1000000/options.parameterSubdivisions)throw std::invalid_argument("Invalid coarse seed budget");
    const auto np=(atlas.knots().size()-1)*options.parameterSubdivisions+1;
    if(np>1000000/options.phaseSeeds)throw std::invalid_argument("Too many coarse seeds");
    const auto ns=np*options.phaseSeeds,nb=histories.size()/n;
    if(ns>50000000/n||ns>50000000/m||nb>10000000/nc||(!times.empty()&&nb>10000000/nc/times.size()))throw std::invalid_argument("Feature workspace/output budget exceeded");
    PeriodicFeatureIdentificationBatch out;out.inputs=n;out.outputs=m;
    auto scaled=map.coefficients;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<m;++j)scaled[i*m+j]*=map.noiseStandardDeviations[i];
    for(double v:scaled)if(!std::isfinite(v))throw std::invalid_argument("Feature/noise product overflow");
    const auto svd=numerics::singularValueDecomposition(scaled,n,m,0,map.relativeTolerance);
    if(!svd.rank)throw std::invalid_argument("Feature map has zero retained noise rank");
    out.rank=svd.rank;out.singularValues=svd.singularValues;
    const auto nr=map.whiten?svd.rank:m;out.effectiveCoefficients=map.coefficients;
    if(nb>50000000/nr)throw std::invalid_argument("Transformed observations exceed workspace budget");
    if(map.whiten){out.effectiveCoefficients.resize(n*nr);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<nr;++j)
            out.effectiveCoefficients[i*nr+j]=svd.left[i*nr+j]/map.noiseStandardDeviations[i];}
    const auto& transform=out.effectiveCoefficients;const math::PeriodicAtlasTable table(atlas,options.tableSamples);
    const double pi=std::acos(-1.),lower=atlas.knots().front().parameter,upper=atlas.knots().back().parameter;
    const auto real=[](std::complex<double> v){if(std::abs(v.imag())>1e-9*(1+std::abs(v.real())))throw std::invalid_argument("Real history/output channels required");return v.real();};
    const auto history=[&](double p,double phase){std::vector<double> h;h.reserve(n);const auto w=table.frequency(p);
        for(double t:map.timeOffsets){const auto v=table.evaluate(p,phase+w*t);for(auto c:map.channels)h.push_back(real(v[c]));}return h;};
    const auto features=[&](const std::vector<double>& h){return numerics::multiply(h,1,n,transform,nr);};
    std::vector<double> seeds,grid;seeds.reserve(ns*2);grid.reserve(ns*n);
    for(std::size_t k=0;k<np;++k){const auto left=std::min(atlas.knots().size()-2,k/options.parameterSubdivisions);
        const double fraction=double(k-left*options.parameterSubdivisions)/options.parameterSubdivisions;
        const double p=(1-fraction)*atlas.knots()[left].parameter+fraction*atlas.knots()[left+1].parameter;
        for(std::size_t j=0;j<options.phaseSeeds;++j){const double phase=2*pi*j/options.phaseSeeds;
            seeds.push_back(p);seeds.push_back(phase);const auto h=history(p,phase);grid.insert(grid.end(),h.begin(),h.end());}}
    const auto candidate=numerics::multiply(grid,ns,n,transform,nr),observed=numerics::multiply(histories,nb,n,transform,nr);
    for(std::size_t b=0;b<nb;++b){
        std::vector<std::pair<double,std::size_t>> costs;costs.reserve(ns);
        for(std::size_t s=0;s<ns;++s){double cost=0;for(std::size_t j=0;j<nr;++j)cost+=std::pow(candidate[s*nr+j]-observed[b*nr+j],2);costs.push_back({cost,s});}
        std::sort(costs.begin(),costs.end());
        std::vector<std::size_t> selected;
        for(auto [cost,s]:costs){bool separated=true;
            for(auto k:selected)if(std::abs(seeds[2*s]-seeds[2*k])/(upper-lower)<.02&&
                std::abs(std::remainder(seeds[2*s+1]-seeds[2*k+1],2*pi))<2*pi/32)separated=false;
            if(separated)selected.push_back(s);if(selected.size()==options.starts)break;}
        const auto residual=[&](const std::vector<double>& p){auto v=features(history(p[0],p[1]));
            for(std::size_t j=0;j<nr;++j)v[j]-=observed[b*nr+j];return v;};
        PeriodicFeatureIdentification result;result.fit.squaredResidual=std::numeric_limits<double>::infinity();
        for(auto s:selected){const auto phase=seeds[2*s+1];auto fit=numerics::boundedNonlinearLeastSquares(residual,
            {seeds[2*s],phase},{lower,phase-pi},{upper,phase+pi},options.optimization);
            if(fit.squaredResidual<result.fit.squaredResidual)result.fit=std::move(fit);}
        result.fit.parameters[1]=std::remainder(result.fit.parameters[1],2*pi);
        const auto p=result.fit.parameters[0],phase=result.fit.parameters[1];result.frequency=table.frequency(p);
        const auto h=history(p,phase);double err=0,norm=0;
        for(std::size_t i=0;i<n;++i){err+=std::pow(h[i]-histories[b*n+i],2);norm+=std::pow(histories[b*n+i],2);}
        result.historyRelativeError=std::sqrt(err/std::max(norm,1e-300));
        for(auto v:table.evaluate(p,phase))result.current.push_back(real(v));
        for(double t:times)for(auto v:table.evaluate(p,phase+result.frequency*t))result.forecast.push_back(real(v));
        std::vector<double> jac(nr*2);const double eps=options.optimization.finiteDifferenceStep;
        const double a=std::max(lower,p-eps*(upper-lower)),z=std::min(upper,p+eps*(upper-lower));
        const auto fa=features(history(a,phase)),fz=features(history(z,phase));
        const auto fm=features(history(p,phase-eps*2*pi)),fp=features(history(p,phase+eps*2*pi));
        for(std::size_t j=0;j<nr;++j){jac[2*j]=(fz[j]-fa[j])/(z-a)*(upper-lower);jac[2*j+1]=(fp[j]-fm[j])/(2*eps);}
        result.jacobianSingularValues=numerics::singularValueDecomposition(jac,nr,2).singularValues;
        out.estimates.push_back(std::move(result));
    }return out;
}
}
