// Measurement-only composition of the existing native least-squares and map APIs.
// Same sampled residual, column scaling and ridge as HistoryInputRegression;
// predictors now include EVERY selected eigenfunction for EVERY output row.
#include "echo_multimode.hpp"
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace echo_multi {
namespace {
std::vector<double> features(const ModeSamples& rows,const std::vector<std::size_t>& ids,std::size_t k) {
    std::vector<double> f{1.};
    for(auto j:ids){const auto z=rows.at(j).at(k).psi;f.push_back(z.real());f.push_back(z.imag());}
    return f;
}
}
JointInputFit fitJointInput(const ModeSamples& rows,const std::vector<std::size_t>& ids,
    const std::vector<Complex>& eigenvalues,double step,double ridge) {
    if(ids.empty()||!(step>0)||ridge<0)throw std::invalid_argument("invalid joint fit configuration");
    const auto n=rows.at(ids[0]).size(),m=ids.size(),p=1+2*m;
    if(!n)throw std::invalid_argument("empty joint fit");
    for(auto j:ids)if(rows.at(j).size()!=n)throw std::invalid_argument("unaligned eigenfunction rows");
    std::vector<double> design(n*p),rhs(n*2*m),scales(p);
    for(std::size_t k=0;k<n;++k) {
        const auto f=features(rows,ids,k);const double u=rows[ids[0]][k].input;
        for(std::size_t c=0;c<p;++c){design[k*p+c]=u*f[c];scales[c]+=design[k*p+c]*design[k*p+c];}
        for(std::size_t j=0;j<m;++j) {
            const auto& s=rows[ids[j]][k];
            if(s.input!=u)throw std::invalid_argument("joint fit has misaligned inputs");
            const auto q=s.nextPsi-std::exp(eigenvalues.at(ids[j])*step)*s.psi;
            rhs[k*2*m+2*j]=q.real();rhs[k*2*m+2*j+1]=q.imag();
        }
    }
    for(std::size_t c=0;c<p;++c) {
        scales[c]=std::sqrt(scales[c]);
        if(!(scales[c]>1e-20)||!std::isfinite(scales[c]))throw std::invalid_argument("unexcited joint predictor");
        for(std::size_t k=0;k<n;++k)design[k*p+c]/=scales[c];
    }
    dynsyslab::numerics::LeastSquaresOptions opts;
    opts.method=dynsyslab::numerics::LeastSquaresMethod::divideAndConquerSvd;opts.ridge=ridge;
    const auto coefficients=dynsyslab::numerics::solveLeastSquares(design,n,p,rhs,2*m,opts);
    const auto singular=dynsyslab::numerics::singularValues(std::vector<Complex>(design.begin(),design.end()),n,p);
    JointInputFit result;result.columns=p;result.samples=n;
    result.rank=std::count_if(singular.begin(),singular.end(),[&](double s){return s>singular.front()*1e-10;});
    result.condition=singular.front()/std::max(singular.back(),1e-300);
    result.coefficients.resize(m*p);result.map.modes=m;result.map.inputs=1;
    result.map.offsets.resize(m);result.map.operators.resize(m*m);result.map.conjugateOperators.resize(m*m);
    result.map.provenance="joint subset-specific echo-history sampled residual; all selected modes and conjugates; ridge=1e-8 by default";
    const Complex i(0,1);
    for(std::size_t j=0;j<m;++j) {
        for(std::size_t c=0;c<p;++c)result.coefficients[j*p+c]=Complex(coefficients[c*2*m+2*j],coefficients[c*2*m+2*j+1])/scales[c];
        const auto l=eigenvalues[ids[j]],integral=std::abs(l)>1e-14?(std::exp(l*step)-1.)/l:Complex(step,0);
        result.map.offsets[j]=result.coefficients[j*p]/integral;
        for(std::size_t a=0;a<m;++a) {
            const auto re=result.coefficients[j*p+1+2*a],im=result.coefficients[j*p+2+2*a];
            result.map.operators[j*m+a]=(re-i*im)/(2.*integral);
            result.map.conjugateOperators[j*m+a]=(re+i*im)/(2.*integral);
        }
    }
    // Exact scalar-method compatibility: same data, feature ordering and ridge.
    if(m==1) {
        const auto scalar=core::fitHistoryInputMap(rows[ids[0]],std::exp(eigenvalues[ids[0]]*step),step,1,{1,1},ridge);
        for(std::size_t c=0;c<p;++c)
            if(std::abs(result.coefficients[c]-scalar.coefficients[c])>1e-10*std::max(1.,std::abs(scalar.coefficients[c])))
                throw std::runtime_error("joint m=1 differs from canonical history-input fitter");
    }
    return result;
}
Json scoreJointInput(const JointInputFit& fit,const ModeSamples& rows,const std::vector<std::size_t>& ids,
    const std::vector<Complex>& eigenvalues,double step) {
    const auto m=ids.size(),p=fit.columns,n=rows.at(ids[0]).size();
    if(!n)throw std::invalid_argument("empty joint score");
    std::vector<double> error(m),energy(m);double canonicalError=0;
    for(std::size_t k=0;k<n;++k) {
        const auto f=features(rows,ids,k);
        for(std::size_t j=0;j<m;++j) {
            const auto& s=rows[ids[j]][k];const auto l=eigenvalues[ids[j]],mu=std::exp(l*step);
            const auto q=s.nextPsi-mu*s.psi;Complex b{};
            for(std::size_t c=0;c<p;++c)b+=fit.coefficients[j*p+c]*f[c];
            error[j]+=std::norm(q-b*s.input);energy[j]+=std::norm(q);
            Complex bc=fit.map.offsets[j];
            for(std::size_t a=0;a<m;++a){const auto z=rows[ids[a]][k].psi;
                bc+=fit.map.operators[j*m+a]*z+fit.map.conjugateOperators[j*m+a]*std::conj(z);}
            const auto integral=std::abs(l)>1e-14?(mu-1.)/l:Complex(step,0);
            canonicalError=std::max(canonicalError,std::abs(b-bc*integral)/std::max(1.,std::abs(b)));
        }
    }
    if(canonicalError>1e-8)throw std::runtime_error("joint canonical conjugate map parity failed");
    Json::Array modes;double mean=0;
    for(std::size_t j=0;j<m;++j){const double e=std::sqrt(error[j]/std::max(energy[j],1e-300));mean+=e;
        modes.emplace_back(Json::Object{{"mode",double(ids[j]+1)},{"relative_residual",e}});}
    return Json::Object{{"mean_mode_relative_residual",mean/m},{"per_mode",modes},
        {"samples",double(n)},{"canonical_map_max_relative_discrepancy",canonicalError}};
}
Json jointInputTests() {
    const std::vector<Complex> lambdas{{-.1,6.},{-.2,12.},{-.3,18.}};
    ModeSamples rows(3);std::vector<Complex> truth(3*7);
    for(std::size_t j=0;j<3;++j)for(std::size_t c=0;c<7;++c)truth[j*7+c]=Complex(.01*(1+j+c),-.02*(1+j)+.003*c);
    for(std::size_t k=0;k<400;++k) {
        const double t=k*.031,u=.2+.3*std::sin(1.31*t);std::vector<Complex> z;std::vector<double> f{1.};
        for(std::size_t j=0;j<3;++j){z.emplace_back(std::sin((j+1.1)*t),std::cos((j+1.7)*t));f.push_back(z.back().real());f.push_back(z.back().imag());}
        for(std::size_t j=0;j<3;++j){Complex b{};for(std::size_t c=0;c<7;++c)b+=truth[j*7+c]*f[c];
            rows[j].push_back({{},z[j],std::exp(lambdas[j]*.01)*z[j]+b*u,u});}
    }
    const auto fit=fitJointInput(rows,{0,1,2},lambdas,.01,0.);
    double maximum=0;for(std::size_t i=0;i<truth.size();++i)maximum=std::max(maximum,std::abs(fit.coefficients[i]-truth[i]));
    const auto perm=fitJointInput(rows,{2,0,1},lambdas,.01,0.);
    const auto s=scoreJointInput(fit,rows,{0,1,2},lambdas),p=scoreJointInput(perm,rows,{2,0,1},lambdas);
    if(maximum>1e-10||p.at("mean_mode_relative_residual").number()>1e-10)throw std::runtime_error("joint cross-mode recovery test failed");
    // Also exercises the explicit scalar compatibility assertion.
    (void)fitJointInput(rows,{0},lambdas);
    return Json::Object{{"synthetic_dense_cross_mode_coefficient_error",maximum},{"synthetic_score",s},
        {"permuted_score",p},{"single_mode_parity","asserted against existing native fitter on synthetic and real records"}};
}
}
