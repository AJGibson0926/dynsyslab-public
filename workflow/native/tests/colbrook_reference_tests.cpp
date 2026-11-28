// Frozen-reference and analytic checks for the native mpEDMD and Rigged DMD stages.
// They replace the retired MATLAB scripts check_native_reference.m (pinned upstream rotation),
// compare_rigged_reference.m (native versus upstream on whitened data) and
// rigged_history_reference.m / test_rigged_history_reference.m (weighted block-QR compression).
// Everything runs through the public koopman-spectral query contract
// (dynsyslab.koopman-spectral-query/v1); the model JSON is passed back opaquely, never parsed.
//
// mpEDMD: M. J. Colbrook, "The mpEDMD algorithm for data-driven computations of measure-preserving
//   dynamical systems", SIAM J. Numer. Anal. 61 (2023) 1585-1608, doi:10.1137/22M1521407.
// Rigged DMD: M. J. Colbrook, C. Drysdale and A. Horning, "Rigged Dynamic Mode Decomposition: data-driven
//   generalized eigenfunction decompositions for Koopman operators", arXiv:2405.00782 (2024).
// Block QR of weighted rows (TSQR): J. Demmel, L. Grigori, M. Hoemmen and J. Langou, SIAM J. Sci.
//   Comput. 34 (2012) A206-A239, doi:10.1137/080731992.
// Frozen values: data/colbrook_rotation_reference_v1.json (pinned upstream commits recorded there).
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/numerics/StreamingQr.hpp>
#include <dynsyslab/workflow/KoopmanSpectralQuery.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace json=dynsyslab::core::json;
using J=json::Value;
using cd=std::complex<double>;
constexpr double pi=3.14159265358979323846;

void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
void close(double actual,double expected,double tolerance,const std::string& label){
    if(!std::isfinite(actual)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(label+": actual="+std::to_string(actual)+" expected="+std::to_string(expected));
}
J array(const std::vector<double>& values){J::Array out;for(double v:values)out.emplace_back(v);return out;}
std::vector<double> numbers(const J& value){std::vector<double> out;for(const auto& v:value.array())out.push_back(v.number());return out;}
std::vector<cd> complexes(const J& value){std::vector<cd> out;
    for(const auto& v:value.array())out.emplace_back(v.array().at(0).number(),v.array().at(1).number());return out;}

// First-order Poisson kernel on the circle, P(t)=(1-q^2)/(2 pi (1-2q cos t+q^2)), q=exp(-bandwidth):
// the order-1 Rigged DMD smoothing kernel (Colbrook, Drysdale & Horning 2024, Sec. 3).
double poisson(double t,double bandwidth){const double q=std::exp(-bandwidth);
    return (1-q*q)/(2*pi*(1-2*q*std::cos(t)+q*q));}

std::vector<double> referenceAngles(){std::vector<double> a;for(int i=0;i<257;++i)a.push_back(-pi+2*pi*i/256);return a;}

struct Pairs { std::vector<double> x,y,w; std::size_t dimension{}; };

J fitRequest(const Pairs& p,const std::vector<double>& coefficients,const std::vector<double>& angles,
             double bandwidth,const std::vector<std::size_t>& orders){
    J::Array measures;
    for(auto order:orders)measures.emplace_back(J::Object{{"name","order"+std::to_string(order)},
        {"coefficients",array(coefficients)},{"angles",array(angles)},{"bandwidth",bandwidth},{"order",double(order)}});
    J::Object data{{"dimension",double(p.dimension)},{"x",array(p.x)},{"y",array(p.y)}};
    if(!p.w.empty())data["weights"]=array(p.w);
    return J::Object{{"schema","dynsyslab.koopman-spectral-query/v1"},{"action","fit"},{"data",data},
        {"settings",J::Object{{"method","mpedmd"},{"dictionary",J::Object{{"kind","linear"},{"constant",false}}},
            {"sample_interval",.1},{"assume_measure_preserving",true},{"rigged_measures",measures}}}};
}
// A fit is a krom: each is written to its own archive, which later actions name.
J fit(const Pairs& p,const std::vector<double>& coefficients,const std::vector<double>& angles,
      double bandwidth,const std::vector<std::size_t>& orders){
    static int count=0;
    const auto path=std::filesystem::temp_directory_path()/("dynsyslab-colbrook-"+std::to_string(++count)+"-"+
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    auto request=fitRequest(p,coefficients,angles,bandwidth,orders);
    request.object()["output_archive"]=path.generic_u8string();
    return dynsyslab::workflow::evaluateKoopmanSpectralQuery(request);
}
// Row i of the fitted operator in the MATLAB/Colbrook convention PX*K = PY: one rollout step from e_i.
std::vector<double> operatorRows(const J& result,std::size_t dimension){
    std::vector<double> rows;
    for(std::size_t i=0;i<dimension;++i){std::vector<double> e(dimension,0.);e[i]=1;
        const auto out=numbers(dynsyslab::workflow::evaluateKoopmanSpectralQuery(J::Object{
            {"schema","dynsyslab.koopman-spectral-query/v1"},{"action","rollout"},{"krom",result.at("archive")},
            {"initial",array(e)},{"steps",1.}}).at("values"));
        require(out.size()==2*dimension,"rollout shape");
        rows.insert(rows.end(),out.begin()+dimension,out.end());}
    return rows;
}
std::vector<double> density(const J& result,std::size_t index){
    return numbers(result.at("spectral").at("measures").array().at(index).at("density"));}
std::vector<cd> packet(const J& result,std::size_t index){
    return complexes(result.at("spectral").at("measures").array().at(index).at("packets"));}  // a measure's generalized eigenfunctions

Pairs rotation(){Pairs p;p.dimension=2;
    for(int i=0;i<128;++i){const double t=2*pi*i/128;
        p.x.insert(p.x.end(),{std::cos(t),std::sin(t)});p.y.insert(p.y.end(),{std::cos(t+.31),std::sin(t+.31)});}
    return p;}

// check_native_reference.m: the pinned upstream mpEDMDqr / riggedDMD on the synthetic rotation,
// frozen at every 16th of 257 angles, plus the exact two-atom answer for this rotation.
void pinnedUpstreamRotation(const std::filesystem::path& file){
    const auto reference=json::load(file.wstring());
    require(reference.at("schema").string()=="dynsyslab.colbrook-rotation-reference/v1","reference schema");
    const double bandwidth=reference.at("measure").at("bandwidth").number();
    const auto angles=referenceAngles();
    const auto result=fit(rotation(),{1,0},angles,bandwidth,{1,4});
    const auto rows=operatorRows(result,2);
    const auto frozen=numbers(reference.at("operator_row_major"));
    const double c=std::cos(.31),s=std::sin(.31);const double exact[]{c,s,-s,c};
    for(std::size_t i=0;i<4;++i){close(rows[i],frozen[i],2e-12,"mpEDMD operator versus pinned upstream");
        close(rows[i],exact[i],1e-12,"mpEDMD operator versus exact rotation");}
    const auto& indices=reference.at("measure").at("angle_indices").array();
    for(std::size_t m=0;m<2;++m){
        const auto order=m?"4":"1";const auto& orderData=reference.at("orders").at(order);
        const auto xi=density(result,m);const auto g=packet(result,m);
        require(xi.size()==angles.size()&&g.size()==2*angles.size(),"rigged result shape");
        const auto frozenDensity=numbers(orderData.at("density"));
        const auto& frozenPackets=orderData.at("packet_coefficients").array();
        for(std::size_t k=0;k<indices.size();++k){
            const auto j=std::size_t(indices[k].number());
            close(xi[j],frozenDensity[k],2e-10,std::string("order-")+order+" density versus pinned upstream");
            const auto upstream=complexes(frozenPackets[k]);
            for(std::size_t f=0;f<2;++f)require(std::abs(g[2*j+f]-upstream[f])<2e-10,
                std::string("order-")+order+" wave packet versus pinned upstream");
        }
        if(m==0)for(std::size_t j=0;j<angles.size();++j){
            // Two atoms of mass 1/4 at +-0.31 (||cos||^2=1/2 under the uniform measure), Poisson-smoothed.
            const double minus=poisson(angles[j]-.31,bandwidth),plus=poisson(angles[j]+.31,bandwidth);
            close(xi[j],.25*(minus+plus),1e-12,"analytic order-1 density");
            require(std::abs(g[2*j]-cd(.5*(minus+plus),0))<1e-12&&std::abs(g[2*j+1]-cd(0,.5*(minus-plus)))<1e-12,
                "analytic order-1 wave packet");
        }
    }
}

// Real weighted dictionary on an irregular circle with two harmonics, advanced by a rotation.
Pairs weightedHarmonics(){Pairs p;p.dimension=4;
    for(int i=0;i<96;++i){const double t=2*pi*(i+.37*std::sin(1.3*i))/96;
        p.x.insert(p.x.end(),{std::cos(t),std::sin(t),std::cos(2*t),std::sin(2*t)});
        const double u=t+.31;p.y.insert(p.y.end(),{std::cos(u),std::sin(u),std::cos(2*u),std::sin(2*u)});
        p.w.push_back(.2+std::fmod(.6180339887*i,1.));}
    return p;}
std::vector<double> gridAngles(){std::vector<double> a;for(int i=0;i<30;++i)a.push_back(-pi+2*pi*i/29);return a;}
double relative(const std::vector<cd>& a,const std::vector<cd>& b){double e=0,n=0;
    for(std::size_t i=0;i<a.size();++i){e+=std::norm(a[i]-b[i]);n+=std::norm(b[i]);}return std::sqrt(e/n);}
double relative(const std::vector<double>& a,const std::vector<double>& b){double e=0,n=0;
    for(std::size_t i=0;i<a.size();++i){e+=(a[i]-b[i])*(a[i]-b[i]);n+=b[i]*b[i];}return std::sqrt(e/n);}

// rigged_history_reference.m / test_rigged_history_reference.m: compressing weighted joint rows
// [sqrt(w) PX, sqrt(w) PY] by block QR keeps every weighted product, so the fit must not change.
void blockQrCompression(){
    const auto raw=weightedHarmonics();const std::vector<double> target{0,0,1,0};
    const auto angles=gridAngles();
    const auto direct=fit(raw,target,angles,.08,{4});
    dynsyslab::numerics::WeightedQrAccumulator qr(8);
    for(std::size_t begin=0;begin<96;begin+=48){std::vector<double> rows,weights;
        for(std::size_t i=begin;i<begin+48;++i){
            rows.insert(rows.end(),raw.x.begin()+4*i,raw.x.begin()+4*i+4);
            rows.insert(rows.end(),raw.y.begin()+4*i,raw.y.begin()+4*i+4);weights.push_back(raw.w[i]);}
        qr.append(rows,48,weights);}
    require(qr.samples()==96&&qr.rows()==8,"block QR bookkeeping");
    // Uniform weights over the R rows are normalized to 1/rows; rescale so the Gram matrices agree.
    const double scale=std::sqrt(double(qr.rows())/qr.totalWeight());
    Pairs compressed;compressed.dimension=4;
    for(std::size_t r=0;r<qr.rows();++r)for(std::size_t c=0;c<8;++c)
        (c<4?compressed.x:compressed.y).push_back(scale*qr.factor()[r*8+c]);
    const auto reduced=fit(compressed,target,angles,.08,{4});
    require(relative(operatorRows(reduced,4),operatorRows(direct,4))<1e-10,"compressed operator");
    require(relative(density(reduced,0),density(direct,0))<1e-10,"compressed density");
    require(relative(packet(reduced,0),packet(direct,0))<1e-10,"compressed wave packet");
}

// compare_rigged_reference.m: the fit is posed in L2(w), so whitening the dictionary by B (with
// B^T PX^T W PX B = I) changes coefficients only: g = B g_white, same operator up to similarity,
// same spectral measure for the same observable (Colbrook 2023, Sec. 3).
void whiteningInvariance(){
    const auto raw=weightedHarmonics();const auto angles=gridAngles();const std::vector<double> target{0,0,1,0};
    double total=0;for(double w:raw.w)total+=w;
    dynsyslab::numerics::WeightedQrAccumulator qr(4);std::vector<double> normalized;
    for(double w:raw.w)normalized.push_back(w/total);
    qr.append(raw.x,96,normalized);
    const auto& r=qr.factor();
    std::vector<double> b(16,0.); // B = R^{-1}, upper triangular, by back substitution
    for(std::size_t col=0;col<4;++col)for(std::size_t i=4;i-->0;){
        double v=(i==col)?1.:0.;for(std::size_t k=i+1;k<4;++k)v-=r[i*4+k]*b[k*4+col];b[i*4+col]=v/r[i*4+i];}
    const auto times=[&](const std::vector<double>& rows){std::vector<double> out(rows.size(),0.);
        for(std::size_t s=0;s<rows.size()/4;++s)for(std::size_t j=0;j<4;++j)for(std::size_t k=0;k<4;++k)
            out[4*s+j]+=rows[4*s+k]*b[k*4+j];return out;};
    Pairs white=raw;white.x=times(raw.x);white.y=times(raw.y);
    double defect=0;
    for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j){double g=0;
        for(std::size_t s=0;s<96;++s)g+=normalized[s]*white.x[4*s+i]*white.x[4*s+j];defect+=std::pow(g-(i==j),2);}
    require(std::sqrt(defect)/2<1e-12,"whitening Gram defect");
    std::vector<double> whiteTarget(4,0.); // B^{-1} g = R g
    for(std::size_t i=0;i<4;++i)for(std::size_t k=0;k<4;++k)whiteTarget[i]+=r[i*4+k]*target[k];
    const auto direct=fit(raw,target,angles,.08,{4}),whitened=fit(white,whiteTarget,angles,.08,{4});
    require(relative(density(whitened,0),density(direct,0))<1e-10,"whitened density");
    const auto gw=packet(whitened,0),gd=packet(direct,0);std::vector<cd> mapped(gd.size());
    for(std::size_t a=0;a<angles.size();++a)for(std::size_t i=0;i<4;++i)for(std::size_t k=0;k<4;++k)
        mapped[4*a+i]+=b[i*4+k]*gw[4*a+k];
    require(relative(mapped,gd)<1e-10,"whitened wave packet maps back through B");
}
}

int main(int argc,char** argv){
    try{
        const std::filesystem::path data=argc>1?std::filesystem::path(argv[1]):
            std::filesystem::path(DYNSYSLAB_REFERENCE_DATA)/"colbrook_rotation_reference_v1.json";
        pinnedUpstreamRotation(data);blockQrCompression();whiteningInvariance();
        std::cout<<"dynsyslab_colbrook_reference_tests: PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<"dynsyslab_colbrook_reference_tests: FAIL: "<<e.what()<<'\n';return 1;}
}
