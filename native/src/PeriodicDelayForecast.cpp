#include <dynsyslab/core/PeriodicSpectralComparison.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
// Fourier translation and Parseval: Katznelson, Introduction to Harmonic
// Analysis, Ch. I. Finite Koopman propagation uses the existing spectral API.
namespace dynsyslab::core {
std::vector<PeriodicDelayForecastAudit> auditPeriodicDelayForecast(const SpectralKROM&model,
    const std::vector<math::PeriodicSnapshotOrbit>&bank,const std::vector<std::size_t>&steps,
    const std::vector<std::vector<std::size_t>>&evaluationSteps){
    model.validate();math::auditPeriodicQuadrature(bank,{});
    if(steps.size()!=bank.size())throw std::invalid_argument("One forecast step count per orbit required");
    const auto contract=json::parse(model.coordinateContract);
    if(contract.at("schema").string()!="dynsyslab.periodic-delay-coordinates/v1"||
       contract.at("sample_interval").number()!=model.sampleInterval||model.dictionary.kind!="linear"||
       model.dictionary.constant||!model.dictionary.center.empty()||!model.dictionary.scale.empty())
        throw std::invalid_argument("Forecast audit needs an unnormalized linear periodic delay model with a verified clock");
    std::vector<std::size_t> channels;std::vector<double> offsets;
    for(const auto&v:contract.at("channels").array()){
        const double c=v.number(-1);if(!std::isfinite(c)||c<0||c!=std::floor(c)||c>=bank.front().channels.size())
            throw std::invalid_argument("Forecast channel outside bank");channels.push_back(std::size_t(c));}
    for(const auto&v:contract.at("time_offsets").array()){
        if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite offsets required");
        offsets.push_back(v.number());}
    const auto nc=channels.size(),n=model.dictionary.size(),r=model.rank;
    if(!nc||offsets.empty()||n!=nc*offsets.size()||model.dictionary.dimension!=n)
        throw std::invalid_argument("Forecast delay dimensions mismatch");
    if(!evaluationSteps.empty()&&evaluationSteps.size()!=bank.size())throw std::invalid_argument("One evaluation-step list per orbit required");
    std::vector<std::vector<std::size_t>> selected(bank.size());std::set<std::size_t> combined;
    std::size_t total=0;
    for(std::size_t j=0;j<bank.size();++j){
        if(steps[j]>1000000)throw std::invalid_argument("Forecast step exceeds one million");
        auto& list=selected[j];
        if(evaluationSteps.empty()){list.resize(steps[j]+1);for(std::size_t k=0;k<list.size();++k)list[k]=k;}
        else list=evaluationSteps[j];
        if(list.empty()||list.front()!=0||list.back()!=steps[j]||
           !std::is_sorted(list.begin(),list.end())||std::adjacent_find(list.begin(),list.end())!=list.end())
            throw std::invalid_argument("Evaluation steps must increase strictly from zero to the requested last step");
        if(list.size()>(10000000-total)/nc)throw std::invalid_argument("Forecast audit exceeds ten million time/channel entries");
        total+=list.size()*nc;combined.insert(list.begin(),list.end());
    }
    const std::vector<std::size_t> globalSteps(combined.begin(),combined.end());
    if(globalSteps.size()>50000000/nc||n>50000000/(globalSteps.size()*nc))
        throw std::invalid_argument("Forecast coefficient workspace exceeds fifty million entries");
    const auto latest=std::size_t(std::max_element(offsets.begin(),offsets.end())-offsets.begin());
    const double cutoff=offsets[latest];for(auto&v:offsets)v-=cutoff;
    const auto&readout=model.reconstructions.at("coordinates");
    if(readout.outputs!=n)throw std::invalid_argument("Expected coordinate readout");
    std::vector<double> target(n*nc);
    for(std::size_t i=0;i<n;++i)for(std::size_t c=0;c<nc;++c)
        target[i*nc+c]=readout.coefficients[i*n+latest*nc+c];
    auto latent=numerics::solveLeastSquares(model.basis,n,r,target,nc);
    const auto columns=globalSteps.size()*nc;std::vector<double> coefficients(n*columns);
    // Binary matrix powers reuse long gaps without an eigenbasis/logarithm.
    // Golub and Van Loan, Matrix Computations: repeated squaring. Small gaps
    // retain ordinary repeated propagation; floating-point roundoff can differ.
    std::vector<std::vector<double>> powers{model.action};std::size_t previous=0;
    for(std::size_t index=0;index<globalSteps.size();++index){
        auto gap=globalSteps[index]-previous;
        if(gap<=4){for(std::size_t k=0;k<gap;++k)latent=numerics::multiply(model.action,r,r,latent,nc);}
        else{std::size_t bit=0;while(gap){
            if(bit==powers.size()){
                if(powers.size()+1>50000000/r/r)throw std::invalid_argument("Forecast power cache exceeds fifty million entries");
                powers.push_back(numerics::multiply(powers.back(),r,r,powers.back(),r));}
            if(gap&1)latent=numerics::multiply(powers[bit],r,r,latent,nc);
            gap>>=1;++bit;}}
        const auto block=numerics::multiply(model.basis,n,r,latent,nc);
        for(std::size_t i=0;i<n;++i)for(std::size_t c=0;c<nc;++c)coefficients[i*columns+index*nc+c]=block[i*nc+c];
        previous=globalSteps[index];
    }
    for(double x:coefficients)if(!std::isfinite(x))throw std::runtime_error("Forecast coefficients overflow");
    std::vector<PeriodicDelayForecastAudit> results;
    for(std::size_t orbit=0;orbit<bank.size();++orbit){
        const auto&data=bank[orbit];const double omega=data.channels.front().angularFrequency;
        std::set<int> harmonics;
        for(auto c:channels)for(const auto&[h,a]:data.channels[c].coefficients){
            const auto it=data.channels[c].coefficients.find(-h);
            const auto opposite=it==data.channels[c].coefficients.end()?std::complex<double>{}:it->second;
            if(std::abs(a-std::conj(opposite))>1e-10*(1+std::abs(a)))throw std::invalid_argument("Real Fourier channels required");
            harmonics.insert(h);}
        if(harmonics.size()>50000000/n||harmonics.size()>50000000/columns)
            throw std::invalid_argument("Forecast Fourier workspace exceeds fifty million entries");
        if(harmonics.empty()){
            PeriodicDelayForecastAudit zero;zero.steps=steps[orbit];zero.channels=nc;zero.stepIndices=selected[orbit];
            zero.referenceMeanSquare.resize(nc);zero.predictionMeanSquare.resize(zero.stepIndices.size()*nc);
            zero.errorMeanSquare.resize(zero.predictionMeanSquare.size());zero.persistenceErrorMeanSquare.resize(zero.predictionMeanSquare.size());results.push_back(std::move(zero));continue;
        }
        std::vector<double> real(harmonics.size()*n),imag(real.size());std::size_t hrow=0;
        for(int h:harmonics){
            for(std::size_t tap=0;tap<offsets.size();++tap){const auto phase=std::exp(std::complex<double>(0,h*omega*offsets[tap]));
                for(std::size_t c=0;c<nc;++c){const auto it=data.channels[channels[c]].coefficients.find(h);
                    const auto a=it==data.channels[channels[c]].coefficients.end()?std::complex<double>{}:it->second*phase;
                    real[hrow*n+tap*nc+c]=a.real();imag[hrow*n+tap*nc+c]=a.imag();}}
            ++hrow;
        }
        const auto pr=numerics::multiply(real,harmonics.size(),n,coefficients,columns);
        const auto pi=numerics::multiply(imag,harmonics.size(),n,coefficients,columns);
        PeriodicDelayForecastAudit out;out.steps=steps[orbit];out.channels=nc;out.stepIndices=selected[orbit];
        std::vector<std::size_t> slots;for(auto k:out.stepIndices)slots.push_back(std::lower_bound(globalSteps.begin(),globalSteps.end(),k)-globalSteps.begin());
        out.referenceMeanSquare.resize(nc);out.predictionMeanSquare.resize(out.stepIndices.size()*nc);out.errorMeanSquare.resize(out.predictionMeanSquare.size());
        out.persistenceErrorMeanSquare.resize(out.predictionMeanSquare.size());
        for(std::size_t c=0;c<nc;++c)out.referenceMeanSquare[c]=math::meanSquare(data.channels[channels[c]]);
        hrow=0;for(int h:harmonics){
            for(std::size_t k=0;k<out.stepIndices.size();++k){const auto phase=std::exp(std::complex<double>(0,h*omega*out.stepIndices[k]*model.sampleInterval));
                for(std::size_t c=0;c<nc;++c){const auto j=k*nc+c,index=hrow*columns+slots[k]*nc+c;
                    const auto it=data.channels[channels[c]].coefficients.find(h);
                    const auto initial=it==data.channels[channels[c]].coefficients.end()?std::complex<double>{}:it->second;
                    const auto truth=initial*phase;
                    const std::complex<double> predicted{pr[index],pi[index]};
                    out.predictionMeanSquare[j]+=std::norm(predicted);out.errorMeanSquare[j]+=std::norm(predicted-truth);
                    out.persistenceErrorMeanSquare[j]+=std::norm(initial-truth);}}
            ++hrow;
        }
        for(const auto* array:{&out.predictionMeanSquare,&out.errorMeanSquare,&out.persistenceErrorMeanSquare})for(double value:*array)
            if(!std::isfinite(value))throw std::runtime_error("Forecast phase mean overflow");
        results.push_back(std::move(out));
    }
    return results;
}
}
