#include <dynsyslab/core/EstimationRealization.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
namespace dynsyslab::core {
RealModalRealization makeRealModalRealization(const KoopmanModel& parent,
    const std::vector<std::size_t>& selected,double tolerance,double realTolerance) {
    const auto* source=parent.scatteredEigenfunctions();
    if(!source||selected.empty()||!(tolerance>0)||!std::isfinite(tolerance)
        ||realTolerance<0||!std::isfinite(realTolerance))
        throw std::invalid_argument("real modal realization requires sampled eigenfunctions and representatives");
    const auto n=parent.eigenfunctionCount(),N=source->samples,m=selected.size();
    ScatteredEigenfunctionData data;data.trainingStates=source->trainingStates;
    data.samples=N;data.stateDimension=source->stateDimension;data.neighbors=source->neighbors;
    data.finiteDifferenceRelativeStep=source->finiteDifferenceRelativeStep;
    data.trainingValues.assign(N*m,{});data.provenance="RMS normalized conjugate-reduced H-DMD realization";
    std::vector<std::size_t> partners;std::vector<std::complex<double>> factors,alpha;
    std::set<std::size_t> consumed;
    for(std::size_t j=0;j<m;++j) {
        const auto index=selected[j];if(index>=n)throw std::out_of_range("modal representative exceeds KROM");
        auto rate=parent.continuousEigenvalues()[index];
        const bool real=std::abs(rate.imag())<=realTolerance*std::max(1.0,std::abs(rate));
        if(real)rate={rate.real(),0.0};
        std::size_t partner=index;
        if(!real) {
            double best=std::numeric_limits<double>::infinity();
            for(std::size_t k=0;k<n;++k)if(k!=index) {
                const auto distance=std::abs(parent.continuousEigenvalues()[k]-std::conj(rate));
                if(distance<best){best=distance;partner=k;}
            }
            if(best>tolerance*std::max(1.0,std::abs(rate)))throw std::invalid_argument("missing conjugate eigenvalue");
        }
        if(!consumed.insert(index).second||(!real&&!consumed.insert(partner).second))
            throw std::invalid_argument("duplicate/conjugate representatives supplied");
        double energy=0,maximum=0;std::complex<double> phase=1.,cross{};
        for(std::size_t k=0;k<N;++k) {
            const auto value=source->trainingValues[k*n+index];energy+=std::norm(value);
            if(std::abs(value)>maximum){maximum=std::abs(value);phase=value/maximum;}
            if(!real)cross+=value*source->trainingValues[k*n+partner];
        }
        if(!(energy>1e-50))throw std::invalid_argument("zero-energy modal coordinate");
        const auto divisor=std::sqrt(energy/double(N))*(real?phase:std::complex<double>{1.,0.});
        const auto ratio=real?std::complex<double>{1.,0.}:cross/energy;
        double error=0,partnerEnergy=0;
        for(std::size_t k=0;k<N;++k) {
            const auto value=source->trainingValues[k*n+index];
            if(real) {
                const auto normalized=value/divisor;error+=normalized.imag()*normalized.imag();
                partnerEnergy+=std::norm(normalized);data.trainingValues[k*m+j]=normalized.real();
            } else {
                const auto paired=source->trainingValues[k*n+partner];
                error+=std::norm(paired-ratio*std::conj(value));partnerEnergy+=std::norm(paired);
                data.trainingValues[k*m+j]=value/divisor;
            }
        }
        if(error>tolerance*tolerance*std::max(partnerEnergy,1e-50))
            throw std::invalid_argument("sampled coordinates do not satisfy the real/conjugate relation");
        data.continuousEigenvalues.push_back(rate);partners.push_back(partner);
        factors.push_back(divisor);alpha.push_back(ratio);
    }
    KoopmanModel model(std::move(data),parent.systemConfiguration());
    for(const auto& item:parent.reconstructions()) {
        ReconstructionResult map;map.outputs=item.second.outputs;map.modes=m;map.matrix.resize(map.outputs*m);
        for(std::size_t row=0;row<map.outputs;++row)for(std::size_t j=0;j<m;++j) {
            const auto index=selected[j],partner=partners[j];
            auto coefficient=item.second.matrix[row*n+index];
            if(partner!=index)coefficient+=std::conj(item.second.matrix[row*n+partner]*alpha[j]);
            coefficient*=factors[j];
            map.matrix[row*m+j]=partner==index?std::complex<double>{coefficient.real(),0.}:coefficient;
        }
        model.setReconstruction(item.first,std::move(map));
    }
    for(const auto& ref:parent.datasetReferences())model.addDatasetReference(ref);
    return {std::move(model),selected,std::move(partners),std::move(factors)};
}
} // namespace dynsyslab::core
