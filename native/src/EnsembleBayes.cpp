#include <dynsyslab/core/EnsembleBayes.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
void finite(const std::vector<double>& x){for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("Finite ensemble values required");}
double normalize(std::vector<double>& logWeights){
    const double top=*std::max_element(logWeights.begin(),logWeights.end());
    if(!std::isfinite(top))throw std::runtime_error("Every ensemble likelihood vanished");
    double sum=0;for(double v:logWeights)sum+=std::exp(v-top);
    const double logZ=top+std::log(sum);for(double& v:logWeights)v-=logZ;return logZ;
}
}
EnsembleBayesFilter::EnsembleBayesFilter(std::vector<State> particles,std::vector<double> prior):particles_(std::move(particles)){
    if(particles_.empty()||particles_[0].empty())throw std::invalid_argument("Nonempty ensemble required");
    if(prior.empty())prior.assign(particles_.size(),1.);
    if(prior.size()!=particles_.size())throw std::invalid_argument("Ensemble prior size mismatch");
    finite(prior);
    for(std::size_t i=0;i<particles_.size();++i){finite(particles_[i]);
        if(particles_[i].size()!=particles_[0].size()||prior[i]<0)throw std::invalid_argument("Invalid ensemble state or prior");
        logWeights_.push_back(prior[i]>0?std::log(prior[i]):-std::numeric_limits<double>::infinity());
    }normalize(logWeights_);
}
void EnsembleBayesFilter::predict(const Transition& f){
    auto next=particles_;for(auto& x:next){x=f(x);finite(x);if(x.size()!=particles_[0].size())throw std::invalid_argument("Ensemble transition changed dimension");}
    particles_=std::move(next);
}
double EnsembleBayesFilter::correct(const State& y,const State& variance,const Readout& h){
    finite(y);finite(variance);if(y.empty()||y.size()!=variance.size())throw std::invalid_argument("Measurement variance shape mismatch");
    double constant=0;for(double v:variance){if(v<=0)throw std::invalid_argument("Positive measurement variances required");constant+=std::log(2*std::acos(-1.)*v);}
    auto next=logWeights_;
    for(std::size_t i=0;i<particles_.size();++i){const auto predicted=h(particles_[i]);finite(predicted);
        if(predicted.size()!=y.size())throw std::invalid_argument("Observation dimension mismatch");
        double squared=0;for(std::size_t j=0;j<y.size();++j){const double delta=(y[j]-predicted[j])/std::sqrt(variance[j]);squared+=delta*delta;}
        next[i]-=.5*(constant+squared);
    }
    const double evidence=normalize(next);logWeights_=std::move(next);return evidence;
}
std::vector<double> EnsembleBayesFilter::weights() const {std::vector<double> w;for(double x:logWeights_)w.push_back(std::exp(x));return w;}
double EnsembleBayesFilter::effectiveSampleSize() const {double sum=0;for(double w:weights())sum+=w*w;return 1/sum;}
EnsembleBayesFilter::State EnsembleBayesFilter::mean(const Readout& h) const {
    State out;const auto w=weights();for(std::size_t i=0;i<particles_.size();++i){const auto y=h(particles_[i]);finite(y);
        if(i==0)out.assign(y.size(),0);if(y.size()!=out.size())throw std::invalid_argument("Ensemble readout changed dimension");
        for(std::size_t j=0;j<y.size();++j)out[j]+=w[i]*y[j];
    }finite(out);return out;
}
}
