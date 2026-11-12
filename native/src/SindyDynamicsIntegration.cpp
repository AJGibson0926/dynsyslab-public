#include <dynsyslab/core/SindyDynamics.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
// Adapter only: integration, tolerances and failure handling belong to the
// existing solver. No SINDy-specific time stepper.
class DynamicsSystem final:public ODESystem {
public:
    DynamicsSystem(const SindyDynamicsModel& model,const std::vector<double>& initial)
        :model_(model),field_(model.vectorField()) {
        if(initial.size()!=model.dimension())throw std::invalid_argument("SINDy rollout IC dimension mismatch");
        descriptor_.identifier="sindy-dynamics";descriptor_.label="Identified dynamics";
        descriptor_.stateDimension=model.dimension();descriptor_.preferredCoordinates="model";
        config_.identifier=descriptor_.identifier;config_.initialState=initial;
    }
    const SystemDescriptor& descriptor() const override{return descriptor_;}
    const SystemConfiguration& configuration() const override{return config_;}
    std::vector<double> initialState() const override{return config_.initialState;}
    void rhs(double,const double* q,double* out) const override {
        const auto v=model_.evaluate(q);std::copy(v.begin(),v.end(),out);
    }
    void rhsWithInput(double t,const double* q,double u,double* out) const override {
        if(u!=0)throw std::invalid_argument("Autonomous SINDy dynamics has no input map");
        rhs(t,q,out);
    }
    std::vector<double> jacobian(double,const double* q) const override {
        const auto complex=field_.jacobian(q);std::vector<double> real;
        for(auto v:complex)real.push_back(v.real());return real;
    }
    std::vector<double> inputMatrix(const double*) const override{return std::vector<double>(model_.dimension());}
    bool invalid(double,const double* q) const override {
        for(std::size_t j=0;j<model_.dimension();++j)
            if(!std::isfinite(q[j])||std::abs(q[j])>1e12)return true;
        return false;
    }
    std::vector<AnalyticExpression> driftExpressions() const override{return field_.expressions();}
private:
    SindyDynamicsModel model_;Observable field_;SystemDescriptor descriptor_;SystemConfiguration config_;
};
}
NativeSolution integrateSindyDynamics(const SindyDynamicsModel& model,const std::vector<double>& initial,
    const IntegrationOptions& options) {
    DynamicsSystem system(model,initial);return integrate(system,initial,1,options);
}
}
