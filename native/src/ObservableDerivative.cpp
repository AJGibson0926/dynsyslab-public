#include <dynsyslab/core/ObservableDerivative.hpp>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
ResolvedObservableDerivative resolveObservableDerivative(const ODESystem& system,
    const ObservableCatalog& catalog,const std::string& name,std::size_t order,
    std::optional<double> constantInput) {
    if (order>8) throw std::invalid_argument("observable derivative order exceeds 8");
    const auto& base=catalog.at(name);
    if (constantInput && !std::isfinite(*constantInput)) throw std::invalid_argument("nonfinite model input");
    auto map=base;
    bool named=true;
    std::string method="none";
    for (std::size_t i=0;i<order;++i) {
        const auto* link=named?catalog.timeDerivative(map.descriptor().name):nullptr;
        if (link && (!link->zeroInputOnly || constantInput==0.)) {
            map=catalog.at(link->target);
            method=link->zeroInputOnly?"catalog-zero-input":"catalog";
            continue;
        }
        const auto drift=system.driftExpressions();
        if (constantInput && !drift.empty() && !map.expressions().empty()) {
            if (*constantInput==0.) {
                map=map.flowDerivative(drift,order-i);
                method="analytic-zero-input";
            } else {
                const auto fields=system.inputExpressions();
                if (!fields.empty()) {
                    map=map.inputFlowDerivative(drift,fields,{*constantInput},{},order-i);
                    method="analytic-constant-input";
                } else if (i+1==order && map.hasAnalyticJacobian()) {
                    return {map,base.descriptor().name+"_dt"+std::to_string(order),"analytic-jacobian-times-trajectory-rhs",true};
                } else throw std::invalid_argument("missing analytic input field for "+name);
            }
            named=false;
            break;
        }
        if (i+1==order && map.hasAnalyticJacobian())
            return {map,base.descriptor().name+"_dt"+std::to_string(order),"analytic-jacobian-times-trajectory-rhs",true};
        throw std::invalid_argument("time derivative unavailable for "+name+" on this trajectory");
    }
    const auto resolvedName=named?map.descriptor().name:base.descriptor().name+"_dt"+std::to_string(order);
    return {std::move(map),resolvedName,method,false};
}
} // namespace dynsyslab::core
