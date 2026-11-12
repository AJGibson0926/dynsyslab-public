#include "BubbleSupport.hpp"

namespace dynsyslab::core::detail {

std::unique_ptr<ODESystem> makeBubbleSystem(const SystemConfiguration& configuration) {
    if (configuration.identifier == "rpe") return makeRPE(configuration);
    if (configuration.identifier == "marmottant") return makeMarmottant(configuration);
    if (configuration.identifier == "reddy_szeri") return makeReddySzeri(configuration);
    return {};
}

}  // namespace dynsyslab::core::detail
