#include <dynsyslab/workflow/SystemResolution.hpp>

#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/PDE.hpp>

#include <stdexcept>

namespace dynsyslab::workflow {

ResolvedSystemDefinition resolveSystemDefinition(
    const dynsyslab::core::SystemConfiguration& requested) {
    if (requested.identifier=="mackey_glass") {
        const auto resolved=dynsyslab::core::resolveMackeyGlass(requested);
        return {dynsyslab::core::mackeyGlassDescriptor(),resolved,{}};
    }
    if (requested.identifier=="heat_equation_1d") {
        const auto resolved=dynsyslab::core::resolveHeatEquation(requested);
        const auto points=static_cast<std::size_t>(resolved.parameters.at("points"));
        return {dynsyslab::core::heatEquationDescriptor(points),resolved,{}};
    }
    if (dynsyslab::core::isLatticeSystemIdentifier(requested.identifier)) {
        const auto resolved=dynsyslab::core::resolveLatticeSystem(requested);
        const auto descriptor=resolved.identifier=="ising"
            ? dynsyslab::core::ising2DDescriptor()
            : resolved.identifier=="ising_3d"?dynsyslab::core::ising3DDescriptor()
            : dynsyslab::core::potts2DDescriptor();
        return {descriptor,resolved,{1.0,1.0,false,false,"sweep","field","field"}};
    }
    auto system=dynsyslab::core::createSystem(requested);
    return {system->descriptor(),system->configuration(),system->coordinates()};
}

}  // namespace dynsyslab::workflow
