#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/System.hpp>

#include <stdexcept>

namespace dsl::compute {

Registry Registry::load(const std::wstring& path) {
    const auto root = dsl::json::load(path);
    if (root.at("schema").string() != "dynsyslab.compute-registry/v1")
        throw std::runtime_error("Unsupported compute registry schema");
    Registry result;
    for (const auto& item : root.at("systems").array()) {
        SystemDefinition system;
        system.name = item.at("name").string(); system.label = item.at("label").string(system.name);
        system.stateDimension = static_cast<int>(item.at("state_dimension").number());
        system.family = item.at("family").string("ode");
        system.preferredCoordinates = item.at("preferred_coordinates").string();
        for (const auto& label : item.at("state_labels").array()) system.stateLabels.push_back(label.string());
        for (const auto& entry : item.at("presets").array())
            system.presets.push_back({entry.at("key").string(), entry.at("description").string()});
        result.systems.push_back(std::move(system));
    }
    for (const auto& item : root.at("kroms").array()) {
        KromDefinition value{item.at("label").string(),item.at("source").string(),
            item.at("system").string()};
        result.kroms.push_back(value);result.spectra.push_back(std::move(value));
    }
    for(const auto& item:root.at("spectra").array())
        result.spectra.push_back({item.at("label").string(),item.at("source").string(),
            item.at("system").string()});
    for(const auto& item:root.at("trajectories").array())
        result.trajectories.push_back({item.at("label").string(),item.at("source").string(),
            item.at("system").string()});
    for (const auto& item : root.at("targets").array()) result.targets.push_back(item.at("name").string());
    for (const auto& item : root.at("postprocessing_methods").array()) result.methods.push_back(item.string());
    return result;
}

Registry Registry::builtin(const dsl::Catalog& catalog) {
    Registry result;
    for (const auto& descriptor : dynsyslab::core::availableSystems()) {
        SystemDefinition system{descriptor.identifier, descriptor.label,
            static_cast<int>(descriptor.stateDimension), descriptor.family,
            descriptor.preferredCoordinates,
            descriptor.stateLabels, {}};
        for (const auto& preset : descriptor.presets)
            system.presets.push_back({preset.key, preset.description});
        result.systems.push_back(std::move(system));
    }

    for (const auto& item : catalog.items)
        if (item.kind == "krom") {
            result.kroms.push_back({item.label, item.source, item.system});
            result.spectra.push_back({item.label, item.source, item.system});
        } else if (item.kind == "koopman_spectrum")
            result.spectra.push_back({item.label,item.source,item.system});
        else if (item.kind == "trajectory")
            result.trajectories.push_back({item.label,item.source,item.system});
    result.targets = {"rpe_multisine_equilibrium_small",
                      "rpe_multisine_equilibrium_moderate_identification",
                      "rpe_multisine_offset_3p5_periodic", "rpe_multisine_offset_4_mixed",
                      "rpe_multisine_offset_4_periodic", "rpe_multisine_offset_2p85_mixed",
                      "rpe_single_sine_offset_2_experimental"};
    result.methods = {"koopman_spectral", "resdmd", "mpedmd", "hermitian_dmd", "tlsdmd", "multdmd", "stochastic_resdmd", "specrkhs",
                      "bagged_pca", "dmd", "dmdc", "edmd", "edmdc", "hankel_pod", "hankel_svd", "hankel_matrix", "hankel_dmd", "havok", "hdmd",
                      "pca", "pod", "robust_pca", "sindy", "svd", "control_regression",
                      "krom_extract", "koopman_validation", "boundary_eigenfunctions",
                      "eigenvalue_optimization", "observability", "klqe"};
    return result;
}

}  // namespace dsl::compute
