#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <string>
#include <vector>

namespace dsl::compute {

struct PresetDefinition {
    std::string key;
    std::string description;
};

struct SystemDefinition {
    std::string name;
    std::string label;
    int stateDimension = 0;
    std::string family{"ode"};
    std::string preferredCoordinates;
    std::vector<std::string> stateLabels;
    std::vector<PresetDefinition> presets;
};

struct KromDefinition {
    std::string label;
    std::string source;
    std::string system;
};

struct Registry {
    std::vector<SystemDefinition> systems;
    std::vector<KromDefinition> kroms;
    std::vector<KromDefinition> spectra;
    std::vector<KromDefinition> trajectories;
    std::vector<std::string> targets;
    std::vector<std::string> methods;
    static Registry load(const std::wstring& path);
    static Registry builtin(const dsl::Catalog& catalog);
};

}  // namespace dsl::compute
