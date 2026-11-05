#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>

namespace dynsyslab::core::detail {
inline json::Value koopmanSystemJson(const SystemConfiguration& source) {
    json::Value::Object parameters;
    for (const auto& item : source.parameters) parameters[item.first] = item.second;
    json::Value::Object units;
    for (const auto& item : source.parameterUnits) units[item.first] = item.second;
    json::Value::Object groups;
    for (const auto& item : source.dimensionlessGroups) groups[item.first] = item.second;
    json::Value::Object groupUnits;
    for (const auto& item : source.dimensionlessGroupUnits) groupUnits[item.first] = item.second;
    json::Value::Array initial;
    for (const double value : source.initialState) initial.emplace_back(value);
    json::Value::Object result;
    result["identifier"] = source.identifier;
    result["preset"] = source.preset;
    result["parameters"] = json::Value(std::move(parameters));
    result["parameter_units"] = json::Value(std::move(units));
    result["dimensionless_groups"] = json::Value(std::move(groups));
    result["dimensionless_group_units"] = json::Value(std::move(groupUnits));
    result["equations_version"] = source.equationsVersion;
    result["initial_state"] = json::Value(std::move(initial));
    return json::Value(std::move(result));
}
}
