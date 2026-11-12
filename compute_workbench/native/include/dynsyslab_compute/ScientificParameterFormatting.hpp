#pragma once

#include <string>

namespace dsl::compute {

std::wstring scientificParameterLatex(const std::string& key);
std::wstring scientificParameterNumber(double value, bool longPrecision);
std::wstring scientificParameterUnit(const std::string& unit);
bool parseScientificParameterInput(const std::wstring& input, double& value);

}  // namespace dsl::compute
