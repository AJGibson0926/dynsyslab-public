#include "dynsyslab_compute/ScientificParameterFormatting.hpp"

#include <algorithm>
#include <cmath>
#include <cwchar>
#include <cwctype>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace dsl::compute {
namespace {

bool parsePlain(const std::wstring& input, double& value) {
    const wchar_t* begin = input.c_str();
    wchar_t* end = nullptr;
    value = std::wcstod(begin, &end);
    while (end && std::iswspace(*end)) ++end;
    return end && end != begin && *end == L'\0' && std::isfinite(value);
}

void replaceToken(std::wstring& input, const std::wstring& from,
                  const std::wstring& to) {
    std::size_t position = 0;
    while ((position = input.find(from, position)) != std::wstring::npos) {
        input.replace(position, from.size(), to);
        position += to.size();
    }
}

std::wstring normalizedScientificInput(const std::wstring& input) {
    std::wstring result;
    for (wchar_t character : input) {
        if (std::iswspace(character)) continue;
        switch (character) {
            case L'\u00d7': case L'\u00b7': case L'*': result.push_back(L'x'); break;
            case L'\u2212': case L'\u207b': result.push_back(L'-'); break;
            case L'\u207a': result.push_back(L'+'); break;
            case L'\u2070': result.push_back(L'0'); break;
            case L'\u00b9': result.push_back(L'1'); break;
            case L'\u00b2': result.push_back(L'2'); break;
            case L'\u00b3': result.push_back(L'3'); break;
            case L'\u2074': result.push_back(L'4'); break;
            case L'\u2075': result.push_back(L'5'); break;
            case L'\u2076': result.push_back(L'6'); break;
            case L'\u2077': result.push_back(L'7'); break;
            case L'\u2078': result.push_back(L'8'); break;
            case L'\u2079': result.push_back(L'9'); break;
            default: result.push_back(static_cast<wchar_t>(std::towlower(character))); break;
        }
    }
    replaceToken(result, L"\\times", L"x");
    replaceToken(result, L"\\cdot", L"x");
    return result;
}

std::wstring unwrapExponent(std::wstring exponent) {
    if (!exponent.empty() && exponent.front() == L'^') exponent.erase(0, 1);
    if (exponent.size() >= 2
        && ((exponent.front() == L'{' && exponent.back() == L'}')
            || (exponent.front() == L'(' && exponent.back() == L')')))
        exponent = exponent.substr(1, exponent.size() - 2);
    return exponent;
}

std::wstring superscript(int exponent) {
    static constexpr wchar_t digits[]{L'\u2070', L'\u00b9', L'\u00b2', L'\u00b3', L'\u2074',
                                      L'\u2075', L'\u2076', L'\u2077', L'\u2078', L'\u2079'};
    if (exponent == 0) return std::wstring(1, digits[0]);
    std::wstring result;
    if (exponent < 0) { result.push_back(L'\u207b'); exponent = -exponent; }
    while (exponent) {
        const std::size_t position = !result.empty() && result.front() == L'\u207b' ? 1 : 0;
        result.insert(result.begin() + static_cast<std::ptrdiff_t>(position),
                      digits[exponent % 10]);
        exponent /= 10;
    }
    return result;
}

}  // namespace

std::wstring scientificParameterLatex(const std::string& key) {
    static const std::map<std::string, std::wstring> names{
        {"equilibrium_radius_m", L"R_{0}"}, {"ambient_pressure_pa", L"p_{0}"},
        {"vapor_pressure_pa", L"p_{V}"}, {"rho_kg_m3", L"\\rho_{L}"},
        {"sigma_n_m", L"\\sigma"}, {"mu_pa_s", L"\\mu_{L}"},
        {"kappa", L"\\kappa"}, {"pressure_input_sign", L"s_{u}"},
        {"liquid_sound_speed_m_s", L"c_{L}"}, {"gas_sound_speed_m_s", L"c_{G}"},
        {"shell_viscosity_kg_s", L"\\mu_{s}"}, {"buckling_radius_m", L"R_{b}"},
        {"critical_surface_tension_n_m", L"\\sigma_{c}"},
        {"shell_elasticity_n_m", L"\\chi_{s}"},
        {"acoustic_pressure_amplitude", L"p_{A}"}, {"drive_frequency_hz", L"f_{d}"},
        {"translation_enabled", L"T_{on}"},
        {"gas_pressure_ratio", L"\\Pi_{G}"}, {"vapor_pressure_ratio", L"\\Pi_{V}"},
        {"ambient_pressure_ratio", L"\\Pi_{0}"}, {"capillarity", L"\\Pi_{\\sigma}"},
        {"viscosity", L"\\Pi_{\\mu}"}, {"natural_frequency_hz", L"f_{0}"},
        {"characteristic_time_s", L"\\tau_{0}"}, {"pressure_scale_pa", L"p_{*}"},
        {"liquid_mach", L"Ma_{L}"}, {"gas_mach", L"Ma_{G}"},
        {"euler", L"Eu"}, {"shell_viscosity", L"\\Pi_{s}"},
        {"buckling_radius_ratio", L"R_{b}/R_{0}"},
        {"rupture_radius_ratio", L"R_{r}/R_{0}"},
        {"breaking_radius_ratio", L"R_{c}/R_{0}"},
        {"pressure_amplitude", L"P_{A}"},
        {"pressure_derivative_amplitude", L"P'_{A}"},
        {"angular_frequency", L"\\omega"}, {"wave_number", L"k"},
        {"delta", L"\\delta"}, {"alpha", L"\\alpha"}, {"beta", L"\\beta"},
        {"gamma", L"\\gamma"}, {"mu", L"\\mu"}, {"sigma", L"\\sigma"},
        {"rho", L"\\rho"}, {"lam", L"\\lambda"}, {"input_gain", L"b_{u}"},
        {"mass_kg", L"m"}, {"mass1_kg", L"m_{1}"}, {"mass2_kg", L"m_{2}"},
        {"length_m", L"l"}, {"length1_m", L"l_{1}"}, {"length2_m", L"l_{2}"},
        {"gravity_m_s2", L"g"}, {"damping_n_m_s", L"c"},
        {"damping1_n_m_s", L"c_{1}"}, {"damping2_n_m_s", L"c_{2}"},
        {"pivot_separation_m", L"d"}, {"spring_constant_n_m", L"k_{s}"},
        {"spring_rest_length_m", L"l_{s}"}, {"frequency_spread", L"\\Delta_{\\omega}"},
        {"oscillators", L"N"}, {"sites", L"N"}, {"points", L"N_{x}"},
        {"diffusivity", L"D"}, {"decay_rate", L"\\gamma"},
        {"amplitude", L"A"}, {"offset", L"u_{0}"}, {"boundary_code", L"b_{c}"},
        {"delay", L"\\tau"}, {"exponent", L"n"}, {"history_value", L"x_{h}"},
        {"size", L"N"}, {"states", L"q"}, {"coupling", L"J"},
        {"field", L"h"}, {"temperature", L"T"}, {"current", L"I"},
        {"bias_current", L"I_{0}"}, {"capacitance", L"C_{m}"},
        {"g_na", L"g_{Na}"}, {"g_k", L"g_{K}"}, {"g_l", L"g_{L}"},
        {"e_na", L"E_{Na}"}, {"e_k", L"E_{K}"}, {"e_l", L"E_{L}"},
        {"x_r", L"x_{r}"}};
    const auto found = names.find(key);
    if (found != names.end()) return found->second;
    const auto separator = key.find('_');
    if (separator == std::string::npos) return std::wstring(key.begin(), key.end());
    const std::string base = key.substr(0, separator);
    std::wstring suffix;
    for (std::size_t index = separator + 1; index < key.size(); ++index)
        if (index == separator + 1 || key[index - 1] == '_') suffix.push_back(key[index]);
    return std::wstring(base.begin(), base.end()) + L"_{" + suffix + L"}";
}

std::wstring scientificParameterNumber(double value, bool longPrecision) {
    if (value == 0.0) return L"0";
    std::wostringstream text;
    text << std::setprecision(longPrecision ? std::numeric_limits<double>::max_digits10 : 4)
         << std::defaultfloat << value;
    const std::wstring encoded = text.str();
    const auto marker = encoded.find_first_of(L"eE");
    if (marker == std::wstring::npos) return encoded == L"-0" ? L"0" : encoded;
    const int exponent = std::stoi(encoded.substr(marker + 1));
    return encoded.substr(0, marker) + L" \u00d7 10" + superscript(exponent);
}

std::wstring scientificParameterUnit(const std::string& unit) {
    if (unit.empty() || unit == "dimensionless") return {};
    std::wstring result(unit.begin(), unit.end());
    const std::pair<const wchar_t*, const wchar_t*> replacements[]{
        {L"m^3", L"m\u00b3"}, {L"m^2", L"m\u00b2"}, {L"s^2", L"s\u00b2"},
        {L"1/s", L"s\u207b\u00b9"}, {L"1/time", L"time\u207b\u00b9"}};
    for (const auto& item : replacements) {
        const auto at = result.find(item.first);
        if (at != std::wstring::npos)
            result.replace(at, std::wcslen(item.first), item.second);
    }
    return result;
}

bool parseScientificParameterInput(const std::wstring& input, double& value) {
    if (parsePlain(input, value)) return true;
    const std::wstring normalized = normalizedScientificInput(input);
    if (parsePlain(normalized, value)) return true;
    const auto product = normalized.find(L"x10");
    const bool powerOnly = product == std::wstring::npos
        && normalized.rfind(L"10", 0) == 0 && normalized.size() > 2;
    if (product == std::wstring::npos && !powerOnly) return false;
    const std::wstring mantissaText = powerOnly ? L"1" : normalized.substr(0, product);
    const std::size_t exponentStart = powerOnly ? 2 : product + 3;
    const std::wstring exponentText = unwrapExponent(normalized.substr(exponentStart));
    double mantissa = 0.0, exponent = 0.0;
    if (!parsePlain(mantissaText, mantissa) || !parsePlain(exponentText, exponent))
        return false;
    value = mantissa * std::pow(10.0, exponent);
    return std::isfinite(value);
}

}  // namespace dsl::compute
