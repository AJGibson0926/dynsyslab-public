#include "CatalogPresentationText.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace dsl {

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::string token(const std::string& value) {
    std::string result;
    for (const unsigned char c : value)
        if (std::isalnum(c)) result.push_back(static_cast<char>(std::tolower(c)));
    return result;
}

bool has(const std::string& text, const std::string& part) {
    return lower(text).find(lower(part)) != std::string::npos;
}

std::string property(const CatalogItem& item, const char* name) {
    const auto found = item.properties.find(name);
    return found == item.properties.end() ? std::string{} : found->second;
}

std::string compactSystem(const std::string& system) {
    const auto key = token(system);
    if (key == "rayleighplesset" || key == "rpe")
        return "Bubble · Rayleigh-Plesset";
    if (key == "marmottant") return "Bubble · Marmottant";
    if (key == "reddyszeri") return "Bubble · Reddy-Szeri";
    return system.empty() ? "Other" : system;
}

std::string titleWords(std::string value) {
    for (char& c : value)
        if (c == '_' || c == '-') c = ' ';
    std::string result, word;
    std::istringstream input(value);
    while (input >> word) {
        const auto key = lower(word);
        if (key == "default" || key == "canonical") continue;
        word[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(word[0])));
        if (!result.empty()) result += ' ';
        result += word;
    }
    return result;
}

std::string compactPhysics(const CatalogItem& item) {
    std::string raw = property(item, "physics preset");
    if (raw.empty()) raw = property(item, "system preset");
    if (raw.empty()) raw = property(item, "preset");
    const auto key = token(raw);
    const auto system = token(item.system);
    if ((system == "rayleighplesset" || system == "rpe") &&
        (key.find("air") != std::string::npos && key.find("water") != std::string::npos)) {
        auto temperature = property(item, "temperature");
        if (temperature.empty() && key.find("20c") != std::string::npos) temperature = "20 C";
        if (!temperature.empty()) return "Air in water · " + temperature;
        return key.find("reference") != std::string::npos
            ? "Air in water · reference" : "Air in water";
    }
    if (raw.empty()) return "General parameters";
    auto result = titleWords(raw);
    return result.empty() ? "General parameters" : result;
}

}  // namespace dsl
