#include "SettingsMerge.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <stdexcept>

namespace dsl::compute {

std::string mergeSettingsObject(
    std::string object,
    const std::vector<std::pair<std::string, std::string>>& replacements) {
    object = jsonObject(object);
    try {
        auto parsed = dsl::json::parse(object);
        if (!std::holds_alternative<dsl::json::Value::Object>(parsed.storage()))
            throw std::invalid_argument("settings must be a JSON object");
        auto& values = parsed.object();
        for (const auto& item : replacements)
            values[item.first] = dsl::json::parse(item.second);
        return dsl::json::dump(parsed);
    } catch (...) {
        // Preserve malformed user text for the normal plan/validation gate.
    }
    const auto start = object.find('{'), end = object.rfind('}');
    if (start == std::string::npos || end == std::string::npos || end < start)
        return object;
    std::string body = object.substr(start + 1, end - start - 1);
    for (const auto& item : replacements) {
        if (!body.empty()) body += ',';
        body += '"' + jsonEscape(item.first) + "\":" + item.second;
    }
    return '{' + body + '}';
}

}  // namespace dsl::compute
