#pragma once
// JSON helpers shared by experiment expansion and execution.
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace dynsyslab::workflow::experiment_detail {
inline bool isObject(const json::Value& value) {
    return std::holds_alternative<json::Value::Object>(value.storage());
}
inline json::Value readJson(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::invalid_argument("experiment base is unavailable: " + path.u8string());
    std::stringstream buffer; buffer << stream.rdbuf();
    return json::parse(buffer.str());
}
}  // namespace dynsyslab::workflow::experiment_detail
