#pragma once

#include <map>
#include <string>
#include <variant>
#include <vector>

namespace dynsyslab::core::json {

class Value {
public:
    using Array = std::vector<Value>;
    using Object = std::map<std::string, Value>;
    using Storage = std::variant<std::nullptr_t, bool, double, std::string, Array, Object>;

    Value() : data_(nullptr) {}
    explicit Value(Storage data) : data_(std::move(data)) {}
    Value(bool value) : data_(value) {}
    Value(double value) : data_(value) {}
    Value(std::string value) : data_(std::move(value)) {}
    Value(const char* value) : data_(std::string(value)) {}
    Value(Array value) : data_(std::move(value)) {}
    Value(Object value) : data_(std::move(value)) {}

    bool isNull() const;
    bool boolean(bool fallback = false) const;
    double number(double fallback = 0.0) const;
    std::string string(std::string fallback = {}) const;
    const Array& array() const;
    const Object& object() const;
    Array& array();
    Object& object();
    const Value& at(const std::string& key) const;
    Value& at(const std::string& key);
    const Storage& storage() const { return data_; }

private:
    Storage data_;
};

Value parse(const std::string& text);
Value load(const std::wstring& path);
std::string dump(const Value& value);
void setPath(Value& root, const std::string& dottedPath, Value value);

}  // namespace dynsyslab::core::json
