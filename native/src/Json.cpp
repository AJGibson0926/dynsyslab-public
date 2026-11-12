#include <dynsyslab/core/Json.hpp>

#include <cctype>
#include <charconv>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>

namespace dynsyslab::core::json {
namespace {

const Value::Array kEmptyArray;
const Value::Object kEmptyObject;
const Value kNull;

void appendUtf8(std::string& target, unsigned value) {
    if (value <= 0x7f) target.push_back(static_cast<char>(value));
    else if (value <= 0x7ff) {
        target.push_back(static_cast<char>(0xc0 | (value >> 6)));
        target.push_back(static_cast<char>(0x80 | (value & 0x3f)));
    } else {
        target.push_back(static_cast<char>(0xe0 | (value >> 12)));
        target.push_back(static_cast<char>(0x80 | ((value >> 6) & 0x3f)));
        target.push_back(static_cast<char>(0x80 | (value & 0x3f)));
    }
}

class Parser {
public:
    explicit Parser(const std::string& source) : source_(source) {}

    Value run() {
        Value result = value();
        whitespace();
        if (position_ != source_.size()) fail("unexpected trailing content");
        return result;
    }

private:
    const std::string& source_;
    std::size_t position_ = 0;

    [[noreturn]] void fail(const char* message) const {
        throw std::runtime_error(std::string("JSON ") + message + " at byte " +
                                 std::to_string(position_));
    }

    void whitespace() {
        while (position_ < source_.size() &&
               std::isspace(static_cast<unsigned char>(source_[position_]))) ++position_;
    }

    bool take(char value) {
        whitespace();
        if (position_ < source_.size() && source_[position_] == value) {
            ++position_;
            return true;
        }
        return false;
    }

    Value value() {
        whitespace();
        if (position_ >= source_.size()) fail("ended unexpectedly");
        switch (source_[position_]) {
        case '{': return object();
        case '[': return array();
        case '"': return Value(string());
        case 't': literal("true"); return Value(true);
        case 'f': literal("false"); return Value(false);
        case 'n': literal("null"); return Value();
        default: return Value(number());
        }
    }

    void literal(const char* expected) {
        for (const char* item = expected; *item; ++item, ++position_) {
            if (position_ >= source_.size() || source_[position_] != *item)
                fail("has an invalid literal");
        }
    }

    std::string string() {
        if (!take('"')) fail("expected a string");
        std::string result;
        while (position_ < source_.size()) {
            char current = source_[position_++];
            if (current == '"') return result;
            if (current != '\\') {
                result.push_back(current);
                continue;
            }
            if (position_ >= source_.size()) fail("has an incomplete escape");
            const char escaped = source_[position_++];
            switch (escaped) {
            case '"': case '\\': case '/': result.push_back(escaped); break;
            case 'b': result.push_back('\b'); break;
            case 'f': result.push_back('\f'); break;
            case 'n': result.push_back('\n'); break;
            case 'r': result.push_back('\r'); break;
            case 't': result.push_back('\t'); break;
            case 'u': {
                unsigned code = 0;
                for (int index = 0; index < 4; ++index) {
                    if (position_ >= source_.size()) fail("has an incomplete unicode escape");
                    const char hex = source_[position_++];
                    code <<= 4;
                    if (hex >= '0' && hex <= '9') code += static_cast<unsigned>(hex - '0');
                    else if (hex >= 'a' && hex <= 'f') code += static_cast<unsigned>(hex - 'a' + 10);
                    else if (hex >= 'A' && hex <= 'F') code += static_cast<unsigned>(hex - 'A' + 10);
                    else fail("has an invalid unicode escape");
                }
                appendUtf8(result, code);
                break;
            }
            default: fail("has an invalid escape");
            }
        }
        fail("has an unterminated string");
    }

    double number() {
        whitespace();
        const std::size_t start = position_;
        if (position_ < source_.size() && source_[position_] == '-') ++position_;
        while (position_ < source_.size() && std::isdigit(static_cast<unsigned char>(source_[position_]))) ++position_;
        if (position_ < source_.size() && source_[position_] == '.') {
            ++position_;
            while (position_ < source_.size() && std::isdigit(static_cast<unsigned char>(source_[position_]))) ++position_;
        }
        if (position_ < source_.size() && (source_[position_] == 'e' || source_[position_] == 'E')) {
            ++position_;
            if (position_ < source_.size() && (source_[position_] == '+' || source_[position_] == '-')) ++position_;
            while (position_ < source_.size() && std::isdigit(static_cast<unsigned char>(source_[position_]))) ++position_;
        }
        if (start == position_) fail("expected a value");
        return std::stod(source_.substr(start, position_ - start));
    }

    Value array() {
        take('[');
        Value::Array result;
        if (take(']')) return Value(std::move(result));
        do { result.push_back(value()); } while (take(','));
        if (!take(']')) fail("expected ]");
        return Value(std::move(result));
    }

    Value object() {
        take('{');
        Value::Object result;
        if (take('}')) return Value(std::move(result));
        do {
            whitespace();
            const std::string key = string();
            if (!take(':')) fail("expected :");
            result.emplace(key, value());
        } while (take(','));
        if (!take('}')) fail("expected }");
        return Value(std::move(result));
    }
};

}  // namespace

bool Value::isNull() const { return std::holds_alternative<std::nullptr_t>(data_); }
bool Value::boolean(bool fallback) const {
    const auto* value = std::get_if<bool>(&data_); return value ? *value : fallback;
}
double Value::number(double fallback) const {
    const auto* value = std::get_if<double>(&data_); return value ? *value : fallback;
}
std::string Value::string(std::string fallback) const {
    const auto* value = std::get_if<std::string>(&data_); return value ? *value : fallback;
}
const Value::Array& Value::array() const {
    const auto* value = std::get_if<Array>(&data_); return value ? *value : kEmptyArray;
}
const Value::Object& Value::object() const {
    const auto* value = std::get_if<Object>(&data_); return value ? *value : kEmptyObject;
}
const Value& Value::at(const std::string& key) const {
    const auto& values = object();
    const auto found = values.find(key); return found == values.end() ? kNull : found->second;
}

Value parse(const std::string& text) { return Parser(text).run(); }

Value load(const std::wstring& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("Unable to open JSON file");
    std::ostringstream buffer; buffer << stream.rdbuf();
    return parse(buffer.str());
}

Value::Array& Value::array() {
    if (!std::holds_alternative<Array>(data_)) data_ = Array{};
    return std::get<Array>(data_);
}

Value::Object& Value::object() {
    if (!std::holds_alternative<Object>(data_)) data_ = Object{};
    return std::get<Object>(data_);
}

Value& Value::at(const std::string& key) { return object()[key]; }

namespace {

std::string escaped(const std::string& text) {
    std::ostringstream out; out << '"';
    for (const unsigned char item : text) {
        switch (item) {
        case '"': out << "\\\""; break; case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break; case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break; case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (item < 0x20) out << "\\u" << std::hex << std::setw(4)
                                 << std::setfill('0') << static_cast<int>(item) << std::dec;
            else out << static_cast<char>(item);
        }
    }
    return out.str() + '"';
}

void encode(std::ostringstream& out, const Value& value) {
    if (value.isNull()) { out << "null"; return; }
    const std::string stringValue = value.string("\x01");
    if (stringValue != "\x01") { out << escaped(stringValue); return; }
    if (!value.array().empty()) {
        out << '[';
        for (std::size_t index = 0; index < value.array().size(); ++index) {
            if (index) out << ','; encode(out, value.array()[index]);
        }
        out << ']'; return;
    }
    if (!value.object().empty()) {
        out << '{'; bool first = true;
        for (const auto& entry : value.object()) {
            if (!first) out << ','; first = false;
            out << escaped(entry.first) << ':'; encode(out, entry.second);
        }
        out << '}'; return;
    }
    // Empty arrays and objects need their actual variant, not their accessor fallback.
    const auto& storage = value.storage();
    if (std::holds_alternative<Value::Array>(storage)) { out << "[]"; return; }
    if (std::holds_alternative<Value::Object>(storage)) { out << "{}"; return; }
    if (const auto* boolean = std::get_if<bool>(&storage)) out << (*boolean ? "true" : "false");
    else if (const auto* number = std::get_if<double>(&storage))
        out << std::setprecision(17) << *number;
    else out << "null";
}

}  // namespace

std::string dump(const Value& value) {
    std::ostringstream out; encode(out, value); return out.str();
}

void setPath(Value& root, const std::string& dottedPath, Value value) {
    if (dottedPath.empty()) throw std::invalid_argument("JSON path is empty");
    const auto index = [](const std::string& component, std::size_t extent) {
        std::size_t result{};
        const auto parsed = std::from_chars(
            component.data(), component.data() + component.size(), result);
        if (component.empty() || parsed.ec != std::errc{}
            || parsed.ptr != component.data() + component.size())
            throw std::invalid_argument(
                "JSON array path component must be a non-negative integer");
        if (result >= extent)
            throw std::out_of_range("JSON array path index is out of range");
        return result;
    };
    const auto child = [&index](Value& parent, const std::string& component) {
        if (std::holds_alternative<Value::Array>(parent.storage())) {
            auto& values = parent.array();
            return &values[index(component, values.size())];
        }
        return &parent.at(component);
    };
    Value* current = &root;
    std::size_t start = 0;
    while (true) {
        const std::size_t dot = dottedPath.find('.', start);
        const std::string key = dottedPath.substr(start, dot - start);
        if (key.empty()) throw std::invalid_argument("JSON path contains an empty component");
        if (dot == std::string::npos) {
            *child(*current, key) = std::move(value);
            return;
        }
        current = child(*current, key);
        start = dot + 1;
    }
}

}  // namespace dynsyslab::core::json
