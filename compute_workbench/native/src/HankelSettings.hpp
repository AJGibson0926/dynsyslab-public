#pragma once
#include <dynsyslab/core/Json.hpp>
#include <string>

namespace dsl::compute {
namespace hankel_ui {
using Value=dynsyslab::core::json::Value;

inline std::string channelText(const Value& settings,const std::string& fallback) {
    const auto& channels=settings.at("observables");
    if (!channels.array().empty()) return dynsyslab::core::json::dump(channels);
    const auto& one=settings.at("observable");
    if (std::holds_alternative<double>(one.storage())) return "["+dynsyslab::core::json::dump(one)+"]";
    return one.string(fallback);
}

inline Value channels(const std::string& text) {
    const auto first=text.find_first_not_of(" \t\r\n");
    if (first==std::string::npos) return Value::Array{};
    if (text[first]=='[') return dynsyslab::core::json::parse(text.substr(first));
    Value::Array result;
    for (std::size_t start=first;start<text.size();) {
        const auto end=text.find(',',start);
        const auto field=text.substr(start,end-start);
        const auto a=field.find_first_not_of(" \t"),b=field.find_last_not_of(" \t");
        if (a!=std::string::npos) result.emplace_back(field.substr(a,b-a+1));
        if (end==std::string::npos) break;
        start=end+1;
    }
    return result;
}

// Editing a previously saved common-window job preserves fields not represented
// by this dialog and writes the visible geometry back in its original convention.
inline std::string preserveWindow(std::string settings,bool dmd) {
    auto obj=dynsyslab::core::json::parse(settings);
    auto& fields=obj.object();
    if (!fields.count("window") || fields.at("window").isNull()) return settings;
    auto& w=fields.at("window").object();
    auto transfer=[&](const char* key,const char* target) {
        auto it=fields.find(key);
        if (it!=fields.end()) { w[target]=it->second; fields.erase(it); }
    };
    transfer("matrix_size",dmd?"columns":"delays");
    transfer("sampling_window","delays");
    transfer("sample_stride","sample_stride"); transfer("snapshot_stride","snapshot_stride");
    transfer("delay_stride","delay_stride");
    return dynsyslab::core::json::dump(obj);
}
}
} // namespace dsl::compute
