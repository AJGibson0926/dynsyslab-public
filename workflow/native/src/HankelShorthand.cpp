// Older jobs (and the Solver dialogs) name physics observables directly in a
// Hankel method, e.g. hdmd {"observables": ["radius", "far_field_acoustic_echo"]}.
// The Hankel module only selects columns, so the orchestrator expands that
// shorthand into the explicit form: a coordinate_transform step that evaluates
// the observables, followed by the Hankel method reading that step's columns.
#include "HankelShorthand.hpp"
#include "CoordinateChannels.hpp"
#include <set>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
using J = json::Value;

bool hankelMethod(const std::string& name) {
    static const std::set<std::string> names{"hdmd", "hankel_dmd", "hankel_matrix", "hankel_pod", "hankel_svd",
        "koopman_spectral", "resdmd", "mpedmd", "hermitian_dmd", "tlsdmd", "multdmd", "stochastic_resdmd", "specrkhs"};
    return names.count(name) > 0;
}

J::Array selectors(const J& settings) {
    if (!settings.at("observables").array().empty()) return settings.at("observables").array();
    const auto& single = settings.at("observable");
    if (single.isNull() || (std::holds_alternative<std::string>(single.storage()) && single.string().empty())) return {};
    return {single};
}

bool hasNamedSelector(const J& job) {
    for (const auto& stage : job.at("methods").array()) {
        if (!hankelMethod(stage.at("name").string())) continue;
        for (const auto& item : selectors(stage.at("settings"))) {
            const bool object = std::holds_alternative<J::Object>(item.storage());
            const auto& selector = object ? item.at("observable") : item;
            if (!std::holds_alternative<double>(selector.storage())) return true;
            if (object && item.at("derivative_order").number(0) > 0) return true;
        }
    }
    return false;
}
}  // namespace

J expandHankelCoordinateShorthand(const J& job) {
    if (!hasNamedSelector(job) || job.at("sources").array().empty()) return job;
    const auto columns = timeSeriesColumnNames(loadSource(std::filesystem::u8path(job.at("sources").array().front().string())));
    J result = job;
    J::Array methods;
    std::string current; // channels of the coordinate change now in effect
    for (const auto& stage : job.at("methods").array()) {
        if (stage.at("name").string() == "coordinate_transform") current.clear();
        const auto& settings = stage.at("settings");
        const auto items = selectors(settings);
        bool needed = false;
        if (hankelMethod(stage.at("name").string()) && stage.at("enabled").boolean(true)
            && settings.at("hankel_archive").isNull() && !settings.at("use_previous_hankel").boolean()
            && !settings.at("use_previous_output").boolean())
            for (const auto& item : items) needed = needed || selectorNeedsCoordinateChange(item, columns);
        if (!needed) { methods.push_back(stage); continue; }
        // Channel specs carry everything physical; the Hankel items keep only taps.
        J::Array channels, rewritten;
        for (const auto& item : items) {
            const bool object = std::holds_alternative<J::Object>(item.storage()) && item.at("artifact").isNull();
            J spec = object ? item.at("observable") : item;
            if (object && item.at("derivative_order").number(0) > 0)
                spec = J::Object{{"observable", item.at("observable")}, {"derivative_order", item.at("derivative_order")}};
            const auto key = json::dump(spec);
            std::size_t index = 0;
            while (index < channels.size() && json::dump(channels[index]) != key) ++index;
            if (index == channels.size()) channels.push_back(spec);
            J::Object column{{"observable", double(index)}};
            if (object && !item.at("taps").isNull()) column["taps"] = item.at("taps");
            rewritten.emplace_back(std::move(column));
        }
        // Consecutive Hankel methods over the same channels share one coordinate change.
        if (json::dump(J(channels)) != current) {
            methods.emplace_back(J::Object{{"name", "coordinate_transform"}, {"enabled", true},
                {"settings", J::Object{{"kind", "channels"}, {"channels", channels},
                                        {"expanded_for", stage.at("name")}}}});
            current = json::dump(J(channels));
        }
        J replaced = stage;
        auto& out = replaced.at("settings");
        out.object().erase("observable");
        out.at("observables") = rewritten;
        out.at("use_previous_output") = true;
        methods.push_back(std::move(replaced));
    }
    result.at("methods") = methods;
    return result;
}

} // namespace dynsyslab::workflow::postprocess_detail
