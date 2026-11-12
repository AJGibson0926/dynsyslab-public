#pragma once
#include <dynsyslab/core/Observable.hpp>
#include <optional>
#include <algorithm>
#include <set>

namespace dynsyslab::core::jet {

inline void readAxes(const json::Value& node,std::set<std::size_t>& used) {
    if (node.at("operation").string()=="state") used.insert(std::size_t(node.at("state_index").number()));
    for (const auto& arg:node.at("arguments").array()) readAxes(arg,used);
}

// If the measured axes form a closed subsystem, varying an excluded axis
// inside the open regular domain preserves every measured derivative.
inline json::Value hiddenSubsystem(const std::vector<json::Value>& field,const std::vector<json::Value>& h) {
    std::set<std::size_t> seen;
    for (const auto& t:h) readAxes(t,seen);
    std::vector<std::size_t> queue(seen.begin(),seen.end());
    for (std::size_t j=0;j<queue.size();++j) {
        std::set<std::size_t> depends;
        readAxes(field.at(queue[j]),depends);
        for (auto i:depends) if (seen.insert(i).second) queue.push_back(i);
    }
    if (seen.size()==field.size()) return {};
    json::Value::Array hidden;
    for (std::size_t i=0;i<field.size();++i) if (!seen.count(i)) hidden.emplace_back(double(i));
    return json::Value::Object{{"method","closed_subsystem_dependencies"},{"hidden_state_axes",std::move(hidden)},
        {"state_recovery","impossible_at_every_order_on_nonempty_open_domain"}};
}

// Rename source axes in analytic rate expressions. Copy the AST verbatim
// apart from axis names/dimension so no parameter arithmetic is rounded.
inline std::optional<json::Value> onAxes(json::Value tree,const std::vector<std::size_t>& axes) {
    tree.at("state_dimension")=double(axes.size());
    if (tree.at("operation").string()=="state") {
        const auto i=std::size_t(tree.at("state_index").number());
        const auto it=std::find(axes.begin(),axes.end(),i);
        if (it==axes.end()) return {};
        tree.at("state_index")=double(it-axes.begin());
    } else if (tree.object().count("arguments")) {
        for (auto& arg:tree.at("arguments").array()) {
            auto renamed=onAxes(arg,axes);
            if (!renamed) return {};
            arg=std::move(*renamed);
        }
    }
    return tree;
}

inline json::Value derivativeCoordinateChart(const std::vector<json::Value>& field,
    std::vector<json::Value> level,std::size_t maxOrder) {
    std::vector<std::size_t> axes;
    const auto width=level.size();
    for (std::size_t k=0;k<=maxOrder;++k) {
        std::vector<json::Value> next;
        for (const auto& node:level) {
            if (node.at("operation").string()!="state") return {};
            const auto i=std::size_t(node.at("state_index").number());
            axes.push_back(i);
            next.push_back(field.at(i));
        }
        json::Value::Array rhs,coordinateAxes;
        for (auto i:axes) coordinateAxes.emplace_back(double(i));
        for (std::size_t i=width;i<axes.size();++i)
            rhs.push_back(AnalyticExpression::state(axes.size(),i).toJson());
        bool closed=true;
        for (const auto& p:next) {
            auto renamed=onAxes(p,axes);
            if (!renamed) { closed=false; break; }
            rhs.push_back(std::move(*renamed));
        }
        if (closed) {
            const std::set<std::size_t> distinct(axes.begin(),axes.end());
            return json::Value::Object{{"order",double(k)},{"coordinate_count",double(axes.size())},
                {"status","certified_closed"},{"method","exact_state_axis_substitution"},
                {"ode_regularity","smooth near each realizable jet"},
                {"jet_state_axes",coordinateAxes},{"state_recovery",distinct.size()==field.size()?
                    "certified_inverse":"impossible_on_open_state_space_by_closed_projection"},
                {"rhs_observable",json::Value::Object{{"schema","dynsyslab.observable-expression/v1"},
                    {"name","jet_rhs"},{"expressions",std::move(rhs)}}}};
        }
        level=std::move(next);
    }
    return {};
}
} // namespace dynsyslab::core::jet
