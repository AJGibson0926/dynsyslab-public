// Sweep and select experiment stages (ExperimentStages.hpp, EXPERIMENTS.md).
#include "ExperimentStages.hpp"
#include "ExperimentJson.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>

namespace dynsyslab::workflow::experiment_detail {
namespace {

using J = json::Value;

// ${sweep.NAME} as the whole string is the point's JSON value; embedded, its text.
// ${sweep.NAME.FIELD...} reads a field of an object-valued variable (one point = one set-up).
J substitute(const J& value, const J::Object& point, const std::map<std::string, std::string>& renamed) {
    static const std::regex variable(R"(\$\{sweep\.([A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*)\})");
    static const std::regex stage(R"(\$\{([A-Za-z][A-Za-z0-9_-]*)\.)");
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        const auto lookup = [&](const std::string& path) -> const J& {
            const auto dot = path.find('.');
            const auto found = point.find(path.substr(0, dot));
            if (found == point.end()) throw std::invalid_argument("undefined sweep variable: " + path);
            const J* current = &found->second;
            for (auto begin = dot; begin != std::string::npos;) {
                const auto end = path.find('.', begin + 1);
                const auto key = path.substr(begin + 1, end == std::string::npos ? std::string::npos : end - begin - 1);
                const auto* fields = std::get_if<J::Object>(&current->storage());
                if (!fields || !fields->count(key)) throw std::invalid_argument("undefined sweep variable field: " + path);
                current = &fields->at(key);
                begin = end;
            }
            return *current;
        };
        std::smatch match;
        if (std::regex_match(*text, match, variable)) return lookup(match[1].str());
        std::string out;
        auto begin = text->cbegin();
        while (std::regex_search(begin, text->cend(), match, variable)) {
            out.append(begin, match[0].first);
            const auto& v = lookup(match[1].str());
            out += std::holds_alternative<std::string>(v.storage()) ? v.string() : json::dump(v);
            begin = match[0].second;
        }
        out.append(begin, text->cend());
        std::string renamedText;  // references to the block's own stages
        begin = out.cbegin();
        while (std::regex_search(begin, out.cend(), match, stage)) {
            renamedText.append(begin, match[0].first);
            const auto found = renamed.find(match[1].str());
            renamedText += "${" + (found == renamed.end() ? match[1].str() : found->second) + ".";
            begin = match[0].second;
        }
        renamedText.append(begin, out.cend());
        return J(renamedText);
    }
    if (const auto* array = std::get_if<J::Array>(&value.storage())) {
        J::Array out;
        for (const auto& item : *array) out.push_back(substitute(item, point, renamed));
        return J(std::move(out));
    }
    if (const auto* object = std::get_if<J::Object>(&value.storage())) {
        J::Object out;
        for (const auto& [key, item] : *object) out[key] = substitute(item, point, renamed);
        return J(std::move(out));
    }
    return value;
}

std::vector<J::Object> points(const J& spec, const std::string& id) {
    const auto& variables = spec.at("variables").object();
    if (variables.empty()) throw std::invalid_argument("sweep '" + id + "' names no variables");
    const auto mode = spec.at("mode").string("product");
    std::vector<std::pair<std::string, J::Array>> axes;
    for (const auto& [name, values] : variables) {
        if (values.array().empty()) throw std::invalid_argument("sweep '" + id + "' variable " + name + " has no values");
        axes.emplace_back(name, values.array());
    }
    // Object keys are sorted; "order" keeps the author's slowest-to-fastest order when given.
    if (!spec.at("order").isNull()) {
        std::vector<std::pair<std::string, J::Array>> ordered;
        for (const auto& name : spec.at("order").array()) {
            const auto at = std::find_if(axes.begin(), axes.end(), [&](const auto& a) { return a.first == name.string(); });
            if (at == axes.end()) throw std::invalid_argument("sweep order names unknown variable " + name.string());
            ordered.push_back(*at);
        }
        if (ordered.size() != axes.size()) throw std::invalid_argument("sweep order must list every variable");
        axes = std::move(ordered);
    }
    std::vector<J::Object> out;
    if (mode == "zip") {
        const auto n = axes.front().second.size();
        for (const auto& a : axes) if (a.second.size() != n) throw std::invalid_argument("zip sweep variables differ in length");
        for (std::size_t k = 0; k < n; ++k) {
            J::Object p;
            for (const auto& a : axes) p[a.first] = a.second[k];
            out.push_back(std::move(p));
        }
    } else if (mode == "product") {
        std::size_t total = 1;
        for (const auto& a : axes) total *= a.second.size();
        if (total > 10000) throw std::invalid_argument("sweep '" + id + "' has more than 10000 points");
        for (std::size_t k = 0; k < total; ++k) {
            J::Object p;
            std::size_t rest = k;
            for (std::size_t a = axes.size(); a-- > 0;) {  // the last variable varies fastest
                p[axes[a].first] = axes[a].second[rest % axes[a].second.size()];
                rest /= axes[a].second.size();
            }
            out.push_back(std::move(p));
        }
    } else {
        throw std::invalid_argument("sweep mode is product or zip");
    }
    return out;
}

// Follows a dotted metadata path (numeric segments index arrays).
const J& follow(const J& root, const std::string& path) {
    const J* node = &root;
    std::size_t begin = 0;
    while (begin <= path.size()) {
        const auto dot = path.find('.', begin);
        const auto segment = path.substr(begin, dot == std::string::npos ? std::string::npos : dot - begin);
        if (std::holds_alternative<J::Array>(node->storage())) {
            if (segment.find_first_not_of("0123456789") != std::string::npos || std::stoul(segment) >= node->array().size())
                throw std::invalid_argument("select metric path " + path + " indexes an array with " + segment);
            node = &node->array()[std::stoul(segment)];
        } else {
            node = &node->at(segment);
        }
        if (node->isNull()) throw std::invalid_argument("select metric path " + path + " is missing at " + segment);
        if (dot == std::string::npos) break;
        begin = dot + 1;
    }
    return *node;
}

void numbers(const J& value, std::vector<double>& out) {
    if (const auto* x = std::get_if<double>(&value.storage())) { out.push_back(*x); return; }
    for (const auto& item : value.array()) numbers(item, out);
}

double reduce(std::vector<double> v, const std::string& how) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (how == "first") return v.front();
    if (how == "last") return v.back();
    if (how == "max") return *std::max_element(v.begin(), v.end());
    if (how == "min") return *std::min_element(v.begin(), v.end());
    if (how == "mean") { double s = 0; for (const auto x : v) s += x; return s / double(v.size()); }
    if (how == "rms") { double s = 0; for (const auto x : v) s += x * x; return std::sqrt(s / double(v.size())); }
    throw std::invalid_argument("select reduce is mean, rms, max, min, first or last");
}

// Validation-only scores: the metric's section must hold held-out results.
bool validationSection(const std::string& path) {
    const auto head = path.substr(0, path.find('.'));
    return head.find("validation") != std::string::npos || head.find("held_out") != std::string::npos;
}

}  // namespace

void expandSweepsAndSelects(J& expanded) {
    J::Array out;
    auto aliases = expanded.at("stage_aliases").object();
    auto sweeps = expanded.at("sweeps").object();
    const auto label = expanded.at("label").string("experiment");
    for (const auto& stage : expanded.at("stages").array()) {
        const auto id = stage.at("id").string();
        if (!stage.at("sweep").isNull()) {
            const auto& spec = stage.at("sweep");
            const auto& block = spec.at("stages").array();
            if (block.empty()) throw std::invalid_argument("sweep '" + id + "' has no stages");
            std::set<std::string> inner;
            for (const auto& s : block) {
                if (!inner.insert(s.at("id").string()).second || s.at("id").string().empty())
                    throw std::invalid_argument("sweep '" + id + "' stages need distinct ids");
                if (!isObject(s.at("job"))) throw std::invalid_argument("sweep '" + id + "' stages need job objects");
            }
            J::Array recorded, generated;
            const auto all = points(spec, id);
            for (std::size_t k = 0; k < all.size(); ++k) {
                const auto prefix = id + "_" + std::to_string(k + 1) + "_";
                std::map<std::string, std::string> renamed;
                for (const auto& name : inner) renamed[name] = prefix + name;
                J::Object members;
                for (const auto& s : block) {
                    auto job = substitute(s.at("job"), all[k], renamed);
                    if (job.at("schema").isNull()) job.at("schema") = J("dynsyslab.compute-job/v1");
                    if (job.at("label").isNull()) job.at("label") = J(label + " / " + prefix + s.at("id").string());
                    out.push_back(J::Object{{"id", J(prefix + s.at("id").string())}, {"job", std::move(job)}});
                    members[s.at("id").string()] = J(prefix + s.at("id").string());
                    generated.emplace_back(prefix + s.at("id").string());
                }
                recorded.push_back(J::Object{{"index", double(k + 1)}, {"values", J(all[k])}, {"stages", J(std::move(members))}});
            }
            sweeps[id] = J::Object{{"points", J(std::move(recorded))}};
            aliases[id] = J(std::move(generated));
            continue;
        }
        if (!stage.at("select").isNull()) {
            const auto& spec = stage.at("select");
            const auto metric = spec.at("metric").string();
            if (metric.empty()) throw std::invalid_argument("select '" + id + "' needs a metric");
            const auto arrayMetric = metric.rfind("array:", 0) == 0;
            if (!arrayMetric && !validationSection(metric) && !spec.at("training_score").boolean(false))
                throw std::invalid_argument("select '" + id + "' must score a validation section (metric path under a "
                                            "'validation' key); training scores need \"training_score\": true");
            J::Array candidates;
            const auto from = spec.at("sweep").string();
            const auto object = spec.at("object").string();
            const auto reference = [&](const std::string& stageId) {
                return J("${" + stageId + (object.empty() ? ".archive}" : ".object:" + object + "}"));
            };
            if (!from.empty()) {
                const auto& sweep = sweeps[from];
                if (sweep.isNull()) throw std::invalid_argument("select '" + id + "' names unknown or later sweep " + from);
                const auto inner = spec.at("stage").string(), standsFor = spec.at("stands_for").string();
                for (const auto& point : sweep.at("points").array()) {
                    const auto member = point.at("stages").at(inner).string();
                    if (member.empty()) throw std::invalid_argument("select '" + id + "': sweep " + from + " has no stage " + inner);
                    J::Object candidate{{"stage", J(member)}, {"values", point.at("values")}, {"archive", reference(member)}};
                    if (!standsFor.empty()) {  // the winner's archive of another stage of its point
                        const auto other = point.at("stages").at(standsFor).string();
                        if (other.empty()) throw std::invalid_argument("select '" + id + "': sweep " + from + " has no stage " + standsFor);
                        candidate["result"] = J("${" + other + ".archive}");
                    }
                    candidates.push_back(J(std::move(candidate)));
                }
            } else {
                for (const auto& member : spec.at("candidates").array())
                    candidates.push_back(J::Object{{"stage", member}, {"values", J()}, {"archive", reference(member.string())}});
            }
            if (candidates.empty()) throw std::invalid_argument("select '" + id + "' has no candidates");
            J job = J::Object{{"schema", J("dynsyslab.compute-job/v1")}, {"kind", J("select")}, {"label", J(label + " / " + id)},
                {"metric", J(metric)}, {"goal", J(spec.at("goal").string("minimize"))}, {"reduce", J(spec.at("reduce").string("mean"))},
                {"training_score", J(spec.at("training_score").boolean(false))}, {"candidates", J(std::move(candidates))}};
            out.push_back(J::Object{{"id", J(id)}, {"job", std::move(job)}});
            continue;
        }
        out.push_back(stage);
    }
    expanded.at("stages") = J(std::move(out));
    expanded.at("stage_aliases") = J(std::move(aliases));
    if (!sweeps.empty()) expanded.at("sweeps") = J(std::move(sweeps));
}

RunResult runSelect(const J& job, const std::filesystem::path& recordPath) {
    const auto metric = job.at("metric").string(), goal = job.at("goal").string("minimize"), how = job.at("reduce").string("mean");
    if (goal != "minimize" && goal != "maximize") throw std::invalid_argument("select goal is minimize or maximize");
    const bool arrayMetric = metric.rfind("array:", 0) == 0;
    J::Array scored;
    std::size_t best = std::numeric_limits<std::size_t>::max();
    double bestScore = 0;
    const auto& candidates = job.at("candidates").array();
    for (std::size_t k = 0; k < candidates.size(); ++k) {
        const auto path = std::filesystem::u8path(candidates[k].at("archive").string());
        const auto archive = dynsyslab::core::readNpz(path);
        std::vector<double> values;
        std::string reason;
        try {
            if (arrayMetric) values = archive.at(metric.substr(6)).doubles();
            else numbers(follow(json::parse(archive.metadataJson()), metric), values);
        } catch (const std::exception& error) { reason = error.what(); }
        const double score = reason.empty() ? reduce(values, how) : std::numeric_limits<double>::quiet_NaN();
        if (reason.empty() && !std::isfinite(score)) reason = "the score is not finite";
        J::Object entry{{"stage", candidates[k].at("stage")}, {"values", candidates[k].at("values")},
                        {"archive", J(path.generic_u8string())}, {"score", std::isfinite(score) ? J(score) : J()}};
        if (!candidates[k].at("result").isNull()) entry["result"] = candidates[k].at("result");
        if (!reason.empty()) entry["rejected"] = J(reason);
        scored.push_back(J(std::move(entry)));
        if (!reason.empty()) continue;
        if (best == std::numeric_limits<std::size_t>::max() || (goal == "minimize" ? score < bestScore : score > bestScore)) {
            best = k; bestScore = score;  // ties keep the earlier point
        }
    }
    if (best == std::numeric_limits<std::size_t>::max()) throw std::runtime_error("select: no candidate has a finite score");
    for (std::size_t k = 0; k < scored.size(); ++k)
        if (k != best && scored[k].at("rejected").isNull()) scored[k].object()["rejected"] = J(std::string("outscored"));
    const J record = J::Object{{"schema", J("dynsyslab.experiment-selection/v1")}, {"label", job.at("label")},
        {"metric", J(metric)}, {"goal", J(goal)}, {"reduce", J(how)}, {"validation_only", J(!job.at("training_score").boolean(false))},
        {"winner", scored[best]}, {"candidates", J(std::move(scored))}};
    (void)dynsyslab::core::prepareOutputFile(recordPath, "experiment selection");
    { std::ofstream stream(recordPath, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to write the selection record");
      stream << json::dump(record) << '\n'; }
    RunResult run{"select", L"Selected " + dynsyslab::core::wide(record.at("winner").at("stage").string()), {}, recordPath, {}};
    const auto& winner = record.at("winner");
    run.archives.push_back(std::filesystem::u8path((winner.at("result").isNull() ? winner.at("archive") : winner.at("result")).string()));
    run.methods.push_back("select");
    return run;
}

}  // namespace dynsyslab::workflow::experiment_detail
