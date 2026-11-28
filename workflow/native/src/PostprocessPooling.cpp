// "pool": true on a fitting stage (PostprocessInternal.hpp): a post-processing job's sources become one
// pooled fit instead of one fit per source.
#include "PostprocessInternal.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {

namespace {
// The fits that read further records of the same system from settings "additional_sources".
bool readsAdditionalSources(const std::string& method) {
    static const std::set<std::string> methods{
        "hdmd", "sindy", "eigenvalue_optimization", "koopman_spectral", "resdmd", "mpedmd", "hermitian_dmd",
        "tlsdmd", "multdmd", "stochastic_resdmd"};
    return methods.count(method) != 0;
}
}  // namespace

std::optional<json::Value> pooledSourcesJob(const json::Value& job) {
    const auto& stages = job.at("methods").array();
    const bool pooled = std::any_of(stages.begin(), stages.end(), [](const json::Value& stage) {
        return stage.at("enabled").boolean(true) && stage.at("settings").at("pool").boolean(false);
    });
    if (!pooled) return std::nullopt;
    std::size_t enabled = 0;
    for (const auto& stage : stages) enabled += stage.at("enabled").boolean(true) ? 1 : 0;
    if (enabled != 1) throw std::invalid_argument("\"pool\": true needs a post-processing job of one stage");
    const auto& sources = job.at("sources").array();
    if (sources.empty()) throw std::invalid_argument("post-processing requires a source");
    auto rewritten = job;
    for (auto& stage : rewritten.at("methods").array()) {
        if (!stage.at("enabled").boolean(true)) continue;
        const auto method = stage.at("name").string();
        if (!readsAdditionalSources(method))
            throw std::invalid_argument("\"pool\": true is for fits that pool records (additional_sources); not " + method);
        auto& settings = stage.at("settings");
        json::Value::Array additional(sources.begin() + 1, sources.end());
        for (const auto& extra : settings.at("additional_sources").array()) additional.push_back(extra);
        settings.at("additional_sources") = json::Value(std::move(additional));
        settings.at("pool") = json::Value(false);
    }
    rewritten.at("sources") = json::Value(json::Value::Array{sources.front()});
    return rewritten;
}

}  // namespace dynsyslab::workflow::postprocess_detail
