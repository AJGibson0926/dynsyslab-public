// Loads a job group: reads child job files, binds symbolic inputs and flattens the tree to an
// ordered list of leaf jobs. Nothing runs here, so validate and plan share it with run.
#include <dynsyslab/workflow/JobGroup.hpp>

#include <algorithm>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace dynsyslab::workflow {
namespace {

constexpr std::size_t kMaximumDepth = 16;

struct ChildOutputs {
    std::map<std::string, std::string> named;  // output name -> token
    std::string leaf;                          // set when the child is a single leaf job
    std::set<std::string> leaves;              // every leaf beneath the child
};

struct Scope {
    std::map<std::string, json::Value> inputs;
    json::Value parameters;
    const std::map<std::string, ChildOutputs>* siblings{};
};

json::Value readJobFile(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::invalid_argument("job group: child job file is unavailable: " + path.u8string());
    std::ostringstream text; text << stream.rdbuf();
    try { return json::parse(text.str()); }
    catch (const std::exception& error) {
        throw std::invalid_argument("job group: " + path.u8string() + ": " + error.what());
    }
}

const json::Value* parameter(const json::Value& parameters, const std::string& path) {
    const json::Value* node = &parameters;
    std::size_t begin = 0;
    while (true) {
        const auto dot = path.find('.', begin);
        const auto* fields = std::get_if<json::Value::Object>(&node->storage());
        const auto key = path.substr(begin, dot == std::string::npos ? std::string::npos : dot - begin);
        if (!fields || !fields->count(key)) return nullptr;
        node = &fields->at(key);
        if (dot == std::string::npos) return node;
        begin = dot + 1;
    }
}

std::string siblingToken(const Scope& scope, const std::string& id, const std::string& field) {
    if (!scope.siblings || !scope.siblings->count(id))
        throw std::invalid_argument("job group: reference to job '" + id + "', which does not run earlier");
    const auto& child = scope.siblings->at(id);
    if (field.rfind("outputs.", 0) == 0) {
        const auto name = field.substr(8);
        if (!child.named.count(name))
            throw std::invalid_argument("job group: job '" + id + "' declares no output '" + name + "'");
        return child.named.at(name);
    }
    if (child.leaf.empty())
        throw std::invalid_argument("job group: '" + id + "' is a group; refer to ${" + id + ".outputs.NAME}");
    return "@{" + child.leaf + "|" + field + "}";
}

// Replaces ${inputs.X}, ${params.X} and, when siblings are in scope, ${id.outputs.X} and the
// pipeline-style ${id.archive} forms. A whole-string reference keeps the value's JSON type.
json::Value bind(const json::Value& value, const Scope& scope, bool keepUnknownParameters) {
    static const std::regex pattern(
        R"(\$\{(?:(inputs|params)\.([A-Za-z0-9_.-]+)|([A-Za-z0-9][A-Za-z0-9_-]*)\.(outputs\.[A-Za-z0-9_-]+|archive|archives\[[0-9]+\]|output:[A-Za-z0-9_]+))\})");
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        const auto lookup = [&](const std::smatch& match, bool& known) -> json::Value {
            known = true;
            if (match[1].matched) {
                const auto name = match[2].str();
                if (match[1] == "inputs") {
                    if (!scope.inputs.count(name))
                        throw std::invalid_argument("job group: job uses undeclared input '" + name + "'");
                    return scope.inputs.at(name);
                }
                if (const auto* found = parameter(scope.parameters, name)) return *found;
                if (keepUnknownParameters) { known = false; return json::Value(match[0].str()); }
                throw std::invalid_argument("job group: undefined parameter '" + name + "'");
            }
            if (!scope.siblings) { known = false; return json::Value(match[0].str()); }
            return json::Value(siblingToken(scope, match[3].str(), match[4].str()));
        };
        std::smatch match; bool known = true;
        if (std::regex_match(*text, match, pattern)) return lookup(match, known);
        std::string out; auto begin = text->cbegin();
        while (std::regex_search(begin, text->cend(), match, pattern)) {
            out.append(begin, match[0].first);
            const auto replaced = lookup(match, known);
            const auto* asText = std::get_if<std::string>(&replaced.storage());
            out += asText ? *asText : json::dump(replaced);
            begin = match[0].second;
        }
        out.append(begin, text->cend());
        return json::Value(out);
    }
    if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        json::Value::Array out;
        for (const auto& item : *array) out.push_back(bind(item, scope, keepUnknownParameters));
        return json::Value(std::move(out));
    }
    if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        json::Value::Object out;
        for (const auto& [key, item] : *object) out[key] = bind(item, scope, keepUnknownParameters);
        return json::Value(std::move(out));
    }
    return value;
}

void tokens(const json::Value& value, std::set<std::string>& leaves) {
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        static const std::regex pattern(R"(@\{([^|}]+)\|[^}]+\})");
        for (std::sregex_iterator it(text->begin(), text->end(), pattern), end; it != end; ++it)
            leaves.insert((*it)[1].str());
    } else if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        for (const auto& item : *array) tokens(item, leaves);
    } else if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        for (const auto& [key, item] : *object) { (void)key; tokens(item, leaves); }
    }
}

void siblingReferences(const json::Value& value, std::set<std::string>& ids) {
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        static const std::regex pattern(
            R"(\$\{([A-Za-z0-9][A-Za-z0-9_-]*)\.(?:outputs\.[A-Za-z0-9_-]+|archive|archives\[[0-9]+\]|output:[A-Za-z0-9_]+)\})");
        for (std::sregex_iterator it(text->begin(), text->end(), pattern), end; it != end; ++it)
            if ((*it)[1] != "inputs" && (*it)[1] != "params") ids.insert((*it)[1].str());
    } else if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        for (const auto& item : *array) siblingReferences(item, ids);
    } else if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        for (const auto& [key, item] : *object) { (void)key; siblingReferences(item, ids); }
    }
}

// Stable topological order: a child runs as early as its declared position allows.
std::vector<std::size_t> dependencyOrder(const json::Value::Array& entries,
                                         std::vector<std::set<std::string>>& requires_) {
    std::map<std::string, std::size_t> position;
    for (std::size_t k = 0; k < entries.size(); ++k) {
        const auto id = entries[k].at("id").string();
        static const std::regex valid("[A-Za-z0-9][A-Za-z0-9_-]*");
        if (!std::regex_match(id, valid))
            throw std::invalid_argument("job group: job ids use letters, digits, '-' and '_': '" + id + "'");
        if (!position.emplace(id, k).second) throw std::invalid_argument("job group: duplicate job id: " + id);
    }
    requires_.assign(entries.size(), {});
    for (std::size_t k = 0; k < entries.size(); ++k) {
        for (const auto& item : entries[k].at("depends_on").array()) requires_[k].insert(item.string());
        siblingReferences(entries[k].at("inputs"), requires_[k]);
        siblingReferences(entries[k].at("parameters"), requires_[k]);
        for (const auto& id : requires_[k])
            if (!position.count(id))
                throw std::invalid_argument("job group: job '" + entries[k].at("id").string()
                                            + "' depends on unknown job '" + id + "'");
    }
    std::vector<std::size_t> order; std::set<std::string> placed;
    while (order.size() < entries.size()) {
        bool progressed = false;
        for (std::size_t k = 0; k < entries.size() && !progressed; ++k) {
            const auto id = entries[k].at("id").string();
            if (placed.count(id)) continue;
            if (std::all_of(requires_[k].begin(), requires_[k].end(),
                            [&](const std::string& item) { return placed.count(item) > 0; })) {
                order.push_back(k); placed.insert(id); progressed = true;
            }
        }
        if (!progressed) {
            std::string names;
            for (const auto& entry : entries)
                if (!placed.count(entry.at("id").string())) names += " " + entry.at("id").string();
            throw std::invalid_argument("job group: dependency cycle among:" + names);
        }
    }
    return order;
}

ChildOutputs loadNode(const json::Value& job, const std::filesystem::path& file, const std::string& id,
                      const Scope& scope, std::vector<std::filesystem::path>& stack, bool keep,
                      JobGroupGraph& graph);

ChildOutputs loadGroup(const json::Value& job, const std::filesystem::path& file, const std::string& prefix,
                       const Scope& scope, std::vector<std::filesystem::path>& stack, bool keep,
                       JobGroupGraph& graph) {
    if (stack.size() > kMaximumDepth)
        throw std::invalid_argument("job group: nesting deeper than " + std::to_string(kMaximumDepth) + " levels");
    const auto& entries = job.at("jobs").array();
    if (entries.empty()) throw std::invalid_argument("job group lists no jobs: " + file.u8string());
    Scope own;
    own.inputs = scope.inputs;
    own.parameters = job.at("parameters").isNull() ? json::Value(json::Value::Object{})
                                                   : bind(job.at("parameters"), scope, false);
    for (const auto& [key, value] : scope.parameters.object()) own.parameters.at(key) = value;  // caller wins
    std::map<std::string, ChildOutputs> done;
    own.siblings = &done;
    std::vector<std::set<std::string>> required;
    ChildOutputs result;
    for (const auto index : dependencyOrder(entries, required)) {
        const auto& entry = entries[index];
        const auto entryId = entry.at("id").string();
        Scope child;
        for (const auto& [name, value] : entry.at("inputs").object()) child.inputs[name] = bind(value, own, false);
        child.parameters = own.parameters;
        for (const auto& [name, value] : entry.at("parameters").object())
            child.parameters.at(name) = bind(value, own, false);
        const auto path = prefix.empty() ? entryId : prefix + "/" + entryId;
        const auto first = graph.leaves.size();
        const bool keepChild = keep || entry.at("keep").boolean(false);
        ChildOutputs loaded;
        if (const auto* name = std::get_if<std::string>(&entry.at("job").storage())) {
            const auto childFile = std::filesystem::weakly_canonical(
                file.parent_path() / std::filesystem::u8path(*name));
            if (std::find(stack.begin(), stack.end(), childFile) != stack.end())
                throw std::invalid_argument("job group: job file includes itself: " + childFile.u8string());
            stack.push_back(childFile);
            loaded = loadNode(readJobFile(childFile), childFile, path, child, stack, keepChild, graph);
            stack.pop_back();
        } else if (std::holds_alternative<json::Value::Object>(entry.at("job").storage())) {
            stack.push_back(file / std::filesystem::u8path(entryId));  // inline jobs count toward the depth limit
            loaded = loadNode(entry.at("job"), file, path, child, stack, keepChild, graph);
            stack.pop_back();
        } else throw std::invalid_argument("job group: job '" + path + "' needs a file path or job object");
        // Explicit and input dependencies apply to every leaf of the child.
        std::set<std::string> upstream;
        for (const auto& id : required[index]) upstream.insert(done.at(id).leaves.begin(), done.at(id).leaves.end());
        for (std::size_t k = first; k < graph.leaves.size(); ++k)
            graph.leaves[k].dependencies.insert(upstream.begin(), upstream.end());
        result.leaves.insert(loaded.leaves.begin(), loaded.leaves.end());
        done[entryId] = std::move(loaded);
    }
    for (const auto& [name, value] : job.at("outputs").object()) {
        const auto bound = bind(value, own, false);
        if (bound.string().rfind("@{", 0) != 0)
            throw std::invalid_argument("job group: output '" + name + "' must name a child output");
        result.named[name] = bound.string();
    }
    return result;
}

ChildOutputs loadNode(const json::Value& job, const std::filesystem::path& file, const std::string& id,
                      const Scope& scope, std::vector<std::filesystem::path>& stack, bool keep,
                      JobGroupGraph& graph) {
    if (job.at("schema").string("dynsyslab.compute-job/v1") != "dynsyslab.compute-job/v1")
        throw std::invalid_argument("job group: " + file.u8string() + " is not a dynsyslab.compute-job/v1 job");
    if (job.at("kind").string() == "job_group") return loadGroup(job, file, id, scope, stack, keep, graph);
    const bool pipeline = job.at("kind").string() == "pipeline";
    Scope leafScope{scope.inputs, scope.parameters, nullptr};
    auto bound = bind(job, leafScope, pipeline);
    if (bound.at("schema").isNull()) bound.at("schema") = json::Value("dynsyslab.compute-job/v1");
    if (bound.at("label").isNull()) bound.at("label") = json::Value(id);
    JobGroupLeaf leaf;
    leaf.id = id; leaf.file = file; leaf.keep = keep || bound.at("keep").boolean(false);
    leaf.estimatedBytes = bound.at("resources").at("peak_bytes").number(0.0);
    ChildOutputs outputs;
    outputs.leaf = id; outputs.leaves.insert(id);
    for (const auto& [name, selector] : bound.at("outputs").object())
        outputs.named[name] = "@{" + id + "|" + selector.string() + "}";
    if (outputs.named.empty()) outputs.named["archive"] = "@{" + id + "|archive}";
    for (const auto* field : {"outputs", "keep", "resources"}) bound.object().erase(field);
    tokens(bound, leaf.dependencies);
    leaf.dependencies.erase(id);
    leaf.job = std::move(bound);
    if (std::any_of(graph.leaves.begin(), graph.leaves.end(),
                    [&](const JobGroupLeaf& item) { return item.id == id; }))
        throw std::invalid_argument("job group: duplicate leaf id: " + id);
    graph.leaves.push_back(std::move(leaf));
    return outputs;
}

}  // namespace

json::Value bindJobFile(json::Value job, const std::filesystem::path& file) {
    if (job.at("kind").string() == "job_group")
        job.at("job_file") = json::Value(std::filesystem::weakly_canonical(file).generic_u8string());
    return job;
}

JobGroupGraph loadJobGroup(const json::Value& job) {
    if (job.at("kind").string() != "job_group") throw std::invalid_argument("loadJobGroup expects a job_group");
    const auto recorded = job.at("job_file").string();
    // An inline group without a file resolves children against the working directory.
    const auto file = recorded.empty() ? std::filesystem::current_path() / "inline.job.json"
                                       : std::filesystem::u8path(recorded);
    JobGroupGraph graph;
    std::vector<std::filesystem::path> stack{file};
    Scope root;
    root.parameters = json::Value(json::Value::Object{});
    const auto outputs = loadGroup(job, file, "", root, stack, false, graph);
    graph.outputs = outputs.named;
    return graph;
}

}  // namespace dynsyslab::workflow
