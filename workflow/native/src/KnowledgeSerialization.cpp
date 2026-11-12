#include <algorithm>
#include <cctype>
#include <cmath>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
json::Value strings(const std::vector<std::string> &values) {
    json::Value::Array result;
    for (const auto &value : values)
        result.emplace_back(value);
    return result;
}
std::string text(const json::Value &value, const char *key, bool required = false) {
    const auto &field = value.at(key);
    if (!field.isNull() && !std::holds_alternative<std::string>(field.storage()))
        throw std::invalid_argument(std::string(key) + " must be text");
    auto result = field.string();
    if ((required && result.find_first_not_of(" \t\r\n") == std::string::npos) ||
        result.size() > 16384)
        throw std::invalid_argument(std::string(key) + " is empty or too long");
    return result;
}
std::vector<std::string> stringList(const json::Value &value, const char *key) {
    const auto &field = value.at(key);
    if (field.isNull())
        return {};
    if (!std::holds_alternative<json::Value::Array>(field.storage()) || field.array().size() > 256)
        throw std::invalid_argument(std::string(key) + " must be a bounded text array");
    std::vector<std::string> result;
    for (const auto &item : field.array()) {
        if (!std::holds_alternative<std::string>(item.storage()) || item.string().empty() ||
            item.string().size() > 16384)
            throw std::invalid_argument(std::string(key) + " contains invalid text");
        result.push_back(item.string());
    }
    return result;
}
} // namespace
json::Value researchPaperJson(const dynsyslab::core::ResearchPaper &p) {
    return json::Value::Object{{"id", p.id},
                               {"title", p.title},
                               {"year", double(p.year)},
                               {"authors", strings(p.authors)},
                               {"venue", p.venue},
                               {"doi", p.doi},
                               {"urls", strings(p.urls)},
                               {"topics", strings(p.topics)},
                               {"related_fact_ids", strings(p.relatedFactIds)},
                               {"findings", p.findings},
                               {"applicability", p.applicability},
                               {"coordinates", p.coordinates},
                               {"verification", p.verification},
                               {"notes", p.notes}};
}
dynsyslab::core::ResearchPaper researchPaperFromJson(const json::Value &value) {
    dynsyslab::core::ResearchPaper p;
    p.id = text(value, "id", true);
    p.title = text(value, "title", true);
    if (p.id.size() > 128 || !std::all_of(p.id.begin(), p.id.end(), [](unsigned char c) {
            return std::isalnum(c) || c == '-' || c == '_' || c == '.';
        }))
        throw std::invalid_argument(
            "Paper id must contain 1..128 letters, digits, dots, underscores or hyphens");
    const auto &year = value.at("year");
    if (!year.isNull() &&
        (!std::holds_alternative<double>(year.storage()) || !std::isfinite(year.number()) ||
         year.number() < 0 || year.number() > 9999 || std::floor(year.number()) != year.number()))
        throw std::invalid_argument("Paper year must be an integer in 0..9999; 0 means unknown");
    p.year = int(year.number());
    p.venue = text(value, "venue");
    p.doi = text(value, "doi");
    for (const auto &prefix : {"https://doi.org/", "http://doi.org/", "doi:"})
        if (p.doi.rfind(prefix, 0) == 0) {
            p.doi.erase(0, std::string(prefix).size());
            break;
        }
    std::transform(p.doi.begin(), p.doi.end(), p.doi.begin(),
                   [](unsigned char c) { return char(std::tolower(c)); });
    if (!p.doi.empty() && (p.doi.rfind("10.", 0) != 0 || p.doi.find('/') == std::string::npos ||
                           p.doi.find_first_of(" \t\r\n") != std::string::npos))
        throw std::invalid_argument("Invalid DOI syntax");
    p.authors = stringList(value, "authors");
    p.urls = stringList(value, "urls");
    for (const auto &url : p.urls)
        if (url.rfind("https://", 0) != 0 && url.rfind("http://", 0) != 0)
            throw std::invalid_argument("Paper URLs must use http or https");
    p.topics = stringList(value, "topics");
    p.relatedFactIds = stringList(value, "related_fact_ids");
    p.findings = text(value, "findings");
    p.applicability = text(value, "applicability");
    p.coordinates = text(value, "coordinates");
    p.verification = text(value, "verification");
    p.notes = text(value, "notes");
    if (p.verification.empty())
        p.verification =
            "User supplied; scientific applicability has not been independently checked";
    return p;
}
json::Value dynamicalKnowledgeJson(const dynsyslab::core::DynamicalKnowledge &k) {
    json::Value::Array facts, papers;
    for (const auto &f : k.facts)
        facts.emplace_back(json::Value::Object{{"id", f.id},
                                               {"title", f.title},
                                               {"kind", f.kind},
                                               {"statement", f.statement},
                                               {"formula_latex", f.formula},
                                               {"conditions", f.conditions},
                                               {"coordinates", f.coordinates},
                                               {"references", strings(f.references)},
                                               {"applicable", f.applicable},
                                               {"applicability_reason", f.applicabilityReason},
                                               {"evidence", f.evidence}});
    for (const auto &p : k.papers)
        papers.push_back(researchPaperJson(p));
    return json::Value::Object{
        {"schema", k.schema},
        {"version", k.version},
        {"research_status", k.researchStatus},
        {"facts", facts},
        {"paper_catalogue",
         json::Value::Object{{"schema", "dynsyslab.paper-catalogue/v1"}, {"papers", papers}}}};
}
} // namespace dynsyslab::workflow
