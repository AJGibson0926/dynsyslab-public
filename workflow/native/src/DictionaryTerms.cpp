// Dictionary terms written in a job (DictionaryTerms.hpp).
#include "DictionaryTerms.hpp"
#include <dynsyslab/core/AnalyticExpression.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;
using dictionary_terms_detail::numbers;

namespace dictionary_terms_detail {
std::vector<double> numbers(const J& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number());
    return out;
}
}  // namespace dictionary_terms_detail

c::DictionaryTerm dictionaryTerm(const J& spec, std::size_t dimension, const std::vector<std::string>& variables) {
    c::DictionaryTerm t;
    if (!spec.at("center").isNull()) t.center = numbers(spec.at("center"));
    if (!spec.at("scale").isNull()) t.scale = numbers(spec.at("scale"));
    if (!spec.at("expression").isNull()) {
        t.kind = "expression";
        t.expression = std::make_shared<const c::AnalyticExpression>(c::AnalyticExpression::parse(spec.at("expression").string(), variables));
        return t;
    }
    for (const char* kind : {"monomial", "chebyshev", "laguerre", "hermite", "fourier_hermite"})
        if (!spec.at(kind).isNull()) {
            t.kind = kind;
            for (const auto& n : spec.at(kind).array()) t.orders.push_back(std::size_t(n.number()));
            return t;
        }
    for (const char* kind : {"exponential", "sine", "cosine"})
        if (!spec.at(kind).isNull()) {
            t.kind = kind;
            t.vector = numbers(spec.at(kind));
            t.phase = spec.at("phase").number(0);
            return t;
        }
    for (const char* kind : {"gaussian", "laplacian", "lorentzian", "matern"})
        if (!spec.at(kind).isNull()) {
            const auto& radial = spec.at(kind);
            t.kind = kind;
            t.vector = numbers(radial.at("node"));
            t.width = radial.at("width").number(1);
            t.smoothness = radial.at("smoothness").number(2);
            t.amplitude = radial.at("amplitude").number(1);
            return t;
        }
    if (!spec.at("kind").isNull())  // a full term in the archived form
        return c::ExpansionDictionary::fromJson(J::Object{{"schema", "dynsyslab.expansion-dictionary/v1"},
            {"dimension", double(dimension)}, {"terms", J::Array{spec}}}).terms().front();
    throw std::invalid_argument("unrecognized dictionary term " + json::dump(spec));
}

std::shared_ptr<const c::ExpansionDictionary> jobDictionary(const J& spec, std::size_t dimension,
                                                            const std::vector<std::string>& variables) {
    if (spec.isNull()) return nullptr;
    if (std::holds_alternative<J::Object>(spec.storage())) {
        if (!spec.at("schema").isNull())  // an archived dictionary (ExpansionDictionary::describe)
            return std::make_shared<const c::ExpansionDictionary>(c::ExpansionDictionary::fromJson(spec));
        // {"kind": "gaussian" | "laplacian" | "lorentzian" | "matern", "centers": [[...], ...], "width": w}
        const auto kind = spec.at("kind").string("gaussian");
        std::vector<c::DictionaryTerm> terms;
        if (spec.at("include_constant").boolean(false)) {
            c::DictionaryTerm one;
            one.kind = "monomial";
            one.orders.assign(dimension, 0);
            terms.push_back(std::move(one));
        }
        for (const auto& center : spec.at("centers").array()) {
            c::DictionaryTerm t;
            t.kind = kind;
            t.vector = numbers(center);
            if (t.vector.size() != dimension) throw std::invalid_argument("a dictionary center needs one value per coordinate");
            t.width = spec.at("width").number(1);
            t.smoothness = spec.at("smoothness").number(2);
            terms.push_back(std::move(t));
        }
        if (terms.empty()) throw std::invalid_argument("a radial dictionary needs centers");
        return std::make_shared<const c::ExpansionDictionary>(dimension, std::move(terms));
    }
    std::vector<c::DictionaryTerm> terms;
    for (const auto& item : spec.array()) terms.push_back(dictionaryTerm(item, dimension, variables));
    if (terms.empty()) throw std::invalid_argument("a dictionary needs terms");
    return std::make_shared<const c::ExpansionDictionary>(dimension, std::move(terms));
}

}  // namespace dynsyslab::workflow::postprocess_detail
