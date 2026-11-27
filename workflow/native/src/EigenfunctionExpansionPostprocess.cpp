// Postprocess method eigenfunction_expansion: eigenfunctions delivered by hand. Like any other
// method it supplies a dictionary and coefficients; it just takes them from the job.
//
//   {"name": "eigenfunction_expansion", "settings": {
//      "domain": "state" | ["echo", "echo_dot"],        state, or catalog observables of the state
//      "variables": ["r", "v"],                          names used in expressions (default q0, q1, ...)
//      "eigenvalues": [[-0.1, 1.7]],                     continuous rates, [re, im]
//      "dictionary": [{"monomial": [1, 0]}, {"exponential": [0, -0.5]}, {"expression": "sin(r)*v"},
//                     {"gaussian": {"node": [1, 0], "width": 0.3}}, {"cosine": [2, 0], "phase": 0.1}],
//      "coefficients": [[[1, 0], [0.5, -0.2], 0.1, [0, 1], 2]]}}   one row per eigenvalue
//
// Shorthand: "eigenfunctions": [{"real": "r", "imag": "-v/1.7"}] decomposes each text into
// dictionary terms (FunctionExpansion::fromExpressions): polynomial parts become monomials, other
// pieces expression terms, so hand-written and fitted eigenfunctions share terms.
//
// Library: "library": true | ["x_1", "eta"] loads the source system's own hand-derived
// eigenfunctions (DynamicalKnowledge::eigenfunctions) with its parameters; no other settings.
#include "TimeSeriesSource.hpp"
#include "DictionaryTerms.hpp"
#include "ExpansionArchive.hpp"
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

std::complex<double> coefficient(const J& value) {
    if (std::holds_alternative<double>(value.storage())) return value.number();
    const auto& pair = value.array();
    if (pair.size() != 2) throw std::invalid_argument("a complex coefficient is [real, imaginary]");
    return {pair[0].number(), pair[1].number()};
}


}  // namespace

std::filesystem::path runEigenfunctionExpansion(const std::filesystem::path& sourcePath, const J& settings,
                                                const std::filesystem::path& root, const std::string& label) {
    const auto source = loadSource(sourcePath);
    const auto configuration = systemConfiguration(source);
    const auto system = c::createSystem(configuration);
    const auto publishBank = [&](const c::EigenfunctionFitResult& bank, const char* sourceKind) {
        (void)c::KoopmanModel(bank, configuration);  // the eigenfunctions must form a model of this system
        auto meta = json::parse(metadata(source, "eigenfunction-expansion", label,
            std::string("\"eigenfunction_source\":\"") + sourceKind + "\",\"extraction_algorithm\":\"eigenfunction_expansion\""));
        meta.object()["eigenfunction_expansion_settings"] = settings;
        c::NpzRecord archive;
        std::string reason;
        const c::DataOrigin origin{"eigenfunction_expansion", settings, {dataReference(source.path, "time_series")}};
        if (!writeExpansionKrom(archive, bank, origin, {}, {}, {}, json::dump(meta), reason))
            throw std::invalid_argument("eigenfunction_expansion: " + reason);
        return publishObject(root, source, dynsyslab::core::object_kind::krom, "koopman-eigenfunction-expansion", std::move(archive));
    };
    if (const auto& library = settings.at("library"); !library.isNull()) {
        std::vector<std::string> ids;
        if (!std::holds_alternative<bool>(library.storage()))
            for (const auto& id : library.array()) ids.push_back(id.string());
        std::vector<double> reference;  // the orbit whose eigenvalues shell entries take
        if (const auto& at = settings.at("reference_state"); !at.isNull())
            for (const auto& value : at.array()) reference.push_back(value.number());
        return publishBank(c::libraryEigenfunctionBank(*system, ids, reference), "system-library");
    }
    c::FunctionDomain domain;
    const auto& requested = settings.at("domain");
    if (requested.isNull() || (std::holds_alternative<std::string>(requested.storage()) && requested.string() == "state")) {
        domain = c::FunctionDomain::physicalState(system->descriptor().stateDimension);
    } else {
        std::vector<std::string> names;
        for (const auto& name : requested.array()) names.push_back(name.string());
        domain = c::FunctionDomain::fromCatalog(names, names.size(), nullptr);
    }
    std::vector<std::string> variables;
    for (const auto& name : settings.at("variables").array()) variables.push_back(name.string());
    if (variables.empty())
        for (std::size_t k = 0; k < domain.dimension; ++k) variables.push_back("q" + std::to_string(k));
    if (variables.size() != domain.dimension)
        throw std::invalid_argument("eigenfunction_expansion needs one variable name per domain coordinate");

    std::vector<std::complex<double>> rates;
    for (const auto& value : settings.at("eigenvalues").array()) rates.push_back(coefficient(value));
    if (rates.empty()) throw std::invalid_argument("eigenfunction_expansion needs eigenvalues");

    std::vector<c::DictionaryTerm> terms;
    std::vector<std::complex<double>> coefficients;
    if (!settings.at("dictionary").isNull()) {
        for (const auto& spec : settings.at("dictionary").array()) terms.push_back(dictionaryTerm(spec, domain.dimension, variables));
        const auto& rows = settings.at("coefficients").array();
        if (rows.size() != rates.size()) throw std::invalid_argument("eigenfunction_expansion needs one coefficient row per eigenvalue");
        for (const auto& row : rows) {
            if (row.array().size() != terms.size()) throw std::invalid_argument("each coefficient row needs one entry per dictionary term");
            for (const auto& value : row.array()) coefficients.push_back(coefficient(value));
        }
    } else {  // shorthand: each written function decomposes into dictionary terms
        const auto& written = settings.at("eigenfunctions").array();
        if (written.size() != rates.size()) throw std::invalid_argument("eigenfunction_expansion needs one eigenfunction per eigenvalue");
        std::vector<c::AnalyticExpression> parts;
        for (const auto& item : written)
            for (const char* part : {"real", "imag"}) parts.push_back(c::AnalyticExpression::parse(item.at(part).string("0"), variables));
        const auto split = c::FunctionExpansion::fromExpressions(parts, domain);
        terms = split.dictionary().terms();
        const std::size_t n = terms.size();
        coefficients.assign(rates.size() * n, {});
        for (std::size_t j = 0; j < rates.size(); ++j)
            for (std::size_t k = 0; k < n; ++k)
                coefficients[j * n + k] = split.coefficients()[2 * j * n + k]
                                        + std::complex<double>(0.0, 1.0) * split.coefficients()[(2 * j + 1) * n + k];
    }
    return publishBank({c::FunctionExpansion(std::make_shared<const c::ExpansionDictionary>(domain.dimension, std::move(terms)),
                                             std::move(coefficients), domain),
                        c::Spectrum::fromContinuousRates(rates), {}},
                       "hand-written");
}

}  // namespace dynsyslab::workflow::postprocess_detail
