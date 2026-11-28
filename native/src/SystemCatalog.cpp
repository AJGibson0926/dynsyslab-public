// The system catalogue (native/catalog/systems/<id>.json, embedded at build time by
// cmake/EmbedCatalog.cmake): each system's known Koopman eigenfunctions as data. One loader serves
// every system; adding an exact expression is a data change.
//
// A file (schema dynsyslab.system-catalog/v1) holds named references, shared "definitions" (ordered
// [name, expression] pairs) and "eigenfunctions". An entry has an id, title and kind (exact,
// principal or shell); expressions "psi", "psi_imag", "eigenvalue", "eigenvalue_imag" (for a shell
// entry the eigenvalue is its level map, constant on each level set); optional entry definitions;
// conditions and derivation as text; reference keys; and "requires", the checks that decide where
// it applies: lists of expressions that must be "zero", "positive" or "nonnegative", and "near":
// [expression, tolerance] pairs. A shell entry's generalized eigenfunctions are written out as
// "generalized_eigenfunction" (and "_imag") in the state and the band parameter p; by default
// delta(eigenvalue_imag - p) psi. Expressions name the system's state, its parameters, its
// dimensionless groups and the driver (driver_amplitude, driver_frequency_hz, driver_phase_radians).
#include "DynamicalKnowledgeInternal.hpp"

#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/Json.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core::detail {

const std::map<std::string, std::string>& systemCatalogTexts();  // generated: SystemCatalogData.cpp

std::map<std::string, double> catalogConstants(const SystemConfiguration& c) {
    auto out = c.parameters;
    for (const auto& [name, value] : c.dimensionlessGroups) out.emplace(name, value);
    out["driver_amplitude"] = c.driver.amplitude;
    out["driver_frequency_hz"] = c.driver.frequencyHz;
    out["driver_phase_radians"] = c.driver.phaseRadians;
    return out;
}

namespace {

std::string definitionsText(const json::Value& list) {
    std::string out;
    for (const auto& pair : list.array()) {
        const auto& item = pair.array();
        if (item.size() != 2) throw std::invalid_argument("a catalogue definition is [name, expression]");
        out += item[0].string() + " = " + item[1].string() + "; ";
    }
    return out;
}

bool satisfied(const json::Value& requires, const std::string& definitions, const std::vector<std::string>& names,
               const std::map<std::string, double>& constants, const std::string& where) {
    if (requires.isNull()) return true;
    std::vector<double> zero(names.size(), 0.0);
    const auto value = [&](const std::string& text) {
        return AnalyticExpression::parse(definitions + text, names, constants).evaluate(zero.data());
    };
    try {
        for (const auto& e : requires.at("zero").array()) if (value(e.string()) != 0) return false;
        for (const auto& e : requires.at("positive").array()) if (!(value(e.string()) > 0)) return false;
        for (const auto& e : requires.at("nonnegative").array()) if (!(value(e.string()) >= 0)) return false;
        for (const auto& e : requires.at("near").array()) {
            const auto& pair = e.array();
            if (pair.size() != 2) throw std::invalid_argument("a catalogue 'near' check is [expression, tolerance]");
            if (!(std::abs(value(pair[0].string())) <= pair[1].number())) return false;
        }
    } catch (const std::invalid_argument& error) {
        throw std::invalid_argument("catalogue entry " + where + ": " + error.what());
    }
    return true;
}

}  // namespace

std::vector<EigenfunctionEntry> catalogEigenfunctions(const SystemDescriptor& d, const SystemConfiguration& c) {
    const auto& texts = systemCatalogTexts();
    const auto found = texts.find(d.identifier);
    if (found == texts.end() || d.stateNames.empty()) return {};
    const auto catalog = json::parse(found->second);
    if (catalog.at("schema").string() != "dynsyslab.system-catalog/v1" || catalog.at("system").string() != d.identifier)
        throw std::invalid_argument("the catalogue file of " + d.identifier + " is not a dynsyslab.system-catalog/v1 entry for it");
    const auto constants = catalogConstants(c);
    const auto shared = definitionsText(catalog.at("definitions"));
    std::vector<EigenfunctionEntry> out;
    for (const auto& e : catalog.at("eigenfunctions").array()) {
        const auto id = e.at("id").string();
        const auto definitions = shared + (e.at("definitions").isNull() ? std::string() : definitionsText(e.at("definitions")));
        EigenfunctionEntry entry;
        entry.id = id;
        entry.title = e.at("title").string();
        entry.kind = e.at("kind").string("exact");
        const auto text = [&](const char* key) {  // a bare "0" stays bare: the loaders read it as zero
            const auto value = e.at(key).string("0");
            return value == "0" ? value : definitions + value;
        };
        entry.real = text("psi");
        entry.imag = text("psi_imag");
        entry.rate = text("eigenvalue");
        entry.rateImag = text("eigenvalue_imag");
        if (entry.kind == "shell") {  // the explicit family delta(level(x) - p) profile(x), unless stated
            const auto level = e.at("eigenvalue_imag").string("0");
            entry.familyReal = definitions + e.at("generalized_eigenfunction").string("dirac(" + level + " - p)*(" + e.at("psi").string("0") + ")");
            entry.familyImag = definitions + e.at("generalized_eigenfunction_imag").string("dirac(" + level + " - p)*(" + e.at("psi_imag").string("0") + ")");
        }
        entry.conditions = e.at("conditions").string();
        entry.derivation = e.at("derivation").string();
        for (const auto& key : e.at("references").array()) {
            const auto text = catalog.at("references").at(key.string()).string();
            entry.references.push_back(text.empty() ? key.string() : text);
        }
        entry.applicable = satisfied(e.at("requires"), shared, d.stateNames, constants, d.identifier + "/" + id);
        out.push_back(std::move(entry));
    }
    return out;
}

}  // namespace dynsyslab::core::detail
