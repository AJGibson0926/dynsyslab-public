#include <dynsyslab/core/FunctionExpansion.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace dynsyslab::core {

struct FunctionExpansion::Data {
    std::shared_ptr<const ExpansionDictionary> dictionary;
    std::vector<std::complex<double>> coefficients;
    FunctionDomain coordinates;
    std::vector<std::string> outputNames;
};

FunctionExpansion::FunctionExpansion(std::shared_ptr<const ExpansionDictionary> dictionary,
                                         std::vector<std::complex<double>> c, FunctionDomain q,
                                         std::vector<std::string> names) {
    if (!dictionary) throw std::invalid_argument("Expansion dictionary is empty");
    q.validate();
    const auto features = dictionary->size();
    if (q.dimension != dictionary->dimension() || c.empty() || c.size() % features)
        throw std::invalid_argument("Expansion coefficients must be outputs x dictionary terms on the dictionary's domain");
    if (!std::all_of(c.begin(), c.end(), [](auto z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }))
        throw std::invalid_argument("Expansion coefficients must be finite");
    if (!names.empty() && names.size() != c.size() / features)
        throw std::invalid_argument("Expansion output names must match the outputs");
    data_ = std::make_shared<const Data>(Data{std::move(dictionary), std::move(c), std::move(q), std::move(names)});
}

FunctionExpansion::FunctionExpansion(PolynomialDictionary monomials, std::vector<std::complex<double>> c, FunctionDomain q)
    : FunctionExpansion(std::make_shared<const ExpansionDictionary>(ExpansionDictionary::monomials(std::move(monomials))),
                        std::move(c), std::move(q)) {}

FunctionExpansion FunctionExpansion::withCoordinates(FunctionDomain q) const {
    return FunctionExpansion(data_->dictionary, data_->coefficients, std::move(q), data_->outputNames);
}

const ExpansionDictionary& FunctionExpansion::dictionary() const noexcept { return *data_->dictionary; }
std::shared_ptr<const ExpansionDictionary> FunctionExpansion::dictionaryHandle() const noexcept { return data_->dictionary; }
bool FunctionExpansion::monomial() const noexcept { return data_->dictionary->monomialTable() != nullptr; }
const PolynomialDictionary& FunctionExpansion::monomialTable() const {
    const auto* table = data_->dictionary->monomialTable();
    if (!table) throw std::logic_error("this expansion's dictionary is not a single monomial table");
    return *table;
}
const std::vector<std::complex<double>>& FunctionExpansion::coefficients() const noexcept { return data_->coefficients; }
const FunctionDomain& FunctionExpansion::coordinates() const noexcept { return data_->coordinates; }
const std::vector<std::string>& FunctionExpansion::outputNames() const noexcept { return data_->outputNames; }
std::size_t FunctionExpansion::outputCount() const noexcept { return data_->coefficients.size() / data_->dictionary->size(); }

std::complex<double> FunctionExpansion::evaluate(std::size_t output, const double* q) const {
    if (!q || output >= outputCount()) throw std::invalid_argument("Invalid expansion evaluation request");
    return data_->dictionary->combine(coefficients().data() + output * data_->dictionary->size(), q);
}

std::vector<std::complex<double>> FunctionExpansion::evaluate(const double* q) const {
    std::vector<std::complex<double>> values;
    values.reserve(outputCount());
    if (monomial()) {  // the established polynomial arithmetic, output by output
        for (std::size_t r = 0; r < outputCount(); ++r) values.push_back(evaluate(r, q));
        return values;
    }
    const auto phi = data_->dictionary->evaluate(q);  // one dictionary evaluation shared by every output
    const std::size_t n = phi.size();
    for (std::size_t r = 0; r < outputCount(); ++r) {
        std::complex<double> sum{};
        for (std::size_t k = 0; k < n; ++k) sum += coefficients()[r * n + k] * phi[k];
        values.push_back(sum);
    }
    return values;
}

std::vector<std::complex<double>> FunctionExpansion::gradient(std::size_t output, const double* q) const {
    if (!q || output >= outputCount()) throw std::invalid_argument("Invalid expansion gradient request");
    return data_->dictionary->combineGradient(coefficients().data() + output * data_->dictionary->size(), q);
}

std::vector<std::complex<double>> FunctionExpansion::jacobian(const double* q) const {
    std::vector<std::complex<double>> out;
    out.reserve(outputCount() * coordinateDimension());
    for (std::size_t r = 0; r < outputCount(); ++r) {
        const auto g = gradient(r, q);
        out.insert(out.end(), g.begin(), g.end());
    }
    return out;
}

namespace {
// Terms are identified by value: their serialized description.
std::vector<std::string> termKeys(const ExpansionDictionary& d) {
    const auto description = d.describe();
    std::vector<std::string> keys;
    for (const auto& t : description.at("terms").array()) keys.push_back(json::dump(t));
    return keys;
}
}  // namespace

std::shared_ptr<const ExpansionDictionary> mergeDictionaries(const ExpansionDictionary& first, const ExpansionDictionary& second) {
    if (first.dimension() != second.dimension()) throw std::invalid_argument("Dictionaries on different domains cannot merge");
    const auto da = first.describe(), db = second.describe();
    const auto a = da.at("terms").array(), b = db.at("terms").array();
    std::map<std::string, bool> seen;
    for (const auto& t : a) seen[json::dump(t)] = true;
    json::Value::Array merged = a;
    for (const auto& t : b)
        if (!seen.count(json::dump(t))) { merged.push_back(t); seen[json::dump(t)] = true; }
    return std::make_shared<const ExpansionDictionary>(ExpansionDictionary::fromJson(json::Value::Object{
        {"schema", "dynsyslab.expansion-dictionary/v1"}, {"dimension", double(first.dimension())}, {"terms", merged}}));
}

FunctionExpansion FunctionExpansion::inDictionary(std::shared_ptr<const ExpansionDictionary> larger) const {
    if (!larger || larger->dimension() != coordinateDimension()) throw std::invalid_argument("Target dictionary has another domain");
    std::map<std::string, std::size_t> index;
    const auto largerKeys = termKeys(*larger), ownKeys = termKeys(*data_->dictionary);
    for (std::size_t k = 0; k < larger->size(); ++k) index[largerKeys[k]] = k;
    const std::size_t n = data_->dictionary->size(), m = larger->size(), outputs = outputCount();
    std::vector<std::complex<double>> c(outputs * m);
    for (std::size_t k = 0; k < n; ++k) {
        const auto found = index.find(ownKeys[k]);
        if (found == index.end()) throw std::invalid_argument("Target dictionary lacks a term of this expansion");
        for (std::size_t r = 0; r < outputs; ++r) c[r * m + found->second] += coefficients()[r * n + k];
    }
    return FunctionExpansion(std::move(larger), std::move(c), coordinates(), outputNames());
}

std::shared_ptr<const FunctionRepresentation> FunctionExpansion::selectOutputs(const std::vector<std::size_t>& outputs) const {
    return std::make_shared<const FunctionExpansion>(select(outputs));
}

FunctionExpansion FunctionExpansion::select(const std::vector<std::size_t>& outputs) const {
    const auto n = data_->dictionary->size();
    std::vector<std::complex<double>> c;
    std::vector<std::string> names;
    for (const auto r : outputs) {
        if (r >= outputCount()) throw std::out_of_range("Expansion output index");
        c.insert(c.end(), coefficients().begin() + std::ptrdiff_t(r * n), coefficients().begin() + std::ptrdiff_t((r + 1) * n));
        if (!outputNames().empty()) names.push_back(outputNames()[r]);
    }
    return FunctionExpansion(data_->dictionary, std::move(c), coordinates(), std::move(names));
}

}  // namespace dynsyslab::core
