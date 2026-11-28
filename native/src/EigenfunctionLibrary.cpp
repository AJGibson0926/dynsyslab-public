// Loads a system's hand-derived eigenfunctions (DynamicalKnowledge::eigenfunctions) as Eigenfunction
// objects (or as their shared expansion plus spectrum) and checks them against the equations.
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include "DynamicalKnowledgeInternal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

namespace {

// The applicable entries named by `ids` (all of `kind` "shell" or all others, per `shells`), in that order;
// by default every applicable entry of that sort.
std::vector<const EigenfunctionEntry*> chosenEntries(const ODESystem& system, const DynamicalKnowledge& knowledge,
                                                     const std::vector<std::string>& ids, bool shells) {
    std::vector<const EigenfunctionEntry*> chosen;
    for (const auto& id : ids) {
        const auto found = std::find_if(knowledge.eigenfunctions.begin(), knowledge.eigenfunctions.end(),
                                        [&](const EigenfunctionEntry& e) { return e.id == id; });
        if (found == knowledge.eigenfunctions.end())
            throw std::invalid_argument(system.descriptor().identifier + " has no library eigenfunction " + id);
        if (!found->applicable)
            throw std::invalid_argument("library eigenfunction " + id + " does not apply here: " + found->conditions);
        if ((found->kind == "shell") == shells) chosen.push_back(&*found);
    }
    if (ids.empty())
        for (const auto& e : knowledge.eigenfunctions)
            if (e.applicable && (e.kind == "shell") == shells) chosen.push_back(&e);
    return chosen;
}

FunctionDomain stateDomain(const ODESystem& system) {
    auto domain = FunctionDomain::physicalState(system.descriptor().stateNames.size());
    domain.names = system.descriptor().stateNames;
    return domain;
}

// One complex function per pair of real expressions (real, imaginary), on one shared dictionary.
FunctionExpansion complexExpansion(const std::vector<AnalyticExpression>& parts, const FunctionDomain& domain,
                                   std::vector<std::string> outputs) {
    const auto split = FunctionExpansion::fromExpressions(parts, domain);
    const std::size_t K = split.dictionary().size(), m = parts.size() / 2;
    std::vector<std::complex<double>> c(m * K);
    for (std::size_t j = 0; j < m; ++j)
        for (std::size_t k = 0; k < K; ++k)
            c[j * K + k] = split.coefficients()[2 * j * K + k] + std::complex<double>(0, 1) * split.coefficients()[(2 * j + 1) * K + k];
    return FunctionExpansion(split.dictionaryHandle(), std::move(c), domain, std::move(outputs));
}

}  // namespace

std::vector<SpectralBand> libraryBands(const ODESystem& system, const std::vector<std::string>& ids,
                                       const std::vector<double>& grid) {
    const auto knowledge = system.dynamicalKnowledge();
    const auto& names = system.descriptor().stateNames;
    const auto parameters = detail::catalogConstants(system.configuration());  // parameters, groups, driver
    const auto domain = stateDomain(system);
    std::vector<SpectralBand> out;
    for (const auto* e : chosenEntries(system, knowledge, ids, true)) {
        // psi_p = delta(level(x) - p) e^{i n phi(x)} with level(x) = n Omega(H(x)), the frequency on the shell
        // through x: a conservative entry's eigenvalue is i level(x) (Mezic 2020, section 11).
        if (e->rate != "0")
            throw std::invalid_argument("library shell entry " + e->id + " has a nonzero growth rate; its band is not on the imaginary axis");
        SpectralBand band;
        band.name = e->id;
        band.support = band_support::imaginaryAxis;
        band.method = "system-library";
        band.assumptions = e->conditions;
        for (const auto& reference : e->references) band.assumptions += (band.assumptions.empty() ? "" : "; ") + reference;
        band.grid = grid;
        band.eigenfunctions.functions = std::make_shared<const FunctionExpansion>(complexExpansion(
            {AnalyticExpression::parse(e->real, names, parameters), AnalyticExpression::parse(e->imag, names, parameters)},
            domain, {e->id}));
        if (!e->familyReal.empty()) {  // psi_p(x) written out on (state, p)
            FunctionDomain both;
            both.dimension = names.size() + 1;
            both.names = names;
            both.names.push_back("p");
            band.eigenfunctions.family = std::make_shared<const FunctionExpansion>(complexExpansion(
                {AnalyticExpression::parse(e->familyReal, both.names, parameters), AnalyticExpression::parse(e->familyImag, both.names, parameters)},
                both, {e->id + "_family"}));
        }
        band.eigenfunctions.level = std::make_shared<const FunctionExpansion>(FunctionExpansion::fromExpressions(
            {AnalyticExpression::parse(e->rateImag, names, parameters)}, domain, {e->id + "_frequency"}));
        out.push_back(std::move(band));
    }
    return out;
}

EigenfunctionFitResult libraryEigenfunctionBank(const ODESystem& system, const std::vector<std::string>& ids) {
    const auto knowledge = system.dynamicalKnowledge();
    const auto& names = system.descriptor().stateNames;
    const auto parameters = detail::catalogConstants(system.configuration());  // parameters, groups, driver
    const auto chosen = chosenEntries(system, knowledge, ids, false);
    if (chosen.empty())
        throw std::invalid_argument(system.descriptor().identifier + " has no applicable library eigenfunctions with constant eigenvalues"
                                    " (shell entries are generalized eigenfunctions on spectral bands)");
    const auto domain = stateDomain(system);
    std::vector<AnalyticExpression> parts;
    std::vector<std::complex<double>> rates;
    std::vector<std::string> outputs;
    const std::vector<double> at(names.size(), 0.0);
    for (const auto* e : chosen) {
        parts.push_back(AnalyticExpression::parse(e->real, names, parameters));
        parts.push_back(AnalyticExpression::parse(e->imag, names, parameters));
        rates.emplace_back(AnalyticExpression::parse(e->rate, names, parameters).evaluate(at.data()),
                           AnalyticExpression::parse(e->rateImag, names, parameters).evaluate(at.data()));
        if (!std::isfinite(rates.back().real()) || !std::isfinite(rates.back().imag()))
            throw std::invalid_argument("library eigenvalue of " + e->id + " is not finite here (outside its region?)");
        outputs.push_back(e->id);
    }
    EigenfunctionFitDiagnostics provenance;
    provenance.algorithm = "system-library";
    // The atomic part, and the shell entries requested (or all, by default) as spectral bands.
    auto spectrum = Spectrum::fromContinuousRates(std::move(rates));
    for (auto& band : libraryBands(system, ids)) spectrum.addBand(std::move(band));
    return {complexExpansion(parts, domain, std::move(outputs)), std::move(spectrum), std::move(provenance)};
}

std::vector<Eigenfunction> libraryEigenfunctions(const ODESystem& system, const std::vector<std::string>& ids) {
    for (const auto& id : ids)
        for (const auto& e : system.dynamicalKnowledge().eigenfunctions)
            if (e.id == id && e.kind == "shell")
                throw std::invalid_argument("library entry " + id + " is a generalized eigenfunction on a spectral band; load it with libraryBands");
    const auto bank = libraryEigenfunctionBank(system, ids);
    const auto shared = std::make_shared<const FunctionExpansion>(bank.functions);
    std::vector<Eigenfunction> out;
    for (std::size_t j = 0; j < shared->outputCount(); ++j) {
        const auto& id = shared->outputNames()[j];
        out.push_back(Eigenfunction::fromRepresentation({id, id}, shared, j, bank.spectrum.eigenvalue(j)));
    }
    return out;
}

std::vector<double> eigenfunctionGeneratorResidual(const ODESystem& system, const std::vector<Eigenfunction>& psi,
                                                   const std::vector<double>& states) {
    const std::size_t n = system.descriptor().stateDimension, m = psi.size();
    for (const auto& e : psi)
        if (e.stateDimension() != n) throw std::invalid_argument("generator residual needs eigenfunctions of the physical state");
    if (states.empty() || states.size() % n) throw std::invalid_argument("residual states must be rows of the physical state");
    std::vector<double> worst(m, 0.0), scale(m, 0.0), f(n);
    for (std::size_t s = 0; s < states.size() / n; ++s) {
        const double* x = states.data() + s * n;
        system.rhs(0.0, x, f.data());
        for (std::size_t j = 0; j < m; ++j) {
            const auto J = psi[j].jacobian(x);
            const auto v = psi[j].evaluate(x).front();
            std::complex<double> lie{};
            for (std::size_t i = 0; i < n; ++i) lie += J[i] * f[i];
            worst[j] = std::max(worst[j], std::abs(lie - psi[j].eigenvalue().continuous * v));
            scale[j] = std::max(scale[j], std::abs(v));
        }
    }
    for (std::size_t j = 0; j < m; ++j) worst[j] /= std::max(scale[j], 1e-300);
    return worst;
}

}  // namespace dynsyslab::core
