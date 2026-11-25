// Loads a system's hand-derived eigenfunctions (DynamicalKnowledge::eigenfunctions) as Eigenfunction
// objects (or as their shared expansion plus spectrum) and checks them against the equations.
#include <dynsyslab/core/KnownEigenfunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

EigenfunctionFitResult libraryEigenfunctionBank(const ODESystem& system, const std::vector<std::string>& ids,
                                                const std::vector<double>& reference) {
    const auto knowledge = system.dynamicalKnowledge();
    const auto& names = system.descriptor().stateNames;
    const auto& parameters = system.configuration().parameters;
    std::vector<const EigenfunctionEntry*> chosen;
    if (ids.empty()) {
        for (const auto& e : knowledge.eigenfunctions) if (e.applicable) chosen.push_back(&e);
    } else {
        for (const auto& id : ids) {
            const auto found = std::find_if(knowledge.eigenfunctions.begin(), knowledge.eigenfunctions.end(),
                                            [&](const EigenfunctionEntry& e) { return e.id == id; });
            if (found == knowledge.eigenfunctions.end())
                throw std::invalid_argument(system.descriptor().identifier + " has no library eigenfunction " + id);
            if (!found->applicable)
                throw std::invalid_argument("library eigenfunction " + id + " does not apply here: " + found->conditions);
            chosen.push_back(&*found);
        }
    }
    if (chosen.empty()) throw std::invalid_argument(system.descriptor().identifier + " has no applicable library eigenfunctions");
    auto domain = FunctionDomain::physicalState(names.size());
    domain.names = names;
    std::vector<AnalyticExpression> parts;
    std::vector<std::complex<double>> rates;
    std::vector<std::string> outputs;
    if (!reference.empty() && reference.size() != names.size())
        throw std::invalid_argument("library reference state must be a physical state");
    const std::vector<double> at = reference.empty() ? std::vector<double>(names.size(), 0.0) : reference;
    for (const auto* e : chosen) {
        if (e->kind == "shell" && reference.empty())
            throw std::invalid_argument("library eigenfunction " + e->id + " has an orbit-dependent eigenvalue; pass a reference state");
        parts.push_back(AnalyticExpression::parse(e->real, names, parameters));
        parts.push_back(AnalyticExpression::parse(e->imag, names, parameters));
        rates.emplace_back(AnalyticExpression::parse(e->rate, names, parameters).evaluate(at.data()),
                           AnalyticExpression::parse(e->rateImag, names, parameters).evaluate(at.data()));
        if (!std::isfinite(rates.back().real()) || !std::isfinite(rates.back().imag()))
            throw std::invalid_argument("library eigenvalue of " + e->id + " is not finite here (outside its region?)");
        outputs.push_back(e->id);
    }
    // Real and imaginary parts share one dictionary; psi_j = row(2j) + i row(2j+1).
    const auto split = FunctionExpansion::fromExpressions(parts, domain);
    const std::size_t K = split.dictionary().size(), m = chosen.size();
    std::vector<std::complex<double>> c(m * K);
    for (std::size_t j = 0; j < m; ++j)
        for (std::size_t k = 0; k < K; ++k)
            c[j * K + k] = split.coefficients()[2 * j * K + k] + std::complex<double>(0, 1) * split.coefficients()[(2 * j + 1) * K + k];
    EigenfunctionFitDiagnostics provenance;
    provenance.algorithm = "system-library";
    return {FunctionExpansion(split.dictionaryHandle(), std::move(c), domain, std::move(outputs)),
            Spectrum::fromContinuousRates(std::move(rates)), std::move(provenance)};
}

std::vector<Eigenfunction> libraryEigenfunctions(const ODESystem& system, const std::vector<std::string>& ids,
                                                 const std::vector<double>& reference) {
    const auto bank = libraryEigenfunctionBank(system, ids, reference);
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
