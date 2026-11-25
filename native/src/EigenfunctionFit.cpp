#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <stdexcept>

namespace dynsyslab::core {
EigenfunctionFitResult polynomialEigenfunctionsFromLegacy(const KoopmanSINDyResult &f) {
    const bool external = f.coordinateProvider == "causal-history" || f.coordinateProvider == "measured";
    FunctionDomain q = external ? FunctionDomain::supplied(f.stateDimension)
        : !f.coordinateObservableDefinition.empty()
            ? FunctionDomain::fromStateMap(Observable::fromExpressionDefinition(json::parse(f.coordinateObservableDefinition)))
        : !f.coordinateObservableNames.empty() ? FunctionDomain::fromCatalog(f.coordinateObservableNames, f.stateDimension)
        : FunctionDomain::physicalState(f.stateDimension);
    if (!q.stateMap && f.dynamicalStateDimension) q.physicalStateDimension = f.dynamicalStateDimension;
    if (!f.coordinateObservableNames.empty()) q.names = f.coordinateObservableNames;
    if (f.coordinateProvider == "causal-history") q.providerIdentity = "causal-history";
    return {FunctionExpansion({f.stateDimension, f.exponents, f.center, f.scale}, f.coefficients, std::move(q)),
            Spectrum::fromStoredRates(f.continuousEigenvalues, f.sampleTime,
                                                   f.eigenvalueOrigin == "discrete-time"
                                                       ? EigenvalueDomain::discreteTime
                                                       : EigenvalueDomain::continuousTime,
                                                   f.logarithmBranches),
            fitDiagnosticsFromLegacy(f)};
}
EigenfunctionFitDiagnostics fitDiagnosticsFromLegacy(const KoopmanSINDyResult &f) {
    // The old DTO was also used for EDMD/manual coefficients. Its C++ type does
    // not establish which algorithm produced them or what its residual meant.
    return {{}, {}, f.samples, f.directFits, f.residualRms, f.residualMaximum};
}
KoopmanSINDyResult legacyPolynomialEigenfunctions(const FunctionExpansion &p, const Spectrum &spectrum,
                                                  const EigenfunctionFitDiagnostics &diagnostics) {
    KoopmanSINDyResult f;
    const auto &q = p.coordinates();
    f.continuousEigenvalues = spectrum.continuousRates(); // explicit failure for zero mu
    f.stateDimension = p.coordinateDimension();
    f.featureCount = p.dictionary().size();
    if (p.monomial()) {  // other bases have no exponent table in this legacy snapshot
        const auto &d = p.monomialTable();
        f.exponents = d.exponents;
        f.center = d.center;
        f.scale = d.scale;
    }
    f.coefficients = p.coefficients();
    f.dynamicalStateDimension = q.stateMap ? q.stateMap->stateDimension() : q.physicalStateDimension.value_or(0);
    f.coordinateObservableNames = q.names;
    f.coordinateObservableDefinition = q.stateMapDefinition();
    f.coordinateProvider = q.requiresExternalValues()
        ? (q.computableFromMeasurements() ? "causal-history" : "measured") : "instantaneous-observable";
    f.sampleTime = spectrum.sampleInterval();
    f.logarithmBranches = spectrum.logarithmBranches();
    f.eigenvalueOrigin = spectrum.domain() == EigenvalueDomain::discreteTime
                             ? "discrete-time"
                             : "continuous-time";
    f.samples = diagnostics.sampleCount;
    f.directFits = diagnostics.directFitCount;
    f.residualRms = diagnostics.residualRms;
    f.residualMaximum = diagnostics.residualMaximum;
    return f;
}
EigenfunctionFitDiagnostics
EigenfunctionFitDiagnostics::extract(const std::vector<std::size_t> &modes) const {
    auto result = *this;
    result.residualRms.clear();
    result.residualMaximum.clear();
    for (auto j : modes) {
        if (!residualRms.empty())
            result.residualRms.push_back(residualRms.at(j));
        if (!residualMaximum.empty())
            result.residualMaximum.push_back(residualMaximum.at(j));
    }
    // A subset is not another fitting run; retain the parent fit count.
    return result;
}
} // namespace dynsyslab::core
