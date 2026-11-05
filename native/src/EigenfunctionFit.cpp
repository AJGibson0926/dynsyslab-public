#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <stdexcept>

namespace dynsyslab::core {
PolynomialEigenfunctionSet polynomialEigenfunctionsFromLegacy(const KoopmanSINDyResult &f) {
    EigenfunctionCoordinates q;
    q.dimension = f.stateDimension;
    if (f.dynamicalStateDimension)
        q.physicalStateDimension = f.dynamicalStateDimension;
    q.names = f.coordinateObservableNames;
    q.expressionDefinition = f.coordinateObservableDefinition;
    q.providerIdentity = f.coordinateProvider;
    q.kind = f.coordinateProvider == "causal-history" ? EigenfunctionCoordinateKind::causalHistory
             : f.coordinateProvider == "measured"     ? EigenfunctionCoordinateKind::measured
             : (!q.names.empty() || !q.expressionDefinition.empty())
                 ? EigenfunctionCoordinateKind::analyticObservable
                 : EigenfunctionCoordinateKind::physicalState;
    if (q.kind == EigenfunctionCoordinateKind::physicalState && !q.physicalStateDimension)
        q.physicalStateDimension = q.dimension;
    return {{f.stateDimension, f.exponents, f.center, f.scale},
            f.coefficients,
            EigenfunctionSpectrum::fromStoredRates(f.continuousEigenvalues, f.sampleTime,
                                                   f.eigenvalueOrigin == "discrete-time"
                                                       ? EigenvalueDomain::discreteTime
                                                       : EigenvalueDomain::continuousTime,
                                                   f.logarithmBranches),
            std::move(q)};
}
EigenfunctionFitDiagnostics fitDiagnosticsFromLegacy(const KoopmanSINDyResult &f) {
    // The old DTO was also used for EDMD/manual coefficients. Its C++ type does
    // not establish which algorithm produced them or what its residual meant.
    return {{}, {}, f.samples, f.directFits, f.residualRms, f.residualMaximum};
}
KoopmanSINDyResult legacyPolynomialEigenfunctions(const PolynomialEigenfunctionSet &p,
                                                  const EigenfunctionFitDiagnostics &diagnostics) {
    KoopmanSINDyResult f;
    const auto &d = p.dictionary();
    const auto &q = p.coordinates();
    f.continuousEigenvalues = p.spectrum().continuousRates(); // explicit failure for zero mu
    f.stateDimension = d.coordinateDimension;
    f.featureCount = d.featureCount();
    f.exponents = d.exponents;
    f.center = d.center;
    f.scale = d.scale;
    f.coefficients = p.coefficients();
    f.dynamicalStateDimension = q.physicalStateDimension.value_or(0);
    f.coordinateObservableNames = q.names;
    f.coordinateObservableDefinition = q.expressionDefinition;
    f.coordinateProvider = q.kind == EigenfunctionCoordinateKind::causalHistory ? "causal-history"
                           : q.kind == EigenfunctionCoordinateKind::measured
                               ? "measured"
                               : "instantaneous-observable";
    f.sampleTime = p.spectrum().sampleInterval();
    f.logarithmBranches = p.spectrum().logarithmBranches();
    f.eigenvalueOrigin = p.spectrum().domain() == EigenvalueDomain::discreteTime
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
