#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <cmath>

namespace dynsyslab::core {
KnownEigenfunctionResult knownPolynomialEigenfunctions(const ODESystem& system) {
    KnownEigenfunctionResult result;
    if (system.descriptor().identifier != "quadratic_triangular") {
        result.reason = "No exact polynomial eigenfunction family registered for this system";
        return result;
    }
    result.identifier = "quadratic-triangular-gradient-closed/v1";
    result.references = {"https://doi.org/10.1371/journal.pone.0150171",
                         "https://arxiv.org/abs/1707.01146"};
    result.conditions = "Eigenfunctions of the zero-input drift on R^2; b != 2*a or c == 0";
    const auto& p = system.configuration().parameters;
    const double a = p.at("a"), b = p.at("b"), c = p.at("c");
    if (b == 2*a && c != 0) {
        result.reason = "Quadratic resonance b=2*a, c!=0: the state-containing polynomial lift "
                        "has a Jordan block; this diagonal eigenfunction family is unavailable";
        return result;
    }
    const double k = c == 0 ? 0 : c/(b-2*a);
    if (!std::isfinite(k) || !std::isfinite(2*k) || !std::isfinite(2*a)) {
        result.reason = "Eigenfunction coefficients or rates exceed finite numeric range";
        return result;
    }
    EigenfunctionCoordinates coordinates;
    coordinates.kind = EigenfunctionCoordinateKind::physicalState;
    coordinates.dimension = 2;
    coordinates.physicalStateDimension = 2;
    coordinates.names = {"x_1", "x_2"};
    // Brunton et al. (2016), eqs. 25-27; Kaiser et al. (2021), eqs. 50-52.
    // Native family x1'=a*x1, x2'=b*x2+c*x1^2; eta=x2+k*x1^2.
    result.eigenfunctions.emplace(
        PolynomialDictionary{2, {0,0, 1,0, 2,0, 0,1}, {}, {}},
        std::vector<std::complex<double>>{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,k,1},
        EigenfunctionSpectrum::fromContinuousRates({0.,a,2*a,b}), coordinates);
    result.names = {"one", "x_1", "x_1_squared", "eta"};
    result.derivativeOperators = {
        0,0,0,0, 1,0,0,0, 0,2,0,0, 0,2*k,0,0,
        0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0};
    result.stateReadout = {0,1,0,0, 0,0,-k,1};
    result.status = "exact";
    result.reason = "Algebraic identities for the configured zero-input drift; no regression";
    return result;
}

json::Value knownEigenfunctionsJson(const KnownEigenfunctionResult& r) {
    const auto strings = [](const auto& values) {
        json::Value::Array out;
        for (const auto& v : values) out.emplace_back(v);
        return out;
    };
    const auto numbers = [](const auto& values) {
        json::Value::Array out;
        for (const auto& v : values) out.emplace_back(static_cast<double>(v));
        return out;
    };
    const auto complex = [](const auto& values) {
        json::Value::Array out;
        for (const auto& v : values)
            out.emplace_back(json::Value::Array{v.real(),v.imag()});
        return out;
    };
    json::Value::Object out{{"schema","dynsyslab.known-eigenfunctions/v1"},
        {"id",r.identifier},{"status",r.status},{"reason",r.reason},
        {"conditions",r.conditions},{"references",strings(r.references)}};
    if (!r.eigenfunctions) return out;
    const auto& bank = *r.eigenfunctions;
    const auto& dictionary = bank.dictionary();
    out["names"] = strings(r.names);
    out["coordinates"] = eigenfunctionCoordinatesJson(bank.coordinates());
    out["coordinate_dimension"] = double(bank.coordinateDimension());
    out["mode_count"] = double(bank.eigenfunctionCount());
    out["feature_count"] = double(dictionary.featureCount());
    out["exponents"] = numbers(dictionary.exponents);
    out["center"] = numbers(dictionary.center);
    out["scale"] = numbers(dictionary.scale);
    out["coefficients"] = complex(bank.coefficients());
    out["continuous_rates"] = complex(bank.spectrum().continuousRates());
    out["derivative_operators"] = numbers(r.derivativeOperators);
    out["state_readout"] = numbers(r.stateReadout);
    out["array_layout"] = "row-major; coefficients mode x feature; complex values [real,imag]; "
                          "derivatives coordinate x mode x mode; readout state x mode";
    return out;
}
} // namespace dynsyslab::core
