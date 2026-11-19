#include "KlqeProcessCovariance.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dynsyslab::workflow {
namespace {
namespace core = dynsyslab::core;
namespace json = core::json;

std::vector<double> numbers(const json::Value& v) {
    std::vector<double> out;
    for (const auto& x : v.array()) out.push_back(x.number());
    return out;
}

std::vector<double> plain(const json::Value& v, std::size_t n) {  // scalar, diagonal or square
    auto values = numbers(v);
    if (values.size() == n * n) return values;
    std::vector<double> m(n * n, 0.0);
    if (!values.empty() && values.size() != n)
        throw std::invalid_argument("process covariance must be scalar, diagonal, square or a records object");
    for (std::size_t k = 0; k < n; ++k) m[k * n + k] = values.empty() ? v.number(1) : values[k];
    return m;
}

std::vector<double> constantMap(const core::KoopmanModel& model, const std::string& representation) {
    if (representation.empty() || representation == "none")
        throw std::invalid_argument("residual or input-direction process covariance needs the estimator's control_representation");
    const std::vector<double> origin(model.stateDimension(), 1.0);
    return model.controlMapSplit(origin.data(), representation);  // constant maps ignore the state
}

// W = E[r r^T] / h from the K-ROM's one-step prediction residual on recorded trajectories.
std::vector<double> residualCovariance(const json::Value& spec, const core::KoopmanModel& model,
                                       const std::string& representation, json::Value::Object& evidence) {
    if (model.requiresExternalCoordinates())
        throw std::invalid_argument("residual process covariance needs a K-ROM evaluated on the plant state");
    const std::size_t n = model.splitDimension(), d = model.stateDimension();
    const auto A = model.generatorSplit();
    const bool forced = representation != "none" && !representation.empty();
    const auto B = forced ? constantMap(model, representation) : std::vector<double>(n, 0.0);
    const auto window = numbers(spec.at("window"));
    std::vector<double> S(n * n, 0.0);
    double residualEnergy = 0, stepEnergy = 0, horizonUsed = 0;
    std::size_t count = 0;
    json::Value::Array records;
    for (const auto& item : spec.at("from_records").array()) {
        const auto path = item.string();
        const auto record = core::readNpz(core::wide(path));
        const auto t = record.at("time").doubles(), x = record.at("values").doubles();
        if (t.size() < 3 || x.size() != t.size() * d) throw std::invalid_argument("record state history does not match the K-ROM: " + path);
        std::vector<double> u(t.size(), 0.0);
        if (record.contains("control_values")) {
            u = record.at("control_values").doubles();
            if (u.size() != t.size()) throw std::invalid_argument("record needs one scalar control per output sample: " + path);
        }
        const double dt = t[1] - t[0];
        const double h = spec.at("residual_horizon").number(dt);
        const auto step = static_cast<std::size_t>(std::llround(h / dt));
        if (step < 1 || std::abs(step * dt - h) > 1e-6 * h) throw std::invalid_argument("residual_horizon must be a multiple of the record step");
        horizonUsed = h;
        const auto Phi = numerics::matrixExponential(A, n, h);
        const auto Gamma = numerics::zeroOrderHold(A, B, n, 1, h).response;
        std::vector<std::vector<double>> psi(t.size());
        for (std::size_t k = 0; k < t.size(); ++k) psi[k] = model.evaluateSplit(&x[k * d]);
        for (std::size_t k = 0; k + step < t.size(); ++k) {
            if (window.size() == 2 && (t[k] < window[0] || t[k + step] > window[1])) continue;
            std::vector<double> r(n);
            for (std::size_t i = 0; i < n; ++i) {
                double pred = Gamma[i] * u[k];
                for (std::size_t j = 0; j < n; ++j) pred += Phi[i * n + j] * psi[k][j];
                r[i] = psi[k + step][i] - pred;
                residualEnergy += r[i] * r[i];
                stepEnergy += (psi[k + step][i] - psi[k][i]) * (psi[k + step][i] - psi[k][i]);
            }
            for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) S[i * n + j] += r[i] * r[j];
            ++count;
        }
        records.push_back(path);
    }
    if (!count) throw std::invalid_argument("residual process covariance found no samples (check records and window)");
    const double scale = spec.at("scale").number(1.0);
    for (auto& v : S) v *= scale / (double(count) * horizonUsed);
    // keep W positive definite: raise eigenvalues below floor * max
    const auto eig = numerics::symmetricEigensystem(S, n);
    const double top = eig.eigenvalues.back(), floor = spec.at("floor").number(1e-6) * top;
    std::vector<double> W(n * n, 0.0);
    json::Value::Array spectrum;
    for (std::size_t k = 0; k < n; ++k) {
        const double lam = std::max(eig.eigenvalues[k], floor);
        spectrum.push_back(lam);
        for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j)
            W[i * n + j] += lam * eig.eigenvectors[i * n + k] * eig.eigenvectors[j * n + k];
    }
    const double add = spec.at("add").number(0.0);
    for (std::size_t k = 0; k < n; ++k) W[k * n + k] += add;
    evidence["source"] = "k-rom one-step prediction residual";
    evidence["records"] = records; evidence["residual_horizon"] = horizonUsed; evidence["samples"] = double(count);
    evidence["scale"] = scale; evidence["add"] = add; evidence["eigenvalues"] = spectrum;
    evidence["relative_residual"] = stepEnergy > 0 ? std::sqrt(residualEnergy / stepEnergy) : 0.0;
    return W;
}
}  // namespace

KlqeProcessCovariance klqeProcessCovariance(const json::Value& spec, const core::KoopmanModel& model) {
    const std::size_t n = model.splitDimension();
    const auto representation = spec.at("control_representation").string("none");
    const auto& source = spec.at("process_covariance");
    KlqeProcessCovariance out;
    json::Value::Object evidence;
    if (std::holds_alternative<json::Value::Object>(source.storage()))
        out.matrix = residualCovariance(source, model, representation, evidence);
    else out.matrix = plain(source, n);
    const double q2 = spec.at("process_covariance_input_weight").number(0.0);
    if (!(q2 >= 0) || !std::isfinite(q2)) throw std::invalid_argument("process_covariance_input_weight must be nonnegative");
    if (q2 > 0) {
        const auto B = constantMap(model, representation);
        for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) out.matrix[i * n + j] += q2 * B[i] * B[j];
        evidence["input_weight"] = q2;
    }
    if (!evidence.empty()) out.evidence = json::Value(std::move(evidence));
    return out;
}

}  // namespace dynsyslab::workflow
