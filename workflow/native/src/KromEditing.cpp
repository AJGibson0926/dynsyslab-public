// Post-processing methods that derive a new K-ROM archive from an existing one.
//
// krom_recenter  {"target_state":[...]}: psi_new(x) = psi(x) - psi(x*), by exact re-expansion of
//   the polynomial eigenfunctions about x* with the constant dropped. Readout intercepts absorb
//   the constant, h0' = h0 + C psi(x*); gradients, and so the input map, are unchanged.
// krom_readout   {"readouts":[{"name":N,"kind":"lifted_split"}]}: adds a named linear readout.
//   lifted_split is the real split of the lifted state, [Re psi_1, Im psi_1, Re psi_2, ...],
//   pairing each mode with its complex conjugate.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <complex>
#include <cstdint>
#include <map>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;
using Complex = std::complex<double>;

double binomial(std::int64_t n, std::int64_t k) {
    double out = 1.0;
    for (std::int64_t j = 1; j <= k; ++j) out = out * static_cast<double>(n - k + j) / static_cast<double>(j);
    return std::round(out);
}

// Complex-by-real product written as (a + ib)(w + 0i).
Complex timesReal(const Complex& a, double w) {
    return {a.real() * w - a.imag() * 0.0, a.real() * 0.0 + a.imag() * w};
}

std::filesystem::path publish(c::NpzRecord record, const std::filesystem::path& source, const std::filesystem::path& root,
                              const char* folder, const char* suffix, json::Value provenance, const std::string& label) {
    auto metadata = json::parse(record.metadataJson());
    metadata.at("label") = json::Value(label);
    metadata.at("source_krom") = json::Value(source.generic_u8string());
    metadata.at("krom_edit") = std::move(provenance);
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(metadata));
    const auto directory = root / L"derived" / c::wide(folder);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(source.stem().u8string() + suffix + timestamp()) + ".npz");
    c::writeNpzAtomicNew(path, record);
    return path;
}

}  // namespace

std::filesystem::path runKromRecenter(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                      const std::filesystem::path& root, const std::string& label) {
    if (sources.size() != 1) throw std::invalid_argument("krom_recenter takes one K-ROM source");
    auto record = c::readNpz(sources.front());
    std::vector<double> target;
    for (const auto& value : settings.at("target_state").array()) target.push_back(value.number());
    const auto center = record.at("state_center").doubles();
    const auto scale = record.at("state_scale").doubles();
    const auto dimension = center.size();
    if (target.size() != dimension) throw std::invalid_argument("krom_recenter: target_state dimension differs from the K-ROM");
    const auto& exponentArray = record.at("polynomial_exponents");
    const auto& coefficientArray = record.at("coefficients");
    if (coefficientArray.shape.size() != 2 || exponentArray.shape.size() != 2 || exponentArray.shape[1] != dimension
        || exponentArray.shape[0] != coefficientArray.shape[1])
        throw std::invalid_argument("krom_recenter needs a polynomial K-ROM");
    const auto modes = coefficientArray.shape[0], features = coefficientArray.shape[1];
    const auto coefficients = coefficientArray.complexDoubles();
    std::vector<std::int64_t> exponents;
    for (const double value : exponentArray.doubles()) exponents.push_back(std::llround(value));
    std::vector<double> shift(dimension);
    for (std::size_t k = 0; k < dimension; ++k) shift[k] = (target[k] - center[k]) / scale[k];
    std::map<std::vector<std::int64_t>, std::size_t> index;
    for (std::size_t i = 0; i < features; ++i)
        index[std::vector<std::int64_t>(exponents.begin() + std::ptrdiff_t(i * dimension),
                                        exponents.begin() + std::ptrdiff_t((i + 1) * dimension))] = i;
    std::vector<Complex> updated(coefficients.size()), constant(modes);
    for (std::size_t i = 0; i < features; ++i) {
        const std::vector<std::int64_t> power(exponents.begin() + std::ptrdiff_t(i * dimension),
                                              exponents.begin() + std::ptrdiff_t((i + 1) * dimension));
        std::vector<std::int64_t> lower(dimension, 0);
        while (true) {  // every monomial below this one in each variable
            double weight = 1.0;
            for (std::size_t k = 0; k < dimension; ++k) weight = weight * binomial(power[k], lower[k]);
            for (std::size_t k = 0; k < dimension; ++k)
                weight = weight * std::pow(shift[k], static_cast<double>(power[k] - lower[k]));
            const bool isConstant = std::all_of(lower.begin(), lower.end(), [](std::int64_t v) { return v == 0; });
            if (isConstant) {
                for (std::size_t m = 0; m < modes; ++m) constant[m] += timesReal(coefficients[m * features + i], weight);
            } else {
                const auto found = index.find(lower);
                if (found == index.end()) throw std::invalid_argument("krom_recenter: the exponent set is not downward closed");
                for (std::size_t m = 0; m < modes; ++m)
                    updated[m * features + found->second] += timesReal(coefficients[m * features + i], weight);
            }
            std::size_t digit = dimension;
            while (digit > 0) {
                --digit;
                if (lower[digit] < power[digit]) { ++lower[digit]; break; }
                lower[digit] = 0;
                if (digit == 0) { digit = dimension + 1; break; }
            }
            if (digit > dimension) break;
        }
    }
    const auto shape = coefficientArray.shape;
    record.arrays["coefficients"] = c::NpyArray::fromComplex(shape, updated);
    record.arrays["state_center"] = c::NpyArray::fromDoubles({dimension}, target);
    const std::string prefix = "reconstruction__", offsetPrefix = "reconstruction_offset__";
    std::vector<std::string> names;
    for (const auto& item : record.arrays)
        if (item.first.rfind(prefix, 0) == 0) names.push_back(item.first.substr(prefix.size()));
    for (const auto& name : names) {
        if (!record.contains(offsetPrefix + name)) continue;
        const auto& readout = record.at(prefix + name);
        const auto rows = readout.shape.at(0), columns = readout.shape.at(1);
        if (columns != modes) throw std::invalid_argument("krom_recenter: readout width differs from the mode count");
        const auto matrix = readout.complexDoubles();
        const auto offsetShape = record.at(offsetPrefix + name).shape;
        auto offsets = record.at(offsetPrefix + name).doubles();
        for (std::size_t r = 0; r < rows; ++r) {
            Complex sum{};
            for (std::size_t k = 0; k < columns; ++k) sum += matrix[r * columns + k] * constant[k];
            offsets[r] += sum.real();
        }
        record.arrays[offsetPrefix + name] = c::NpyArray::fromDoubles(offsetShape, offsets);
    }
    json::Value::Array removed, star;
    for (const auto& value : constant) removed.push_back(json::Value(json::Value::Array{json::Value(value.real()), json::Value(value.imag())}));
    for (const double value : target) star.emplace_back(value);
    return publish(std::move(record), sources.front(), root, "koopman-krom-edit", "__recentred__",
                   json::Value(json::Value::Object{{"operation", json::Value("krom_recenter")},
                       {"target_state", json::Value(std::move(star))}, {"psi_at_target_removed", json::Value(std::move(removed))}}), label);
}

std::filesystem::path runKromReadout(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                     const std::filesystem::path& root, const std::string& label) {
    if (sources.size() != 1) throw std::invalid_argument("krom_readout takes one K-ROM source");
    const auto& requested = settings.at("readouts").array();
    if (requested.empty()) throw std::invalid_argument("krom_readout lists no readouts");
    auto record = c::readNpz(sources.front());
    const auto eigenvalues = record.at("continuous_eigenvalues").complexDoubles();
    const auto m = eigenvalues.size();
    auto names = record.contains("reconstruction_names") ? record.at("reconstruction_names").string() : std::string();
    for (const auto& readout : requested) {
        const auto name = readout.at("name").string();
        if (name.empty() || name.find('\n') != std::string::npos) throw std::invalid_argument("krom_readout needs a single-line name");
        if (record.contains("reconstruction__" + name)) throw std::invalid_argument("krom_readout: readout already exists: " + name);
        if (readout.at("kind").string() != "lifted_split")
            throw std::invalid_argument("krom_readout: unknown readout kind: " + readout.at("kind").string());
        const double relative = readout.at("pair_relative_tolerance").number(1e-5);
        const double absolute = readout.at("pair_absolute_tolerance").number(1e-8);
        std::vector<Complex> rows(2 * m * m);
        for (std::size_t k = 0; k < m; ++k) {
            std::size_t conjugate = m;
            const auto wanted = std::conj(eigenvalues[k]);
            for (std::size_t i = 0; i < m && conjugate == m; ++i)
                if (i != k && std::abs(eigenvalues[i] - wanted) <= absolute + relative * std::abs(wanted)) conjugate = i;
            auto* re = &rows[2 * k * m];
            auto* im = &rows[(2 * k + 1) * m];
            if (conjugate < m) {
                re[k] = 0.5; im[k] = Complex(-0.0, -0.5);
                re[conjugate] = 0.5; im[conjugate] = Complex(0.0, 0.5);
            } else re[k] = 1.0;  // a real mode: Re = psi, Im = 0
        }
        record.arrays["reconstruction__" + name] = c::NpyArray::fromComplex({2 * m, m}, rows);
        record.arrays["reconstruction_offset__" + name] = c::NpyArray::fromDoubles({2 * m}, std::vector<double>(2 * m, 0.0));
        names = names.empty() ? name : names + "\n" + name;
    }
    record.arrays["reconstruction_names"] = c::NpyArray::fromString(names);
    return publish(std::move(record), sources.front(), root, "koopman-krom-edit", "__readout__",
                   json::Value(json::Value::Object{{"operation", json::Value("krom_readout")}, {"readouts", settings.at("readouts")}}), label);
}

}  // namespace dynsyslab::workflow::postprocess_detail
