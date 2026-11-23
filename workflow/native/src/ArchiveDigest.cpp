#include <dynsyslab/workflow/ArchiveDigest.hpp>

#include "PngCompare.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <set>
#include <string>
#include <variant>
#include <vector>

namespace dynsyslab::workflow {
namespace {

using core::ArrayType;
using core::NpyArray;

bool isPng(const std::filesystem::path& path) {
    auto extension = path.extension().u8string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return extension == ".png";
}

bool endsWith(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool pathLike(const std::string& text) {
    if (text.find(".npz") != std::string::npos || text.find(".json") != std::string::npos) return true;
    return text.size() > 2 && std::isalpha(static_cast<unsigned char>(text[0])) && text[1] == ':'
        && (text[2] == '/' || text[2] == '\\');
}

json::Value withoutVolatile(const json::Value& value) {
    static const std::set<std::string> dropped{
        "created_utc", "updated_utc", "hierarchy_manifest", "batch_id", "system_run_id", "cluster_id",
        "trajectory_id", "parameter_set_key", "source_sha256", "archive_sha256", "dynsyslab_version",
        "source_commit"};
    if (const auto* text = std::get_if<std::string>(&value.storage()))
        return pathLike(*text) ? json::Value("<path>") : value;
    if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        json::Value::Array out;
        for (const auto& item : *array) out.push_back(withoutVolatile(item));
        return json::Value(std::move(out));
    }
    if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        json::Value::Object out;
        for (const auto& [key, item] : *object)
            if (!dropped.count(key)) out[key] = withoutVolatile(item);
        return json::Value(std::move(out));
    }
    return value;
}

// Text form hashed for a JSON text array; other text is used as it is.
std::string canonicalText(const std::string& name, const NpyArray& array) {
    const auto text = array.string();
    if (!endsWith(name, "_json")) return pathLike(text) ? std::string("<path>") : text;
    try { return json::dump(withoutVolatile(json::parse(text))); }
    catch (const std::exception&) { return text; }
}

std::string shapeText(const NpyArray& array) {
    std::string out;
    for (const auto size : array.shape) out += std::to_string(size) + ",";
    return out;
}

bool numeric(const NpyArray& array) {
    return array.type == ArrayType::float64 || array.type == ArrayType::complex128;
}

std::string arrayHash(const std::string& name, const NpyArray& array) {
    core::Sha256Stream digest;
    digest.update(array.descriptor + "|" + shapeText(array) + "|");
    if (array.type == ArrayType::unicode) digest.update(canonicalText(name, array));
    else digest.update(array.data.data(), array.data.size());
    return digest.finish();
}

// View of one float component; complex arrays are compared part by part.
const double* doublesOf(const NpyArray& array, std::size_t& count) {
    count = array.data.size() / sizeof(double);
    return reinterpret_cast<const double*>(array.data.data());
}

json::Value arrayDigest(const std::string& name, const NpyArray& array) {
    json::Value::Object out;
    out["sha256"] = arrayHash(name, array);
    out["dtype"] = array.descriptor;
    json::Value::Array shape;
    for (const auto size : array.shape) shape.emplace_back(static_cast<double>(size));
    out["shape"] = json::Value(std::move(shape));
    if (numeric(array) && !array.data.empty()) {
        std::size_t count = 0;
        const auto* values = doublesOf(array, count);
        double low = HUGE_VAL, high = -HUGE_VAL, squares = 0.0, nonFinite = 0.0;
        for (std::size_t k = 0; k < count; ++k) {
            if (!std::isfinite(values[k])) { nonFinite += 1.0; continue; }
            low = std::min(low, values[k]); high = std::max(high, values[k]);
            squares += values[k] * values[k];
        }
        if (nonFinite < static_cast<double>(count)) {
            out["min"] = low; out["max"] = high; out["l2"] = std::sqrt(squares);
        }
        out["non_finite"] = nonFinite;
        constexpr std::size_t kSamples = 8;
        json::Value::Array samples;
        for (std::size_t k = 0; k < kSamples && count; ++k) {
            const auto value = values[k * (count - 1) / (kSamples - 1)];
            samples.push_back(std::isfinite(value) ? json::Value(value) : json::Value());
        }
        out["samples"] = json::Value(std::move(samples));
    }
    return json::Value(std::move(out));
}

struct Worst {
    double absolute{};
    double relative{};
    std::size_t index{};
    bool structural{};  // NaN pattern or infinity disagreement
};

Worst worstElement(const NpyArray& a, const NpyArray& b) {
    Worst worst;
    std::size_t count = 0, other = 0;
    const auto* x = doublesOf(a, count);
    const auto* y = doublesOf(b, other);
    for (std::size_t k = 0; k < count && k < other; ++k) {
        if (std::isnan(x[k]) && std::isnan(y[k])) continue;
        if (!std::isfinite(x[k]) || !std::isfinite(y[k])) {
            if (x[k] == y[k]) continue;
            worst.structural = true; worst.index = k; continue;
        }
        const double difference = std::abs(x[k] - y[k]);
        if (difference > worst.absolute) {
            worst.absolute = difference; worst.index = k;
            worst.relative = difference / std::max(std::abs(x[k]), std::abs(y[k]));
        }
    }
    return worst;
}

const char* combine(const char* status, const char* next) {
    if (std::string(status) == "failed" || std::string(next) == "failed") return "failed";
    if (std::string(status) == "numerically_equivalent" || std::string(next) == "numerically_equivalent")
        return "numerically_equivalent";
    return "exact";
}

json::Value comparePng(const std::filesystem::path& actual, const std::filesystem::path& reference) {
    json::Value::Object out;
    out["kind"] = "png";
    if (core::sha256File(actual) == core::sha256File(reference)) {
        out["status"] = "exact"; out["bytes_identical"] = true; out["differing_pixels"] = 0.0;
        return json::Value(std::move(out));
    }
    out["bytes_identical"] = false;
    const auto pixels = detail::comparePngPixels(actual, reference);
    out["same_size"] = pixels.sameSize;
    out["differing_pixels"] = static_cast<double>(pixels.differing);
    out["maximum_channel_difference"] = static_cast<double>(pixels.maximumChannelDifference);
    // Identical pixels in a different encoding are the same image.
    out["status"] = pixels.sameSize && pixels.differing == 0 ? "numerically_equivalent" : "failed";
    return json::Value(std::move(out));
}

}  // namespace

json::Value resultDigest(const std::filesystem::path& path) {
    json::Value::Object out;
    if (isPng(path)) {
        out["kind"] = "png"; out["sha256"] = core::sha256File(path);
        return json::Value(std::move(out));
    }
    const auto record = core::readNpz(path);
    json::Value::Object arrays;
    core::Sha256Stream combined;
    for (const auto& [name, array] : record.arrays) {
        if (name == "metadata_json") continue;
        arrays[name] = arrayDigest(name, array);
        combined.update(name + "=" + arrays[name].at("sha256").string() + "\n");
    }
    out["kind"] = "npz"; out["sha256"] = combined.finish(); out["arrays"] = json::Value(std::move(arrays));
    return json::Value(std::move(out));
}

json::Value compareResults(const std::filesystem::path& actual, const std::filesystem::path& reference,
                           const AgreementTolerance& tolerance) {
    if (isPng(actual)) return comparePng(actual, reference);
    const auto a = core::readNpz(actual);
    const auto b = core::readNpz(reference);
    const char* status = "exact";
    json::Value::Object differing;
    json::Value::Array missing, added;
    double worstRelative = 0.0, worstAbsolute = 0.0;
    for (const auto& [name, array] : b.arrays) {
        if (name == "metadata_json") continue;
        if (!a.contains(name)) { missing.emplace_back(name); status = "failed"; }
    }
    for (const auto& [name, array] : a.arrays) {
        if (name == "metadata_json") continue;
        if (!b.contains(name)) {  // written by newer code; shared arrays decide the numerics
            added.emplace_back(name);
            status = combine(status, "numerically_equivalent");
            continue;
        }
        const auto& other = b.at(name);
        if (arrayHash(name, array) == arrayHash(name, other)) continue;
        json::Value::Object item;
        const char* verdict = "failed";
        if (array.shape != other.shape || array.descriptor != other.descriptor) {
            item["reason"] = "shape or type differs";
        } else if (numeric(array)) {
            const auto worst = worstElement(array, other);
            item["maximum_absolute_difference"] = worst.absolute;
            item["relative_difference_at_maximum"] = worst.relative;
            item["flat_index"] = static_cast<double>(
                array.type == ArrayType::complex128 ? worst.index / 2 : worst.index);
            worstAbsolute = std::max(worstAbsolute, worst.absolute);
            worstRelative = std::max(worstRelative, worst.relative);
            if (worst.structural) item["reason"] = "NaN or infinity pattern differs";
            else if (worst.relative <= tolerance.relative || worst.absolute <= tolerance.absolute)
                verdict = "numerically_equivalent";
        } else item["reason"] = array.type == ArrayType::unicode ? "text differs" : "integer data differs";
        item["status"] = verdict;
        status = combine(status, verdict);
        differing[name] = json::Value(std::move(item));
    }
    json::Value::Object out;
    out["kind"] = "npz"; out["status"] = status;
    out["maximum_absolute_difference"] = worstAbsolute;
    out["maximum_relative_difference"] = worstRelative;
    if (!differing.empty()) out["differing_arrays"] = json::Value(std::move(differing));
    if (!missing.empty()) out["missing_arrays"] = json::Value(std::move(missing));
    if (!added.empty()) out["added_arrays"] = json::Value(std::move(added));
    return json::Value(std::move(out));
}

json::Value compareDigests(const json::Value& actual, const json::Value& expected,
                           const AgreementTolerance& tolerance) {
    json::Value::Object out;
    if (actual.at("sha256").string() == expected.at("sha256").string() && !expected.at("sha256").string().empty()) {
        out["status"] = "exact";
        return json::Value(std::move(out));
    }
    if (expected.at("kind").string() == "png") {
        out["status"] = "failed"; out["reason"] = "image bytes differ";
        return json::Value(std::move(out));
    }
    const char* status = "exact";
    json::Value::Object differing;
    for (const auto& [name, want] : expected.at("arrays").object()) {
        const auto& got = actual.at("arrays").at(name);
        if (got.isNull()) { differing[name] = "missing"; status = "failed"; continue; }
        if (got.at("sha256").string() == want.at("sha256").string()) continue;
        const char* verdict = "numerically_equivalent";
        double worst = 0.0;
        if (json::dump(got.at("shape")) != json::dump(want.at("shape")) || want.at("samples").isNull()) verdict = "failed";
        else {
            std::vector<std::pair<double, double>> pairs;
            for (const auto* key : {"min", "max", "l2"}) pairs.emplace_back(got.at(key).number(), want.at(key).number());
            const auto& x = got.at("samples").array();
            const auto& y = want.at("samples").array();
            for (std::size_t k = 0; k < x.size() && k < y.size(); ++k) {
                if (x[k].isNull() != y[k].isNull()) verdict = "failed";
                else if (!x[k].isNull()) pairs.emplace_back(x[k].number(), y[k].number());
            }
            for (const auto& [p, q] : pairs) {
                const double d = std::abs(p - q), scale = std::max(std::abs(p), std::abs(q));
                worst = std::max(worst, scale > 0 ? d / scale : d);
                if (d > tolerance.absolute && d > tolerance.relative * scale) verdict = "failed";
            }
        }
        differing[name] = json::Value(json::Value::Object{{"status", json::Value(verdict)},
                                                          {"worst_relative_in_summary", json::Value(worst)}});
        status = combine(status, verdict);
    }
    if (std::string(status) == "exact") status = "failed";  // the combined hash differed
    out["status"] = status;
    out["basis"] = "stored ranges, norms and samples; full arrays were not compared";
    out["differing_arrays"] = json::Value(std::move(differing));
    return json::Value(std::move(out));
}

}  // namespace dynsyslab::workflow
