#pragma once

#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>

namespace dynsyslab::workflow {

// Canonical identity of a result, independent of wall-clock names and absolute paths.
//
// NPZ: one SHA-256 per array over its type, shape and raw bytes, plus range, norm and a few
// evenly spaced samples, so a later mismatch can say where and how large. metadata_json is
// excluded (it holds times and paths); other *_json arrays are hashed after volatile fields
// (times, run ids, file paths) are removed. PNG: SHA-256 of the file bytes.
[[nodiscard]] json::Value resultDigest(const std::filesystem::path& path);

struct AgreementTolerance {
    double relative{1e-9};
    double absolute{1e-12};
};

// Status is "exact", "numerically_equivalent" or "failed". Exact means every canonical hash
// matches. Equivalent means the hashes differ but all values agree within tolerance.
// compareResults reads both files and reports the worst element of each differing array;
// PNGs compare by bytes first, then decoded pixels.
[[nodiscard]] json::Value compareResults(const std::filesystem::path& actual,
                                         const std::filesystem::path& reference,
                                         const AgreementTolerance& tolerance = {});

// The same categories from two stored digests, when the reference data no longer exists.
// Equivalence here rests only on ranges, norms and samples, and the report says so.
[[nodiscard]] json::Value compareDigests(const json::Value& actual, const json::Value& expected,
                                         const AgreementTolerance& tolerance = {});

}  // namespace dynsyslab::workflow
