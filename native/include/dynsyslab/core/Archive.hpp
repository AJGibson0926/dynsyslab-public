#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::core {

enum class ArrayType {
    float64, float32, int64, int32, uint64, uint32, uint8, boolean, complex128, unicode, unknown
};

struct NpyArray {
    ArrayType type{ArrayType::unknown};
    std::string descriptor;
    std::vector<std::size_t> shape;
    std::vector<std::uint8_t> data;

    [[nodiscard]] std::size_t elementCount() const;
    [[nodiscard]] std::vector<double> doubles() const;
    [[nodiscard]] std::vector<std::complex<double>> complexDoubles() const;
    [[nodiscard]] std::string string() const;

    static NpyArray fromDoubles(std::vector<std::size_t> shape,
                                const std::vector<double>& values);
    static NpyArray fromDoubles(std::vector<std::size_t> shape,
                                const double* values, std::size_t count);
    static NpyArray fromInt64s(std::vector<std::size_t> shape,
                               const std::vector<std::int64_t>& values);
    static NpyArray fromComplex(std::vector<std::size_t> shape,
                                const std::vector<std::complex<double>>& values);
    static NpyArray fromString(const std::string& value);
};

struct NpzRecord {
    std::map<std::string, NpyArray> arrays;

    [[nodiscard]] bool contains(const std::string& name) const;
    [[nodiscard]] const NpyArray& at(const std::string& name) const;
    [[nodiscard]] std::string metadataJson() const;
};

// Lightweight archive index for catalog/UI discovery.  It enumerates member
// names and inflates only metadata_json, leaving large scientific arrays lazy.
struct NpzIndex {
    std::vector<std::string> arrays;
    std::string metadataJson;
    [[nodiscard]] bool contains(const std::string& name) const;
};

// Reads an archive. A recorded archive (RecordedArchive below) is replayed: its arrays are
// recomputed from the record, and its stored metadata is kept.
[[nodiscard]] NpzRecord readNpz(const std::filesystem::path& path);
// Array names and metadata only; a recorded archive lists the arrays its replay produces.
[[nodiscard]] NpzIndex readNpzIndex(const std::filesystem::path& path);
// The file exactly as stored (no replay).
[[nodiscard]] NpzRecord readNpzStored(const std::filesystem::path& path);
[[nodiscard]] NpzIndex readNpzIndexStored(const std::filesystem::path& path);

// Recorded archives. A derived archive may store how to recompute its arrays instead of the arrays:
// its metadata then holds "recorded_archive": {"schema", "arrays": [names], ...}. Cheap,
// deterministic derivations (e.g. coordinate changes) are recorded this way, so a chain of them
// keeps one copy of the data, at its root. The library that writes a schema registers its replay;
// readNpz calls it with the stored record and returns the recomputed archive.
using RecordedArchiveReader = std::function<NpzRecord(const std::filesystem::path& path, NpzRecord stored)>;
void registerRecordedArchiveReader(const std::string& schema, RecordedArchiveReader reader);
[[nodiscard]] bool isRecordedArchive(const std::string& metadataJson);
void writeNpzAtomic(const std::filesystem::path& path, const NpzRecord& record);
// Atomically create a completed archive; fail if the destination already exists.
void writeNpzAtomicNew(const std::filesystem::path& path, const NpzRecord& record);
// Writes a ZIP of text entries in the given order, with fixed timestamps (byte-reproducible).
void writeTextZip(const std::filesystem::path& path,
                  const std::vector<std::pair<std::string, std::string>>& entries);

}  // namespace dynsyslab::core
