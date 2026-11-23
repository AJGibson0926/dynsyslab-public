#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

namespace dynsyslab::core {

// Lower-case SHA-256 of the exact file bytes. Kept in the archive layer so
// imported evidence has the same durable identity regardless of UI or adapter.
[[nodiscard]] std::string sha256Text(const std::string& text);
[[nodiscard]] std::string sha256File(const std::filesystem::path& path);

// Incremental SHA-256, for digests of data too large to copy into one buffer.
class Sha256Stream {
public:
    void update(const void* data, std::size_t size);
    void update(const std::string& text) { update(text.data(), text.size()); }
    [[nodiscard]] std::string finish();  // lower-case hex; the stream is then spent

private:
    std::array<std::uint32_t, 8> hash_{0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                       0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::array<std::uint8_t, 64> block_{};
    std::size_t filled_{};
    std::uint64_t total_{};
};

}  // namespace dynsyslab::core
