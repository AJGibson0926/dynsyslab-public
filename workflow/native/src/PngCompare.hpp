#pragma once

#include <cstddef>
#include <filesystem>

namespace dynsyslab::workflow::detail {

struct PixelComparison {
    bool sameSize{};
    std::size_t differing{};
    int maximumChannelDifference{};
};

// Decodes both images to 8-bit BGRA with the Windows Imaging Component and compares pixels.
[[nodiscard]] PixelComparison comparePngPixels(const std::filesystem::path& first,
                                               const std::filesystem::path& second);

}  // namespace dynsyslab::workflow::detail
