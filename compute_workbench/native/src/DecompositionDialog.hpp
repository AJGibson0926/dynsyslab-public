#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

struct DecompositionSettings {
    std::string dmd;
    std::string havok;
    std::string edmd;
    std::string svd;
    std::string pca;
    std::string pod;
    std::string baggedPca;
    std::string robustPca;
    std::string hankelPod;
};

// Typed front end for DMD, HAVOK, EDMD, and the PCA/POD/SVD family.
std::optional<DecompositionSettings> showDecompositionDialog(
    HWND owner, const std::filesystem::path& trajectorySource,
    const DecompositionSettings& current, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
