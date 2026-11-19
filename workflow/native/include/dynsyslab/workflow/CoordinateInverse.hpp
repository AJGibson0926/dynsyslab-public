#pragma once

#include <filesystem>

namespace dynsyslab::workflow {

// Apply the inverse a coordinate change retained (exact, or a fitted readout)
// to INPUT.npz ("values": samples x coordinates, optional "time") and write
// the previous coordinates to OUTPUT.npz. Throws when no inverse is retained.
void invertCoordinates(const std::filesystem::path& coordinateArchive,
                       const std::filesystem::path& input, const std::filesystem::path& output);

}  // namespace dynsyslab::workflow
