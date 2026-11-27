#pragma once
// Storage locations of an archive root (dynsyslab-cli storage ...). A registry names extra
// folders (other disks, shared drives) where archives may be kept; it never moves anything by
// itself. Files move only by verified copy: the copy is hashed (SHA-256, FIPS 180-4) against its source before
// it is published, an existing different file is never overwritten, and the source is removed only
// when asked and only after verification. A campaign relocation plan records the copies of one
// campaign's archives to a registered location, keeping their layout relative to the archive root.
//
// Files (JSON, written atomically):
//   <root>/management/storage_locations.json        dynsyslab.storage-registry/v1
//   <root>/management/relocations/campaign-relocation__<campaign>__<stamp>.json
//                                                    dynsyslab.campaign-relocation/v1
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::core {

inline constexpr const char* kStorageRegistrySchema = "dynsyslab.storage-registry/v1";
inline constexpr const char* kCampaignRelocationSchema = "dynsyslab.campaign-relocation/v1";

struct StorageLocation {
    std::string name;
    std::filesystem::path root;
    bool writable{true};
    std::string purpose{"solutions"};
};

[[nodiscard]] std::filesystem::path storageRegistryPath(const std::filesystem::path& archiveRoot);
// Registered locations; "primary" (the archive root itself) is always first.
[[nodiscard]] std::vector<StorageLocation> storageLocations(const std::filesystem::path& archiveRoot);
// Adds or replaces a location and creates its folder. The name "primary" is reserved.
void registerStorageLocation(const std::filesystem::path& archiveRoot, const StorageLocation& location);

struct VerifiedCopy {
    std::filesystem::path source;
    std::filesystem::path destination;
    std::string sha256;
    std::uintmax_t bytes{};
    bool sourceRemoved{};
    bool alreadyAtDestination{};  // the source already lives under the destination root
    std::string copiedUtc;
    std::string relativePath;     // set by relocation copies
};

// Copies `source` into `destinationFolder` under its own name, verifying size and SHA-256.
// An identical existing copy is accepted; a different one is an error. Active *.partial files are
// refused.
VerifiedCopy copyVerified(const std::filesystem::path& source, const std::filesystem::path& destinationFolder,
                          bool removeSource = false);

// The archives of one campaign: its manifests in <root>/campaigns (with every hierarchy manifest
// and solution archive they list) and every archive under <root>/solutions whose metadata has
// that campaign_id, plus the shared time grids those archives reference. Existing files only, sorted.
[[nodiscard]] std::vector<std::filesystem::path> campaignArchives(const std::filesystem::path& archiveRoot,
                                                                  const std::string& campaignId);

// Creates a relocation plan for a campaign to a registered writable location. Moves nothing.
std::filesystem::path createRelocationPlan(const std::filesystem::path& archiveRoot,
                                           const std::string& campaignId, const std::string& destination);
// Copies one completed file for a plan (sources are kept) and records the copy in the plan. A file
// already under the destination root is recorded without copying.
VerifiedCopy appendRelocationCopy(const std::filesystem::path& plan, const std::filesystem::path& source);

[[nodiscard]] std::string storageLocationsJson(const std::vector<StorageLocation>& locations);
[[nodiscard]] std::string verifiedCopyJson(const VerifiedCopy& copy);

}  // namespace dynsyslab::core
