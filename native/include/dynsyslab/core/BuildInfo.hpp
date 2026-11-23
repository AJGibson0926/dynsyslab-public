#pragma once

#include <string>

namespace dynsyslab::core {

// Identity of this compiled suite, fixed when the build is configured.
struct BuildIdentity {
    std::string version;
    std::string sourceCommit;   // full git hash, or "unknown" outside a checkout
    bool sourceModified{};      // the checkout had uncommitted changes when configured
    std::string compiler;       // for example "MSVC 19.44.35211.0"
    std::string buildType;
    std::string platform;       // for example "windows-x64"
};

[[nodiscard]] const BuildIdentity& buildIdentity();
[[nodiscard]] std::string buildIdentityJson();

// While alive, every NPZ this thread publishes records {"run_provenance":<json>} in its
// metadata_json. Only metadata changes; numerical arrays are written untouched.
class ScopedRunProvenance {
public:
    explicit ScopedRunProvenance(std::string provenanceJson);
    ~ScopedRunProvenance();
    ScopedRunProvenance(const ScopedRunProvenance&) = delete;
    ScopedRunProvenance& operator=(const ScopedRunProvenance&) = delete;

private:
    std::string previous_;
};

[[nodiscard]] const std::string& activeRunProvenance();

}  // namespace dynsyslab::core
