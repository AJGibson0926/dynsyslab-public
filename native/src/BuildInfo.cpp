#include <dynsyslab/core/BuildInfo.hpp>

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Version.hpp>

#include <utility>

#ifndef DYNSYSLAB_SOURCE_COMMIT
#define DYNSYSLAB_SOURCE_COMMIT "unknown"
#endif
#ifndef DYNSYSLAB_SOURCE_MODIFIED
#define DYNSYSLAB_SOURCE_MODIFIED 0
#endif
#ifndef DYNSYSLAB_COMPILER_ID
#define DYNSYSLAB_COMPILER_ID "unknown"
#endif
#ifndef DYNSYSLAB_BUILD_TYPE
#define DYNSYSLAB_BUILD_TYPE "unknown"
#endif

namespace dynsyslab::core {
namespace {
thread_local std::string activeProvenance;

std::string platformName() {
#if defined(_WIN64)
    return "windows-x64";
#elif defined(_WIN32)
    return "windows-x86";
#elif defined(__linux__)
    return "linux";
#elif defined(__APPLE__)
    return "macos";
#else
    return "unknown";
#endif
}
}  // namespace

const BuildIdentity& buildIdentity() {
    static const BuildIdentity identity{std::string(version), DYNSYSLAB_SOURCE_COMMIT, DYNSYSLAB_SOURCE_MODIFIED != 0,
                                        DYNSYSLAB_COMPILER_ID, DYNSYSLAB_BUILD_TYPE, platformName()};
    return identity;
}

std::string buildIdentityJson() {
    const auto& identity = buildIdentity();
    return json::dump(json::Value(json::Value::Object{
        {"dynsyslab_version", json::Value(identity.version)}, {"source_commit", json::Value(identity.sourceCommit)},
        {"source_modified", json::Value(identity.sourceModified)}, {"compiler", json::Value(identity.compiler)},
        {"build_type", json::Value(identity.buildType)}, {"platform", json::Value(identity.platform)}}));
}

ScopedRunProvenance::ScopedRunProvenance(std::string provenanceJson)
    : previous_(std::exchange(activeProvenance, std::move(provenanceJson))) {}

ScopedRunProvenance::~ScopedRunProvenance() { activeProvenance = std::move(previous_); }

const std::string& activeRunProvenance() { return activeProvenance; }

}  // namespace dynsyslab::core
