#pragma once
// Command-line options of dynsyslab_cli: one parsed record, filled by arguments().
#include <cstddef>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::cli {
struct Arguments {
    std::vector<std::wstring> positional;
    std::filesystem::path archiveRoot;
    std::string archiveRootSource{"argument"};  // see CliRuntimePaths.hpp when not given
    bool help{};
    std::string format{"markdown"};  // quickstart, hub, agent-setup: markdown | json
    bool hook{}, peek{};             // hub
    bool remove{}, agentCommands{true};  // agent-setup (--remove, --no-agent-cli)
    std::filesystem::path workspace, home;  // agent-setup --workspace, --home
    std::filesystem::path progress;
    double auditDays{2.0};
    std::vector<std::string> auditSections;
    std::vector<std::string> auditSystems;
    bool auditDerived{};
    std::string auditCampaign;           // archive --campaign ID
    std::filesystem::path report;        // archive audit --report DIR
    std::string compression;             // --compression deflate|lzma|stored (archives written)
    std::string metadataPatch;           // archive edit-metadata --metadata-json JSON
    bool readOnly{}, removeSource{};     // storage register --read-only, storage copy --remove-source
    std::string purpose{"solutions"};    // storage register --purpose
    bool apply{};
    std::string importSystem,importPreset,importCoordinates{"unknown"};
    std::string importLabel{"legacy-dynbatch-import"},importParameters{"{}"};
    double inputScale{1.0};
    bool allowUnknownRecords{};
    std::size_t fourierComponent{};
    std::string fourierWindow{"hann"};
    bool fourierRemoveMean{true};
    bool fourierOneSided{true};
    std::vector<std::pair<std::string, std::string>> params, sets;
    // job groups and data packs
    bool allowVersionMismatch{}, noCleanup{};
    double minFreeGb{-1.0};
    std::filesystem::path referenceMap, output;
};

// Parses argv (argv[0] is skipped). Unknown tokens become positional arguments; --help/-h set
// help. Without --archive-root, archiveRoot is defaultArchiveRoot() (never the working directory).
Arguments arguments(int count, wchar_t** values);
}  // namespace dynsyslab::cli
