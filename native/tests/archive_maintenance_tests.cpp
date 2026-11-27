// Archive storage maintenance: codecs, shared time grids, recompression, metadata edits, audit
// campaign scope and report files, storage locations with verified copy and campaign relocation,
// and the campaign catalog (docs/ARCHIVE_MANAGEMENT.md).
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveAudit.hpp>
#include <dynsyslab/core/ArchiveMaintenance.hpp>
#include <dynsyslab/core/CampaignCatalog.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/StorageLocations.hpp>
#include <dynsyslab/core/TimeGrid.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace fs = std::filesystem;
namespace core = dynsyslab::core;
namespace json = dynsyslab::core::json;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

fs::path freshRoot(const char* name) {
    const auto root = fs::temp_directory_path() / (std::string(name) + "-"
        + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(root);
    return root;
}

core::NpzRecord solution(double scale, const std::string& campaign = {}) {
    std::vector<double> time(200), values(400);
    for (std::size_t index = 0; index < time.size(); ++index) {
        time[index] = 0.01 * static_cast<double>(index);
        values[2 * index] = scale * std::sin(time[index]);
        values[2 * index + 1] = scale * std::cos(time[index]);
    }
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({time.size()}, time);
    record.arrays["values"] = core::NpyArray::fromDoubles({time.size(), 2}, values);
    std::string metadata = R"({"schema":"dynsyslab.solution/v1","system":"duffing")";
    if (!campaign.empty()) metadata += R"(,"campaign_id":")" + campaign + "\"";
    record.arrays["metadata_json"] = core::NpyArray::fromString(metadata + "}");
    return record;
}

bool sameArrays(const core::NpzRecord& a, const core::NpzRecord& b, const std::string& skip) {
    for (const auto& [name, array] : a.arrays) {
        if (name == skip) continue;
        if (!b.contains(name)) return false;
        const auto& other = b.at(name);
        if (array.descriptor != other.descriptor || array.shape != other.shape || array.data != other.data) return false;
    }
    return true;
}

void codecs() {
    const auto root = freshRoot("dynsyslab-codecs");
    const auto record = solution(1.0);
    for (const char* name : {"deflate", "lzma", "stored"}) {
        const auto codec = core::parseArchiveCompression(name);
        require(core::archiveCompressionName(codec) == name, "codec name");
        const auto path = root / (std::string(name) + ".npz");
        core::writeNpz(path, record, {codec, false, false});
        require(sameArrays(record, core::readNpz(path), {}), "codec round trip");
    }
    require(fs::file_size(root / "lzma.npz") < fs::file_size(root / "stored.npz"), "lzma compresses");
    bool rejected = false;
    try { (void)core::parseArchiveCompression("zip"); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "unknown codec rejected");
    core::setArchiveCompression(core::ArchiveCompression::lzma);
    require(core::archiveCompression() == core::ArchiveCompression::lzma, "process codec setting");
    core::writeNpzAtomic(root / "setting.npz", record);
    require(fs::file_size(root / "setting.npz") == fs::file_size(root / "lzma.npz"), "setting used by writeNpzAtomic");
    core::setArchiveCompression(core::ArchiveCompression::deflate);
    fs::remove_all(root);
}

void timeGrids() {
    const auto root = freshRoot("dynsyslab-time-grids");
    const auto folder = root / "solutions" / "duffing";
    fs::create_directories(folder);
    std::vector<fs::path> paths{folder / "a.npz", folder / "b.npz", folder / "c.npz"};
    const core::NpzRecord originals[]{solution(1.0), solution(2.0), solution(3.0)};
    for (std::size_t index = 0; index < paths.size(); ++index) core::writeNpzAtomicNew(paths[index], originals[index]);
    auto single = solution(4.0);
    single.arrays["time"] = core::NpyArray::fromDoubles({200}, std::vector<double>(200, 1.0));
    core::writeNpzAtomicNew(folder / "single.npz", single);
    paths.push_back(folder / "single.npz");

    const auto dry = core::shareTimeGrids(root, paths, false);
    require(dry.dryRun && dry.grids.size() == 1 && dry.files.size() == 3, "dry-run groups the three shared meshes");
    require(core::readNpzStored(paths[0]).contains("time"), "dry run changes nothing");
    const auto applied = core::shareTimeGrids(root, paths, true);
    require(applied.files.size() == 3, "three archives share one grid");
    for (const auto& item : applied.files) require(item.rewritten && item.error.empty(), "grid sharing rewrote");
    require(fs::is_regular_file(applied.grids.front()), "grid file written");
    for (std::size_t index = 0; index < 3; ++index) {
        require(!core::readNpzStored(paths[index]).contains("time"), "shared archive stores no time");
        require(sameArrays(originals[index], core::readNpz(paths[index]), "metadata_json"), "readNpz restores time");
        require(core::readNpzIndex(paths[index]).contains("time"), "index lists the shared time");
    }
    require(core::readNpzStored(paths[3]).contains("time"), "a singleton mesh stays with its archive");
    require(core::shareTimeGrids(root, paths, true).files.empty(), "sharing is idempotent");

    // A relocated tree still reads: the recorded root is gone, the grid is found above the archive.
    const auto moved = freshRoot("dynsyslab-time-grids-moved");
    fs::copy(root, moved / "tree", fs::copy_options::recursive);
    fs::remove_all(root);
    const auto relocated = moved / "tree" / "solutions" / "duffing" / "a.npz";
    require(sameArrays(originals[0], core::readNpz(relocated), "metadata_json"), "relocated archive reads its grid");

    // A Python-style archive (reference, no sample count) and the guards.
    auto python = originals[1];
    python.arrays.erase("time");
    python.arrays["metadata_json"] = core::NpyArray::fromString(
        R"({"time_grid_ref":")"
        + json::parse(core::readNpzIndexStored(relocated).metadataJson).at("time_grid_ref").string() + "\"}");
    const auto pythonPath = moved / "tree" / "solutions" / "python.npz";
    core::writeNpz(pythonPath, python, {core::ArchiveCompression::lzma, false, false});
    require(sameArrays(originals[1], core::readNpz(pythonPath), "metadata_json"), "reference-only metadata reads");
    python.arrays["metadata_json"] = core::NpyArray::fromString(R"({"time_grid_ref":"../escape.npz"})");
    core::writeNpz(pythonPath, python, {core::ArchiveCompression::deflate, true, false});
    bool rejected = false;
    try { (void)core::readNpz(pythonPath); } catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "a reference leaving its root is rejected");
    fs::remove_all(moved);
}

void recompression() {
    const auto root = freshRoot("dynsyslab-recompress");
    const auto path = root / "solutions" / "a.npz";
    const auto record = solution(1.0);
    core::writeNpz(path, record, {core::ArchiveCompression::stored, false, false});
    const auto before = fs::file_size(path);
    const auto storedMetadata = core::readNpzStored(path).metadataJson();
    const auto dry = core::recompressArchives({path}, core::ArchiveCompression::lzma, false);
    require(dry.dryRun && !dry.files.front().rewritten && fs::file_size(path) == before, "recompress dry run");
    const auto result = core::recompressArchives({path}, core::ArchiveCompression::lzma, true);
    require(result.files.front().rewritten && result.files.front().error.empty(), "recompress applied");
    require(fs::file_size(path) < before, "lzma is smaller than stored");
    const auto after = core::readNpzStored(path);
    require(sameArrays(record, after, {}) && after.metadataJson() == storedMetadata, "recompression is lossless");
    require(core::archiveMaintenanceJson(result).find("dynsyslab.archive-recompress/v1") != std::string::npos,
            "recompress JSON contract");
    fs::remove_all(root);
}

void metadataEdits() {
    const auto root = freshRoot("dynsyslab-metadata-edit");
    const auto source = root / "a.npz";
    core::writeNpzAtomicNew(source, solution(1.0));
    const auto digest = core::sha256File(source);
    const auto output = core::editArchiveMetadata(source, root / "b.npz", R"({"label":"edited","system":null})");
    require(core::sha256File(source) == digest, "the source archive is unchanged");
    const auto edited = core::readNpz(output);
    const auto metadata = json::parse(edited.metadataJson());
    require(metadata.at("label").string() == "edited" && metadata.at("system").isNull(), "patch applied");
    require(metadata.at("metadata_revision").number() == 1.0, "revision counted");
    require(metadata.at("metadata_edits").array().front().at("source_sha256").string() == digest, "edit provenance");
    require(sameArrays(solution(1.0), edited, "metadata_json"), "arrays copied unchanged");
    bool exists = false, protectedKey = false;
    try { (void)core::editArchiveMetadata(source, output, R"({"a":1})"); } catch (const std::invalid_argument&) { exists = true; }
    try { (void)core::editArchiveMetadata(source, root / "c.npz", R"({"time_grid_ref":"x"})"); }
    catch (const std::invalid_argument&) { protectedKey = true; }
    require(exists && protectedKey && !fs::exists(root / "c.npz"), "edit guards");
    fs::remove_all(root);
}

void auditCampaignAndReport() {
    const auto root = freshRoot("dynsyslab-audit-campaign");
    core::writeNpzAtomicNew(root / "solutions" / "a.npz", solution(1.0, "campaign-1"));
    core::writeNpzAtomicNew(root / "solutions" / "b.npz", solution(1.0, "campaign-1"));
    core::writeNpzAtomicNew(root / "solutions" / "c.npz", solution(1.0, "campaign-2"));
    core::ArchiveAuditScope scope;
    scope.modifiedWithinDays = -1.0;
    require(core::auditArchives(root, scope).scannedFiles == 3, "unscoped audit");
    scope.campaignId = "campaign-1";
    const auto report = core::auditArchives(root, scope);
    require(report.scannedFiles == 2 && report.exactDuplicateGroups.size() == 1, "campaign-scoped audit");
    const auto files = core::writeArchiveAuditReport(report, scope, root / "reports");
    require(fs::is_regular_file(files.json) && fs::is_regular_file(files.markdown), "report files");
    std::string text;
    {
        std::ifstream markdown(files.markdown);
        text.assign(std::istreambuf_iterator<char>(markdown), std::istreambuf_iterator<char>());
    }
    require(text.find("campaign-1") != std::string::npos && text.find("Group 1: 2 files") != std::string::npos,
            "markdown report");
    require(json::load(files.json.wstring()).at("scope").at("campaign_id").string() == "campaign-1", "json report");
    fs::remove_all(root);
}

void storageAndRelocation() {
    const auto root = freshRoot("dynsyslab-storage");
    const auto second = freshRoot("dynsyslab-storage-second");
    require(core::storageLocations(root).size() == 1, "primary only");
    core::registerStorageLocation(root, {"cold", second / "cold", true, "solutions"});
    core::registerStorageLocation(root, {"archive-disk", second / "readonly", false, "solutions"});
    const auto locations = core::storageLocations(root);
    require(locations.size() == 3 && locations.front().name == "primary", "registered locations");
    require(fs::is_directory(second / "cold"), "location folder created");
    bool reserved = false;
    try { core::registerStorageLocation(root, {"primary", second, true, "x"}); } catch (const std::invalid_argument&) { reserved = true; }
    require(reserved, "primary is reserved");

    const auto file = root / "solutions" / "duffing" / "a.npz";
    core::writeNpzAtomicNew(file, solution(1.0, "camp"));
    const auto copy = core::copyVerified(file, second / "copies");
    require(fs::is_regular_file(copy.destination) && copy.sha256 == core::sha256File(file), "verified copy");
    (void)core::copyVerified(file, second / "copies");  // an identical copy is accepted
    core::writeNpzAtomicNew(second / "other" / "a.npz", solution(2.0));
    bool conflict = false;
    try { (void)core::copyVerified(file, second / "other"); } catch (const std::runtime_error&) { conflict = true; }
    require(conflict, "a different destination file is never overwritten");

    fs::create_directories(root / "campaigns");
    std::ofstream(root / "campaigns" / "camp.campaign.json")
        << R"({"schema":"dynsyslab.campaign/v3","id":"camp","label":"sweep","kind":"control_sweep","points":[)"
        << R"({"index":0,"hierarchy_manifest":"","solution_archives":[")" << file.generic_u8string()
        << R"(","missing.npz"]}]})";
    require(core::campaignArchives(root, "camp").size() == 2, "campaign manifest and solution");
    const auto shared = root / "solutions" / "duffing" / "b.npz";
    core::writeNpzAtomicNew(shared, solution(3.0, "camp"));
    (void)core::shareTimeGrids(root, {file, shared}, true);
    const auto members = core::campaignArchives(root, "camp");
    require(members.size() == 4, "campaign manifest, solutions and their shared time grid");
    bool notWritable = false;
    try { (void)core::createRelocationPlan(root, "camp", "archive-disk"); } catch (const std::invalid_argument&) { notWritable = true; }
    require(notWritable, "read-only locations are refused");
    const auto plan = core::createRelocationPlan(root, "camp", "cold");
    for (const auto& member : members) (void)core::appendRelocationCopy(plan, member);
    require(fs::is_regular_file(second / "cold" / "solutions" / "duffing" / "a.npz"), "layout kept");
    require(fs::is_regular_file(file), "relocation keeps sources");
    const auto routed = core::appendRelocationCopy(plan, second / "cold" / "solutions" / "duffing" / "a.npz");
    require(routed.alreadyAtDestination, "a routed completion is recorded without copying");
    require(fs::is_directory(second / "cold" / "time_grids"), "the shared grid travels with the campaign");
    fs::remove_all(root / "time_grids");  // the copied tree reads without the original grids
    require(core::readNpz(second / "cold" / "solutions" / "duffing" / "b.npz").contains("time"), "relocated copy reads");
    const auto planJson = json::load(plan.wstring());
    require(planJson.at("copies").array().size() == 5 && planJson.at("state").string() == "copying", "plan records copies");

    const auto catalog = json::parse(core::campaignCatalogJson(root));
    const auto& campaign = catalog.at("campaigns").array().front();
    require(campaign.at("id").string() == "camp" && campaign.at("label").string() == "sweep", "catalog campaign");
    require(campaign.at("solution_archives").array().size() == 2, "catalog finds solutions by metadata");
    require(campaign.at("missing_archives").array().size() == 1, "catalog reports missing archives");
    fs::remove_all(root);
    fs::remove_all(second);
}

}  // namespace

int main() {
    try {
        codecs();
        timeGrids();
        recompression();
        metadataEdits();
        auditCampaignAndReport();
        storageAndRelocation();
        std::cout << "dynsyslab_archive_maintenance_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_archive_maintenance_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
