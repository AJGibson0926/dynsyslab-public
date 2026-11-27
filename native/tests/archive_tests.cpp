#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveAudit.hpp>
#include <dynsyslab/core/ArchiveDeduplication.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void roundTrip() {
    const auto path = std::filesystem::temp_directory_path() / "dynsyslab-native-archive-test.npz";
    dynsyslab::core::NpzRecord source;
    source.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({3}, {0.0, 0.5, 1.0});
    source.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles({3, 2},
        {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    source.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"schema":"dynsyslab.solution/v1","system":"duffing"})");
    dynsyslab::core::writeNpzAtomic(path, source);
    const auto index=dynsyslab::core::readNpzIndex(path);
    require(index.contains("time") && index.contains("values")
        && index.contains("metadata_json"),"lightweight NPZ index names");
    require(index.metadataJson.find("duffing")!=std::string::npos,
        "lightweight NPZ metadata");
    const auto loaded = dynsyslab::core::readNpz(path);
    require(loaded.at("time").shape == std::vector<std::size_t>{3}, "time shape");
    require(loaded.at("values").shape == std::vector<std::size_t>({3, 2}), "value shape");
    require(loaded.at("values").doubles()[5] == 6.0, "numeric round trip");
    require(loaded.metadataJson().find("duffing") != std::string::npos, "metadata round trip");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}

void longPublicationPath() {
    const auto root=std::filesystem::temp_directory_path()/"dynsyslab-long-publication";
    std::filesystem::create_directories(root);
    const auto name=std::string(242-root.wstring().size()-5,'x')+".npz";
    const auto file=root/name;
    dynsyslab::core::NpzRecord source;
    source.arrays["value"]=dynsyslab::core::NpyArray::fromDoubles({1},{42});
    dynsyslab::core::writeNpzAtomic(file,source);
    require(dynsyslab::core::readNpz(file).at("value").doubles()[0]==42,
        "long Unicode-path atomic archive round trip");
    std::filesystem::remove(file);std::filesystem::remove(root);
}

void readExisting(const char* argument) {
    if (!argument) return;
    const auto record = dynsyslab::core::readNpz(std::filesystem::u8path(argument));
    require(record.contains("metadata_json"), "existing archive metadata");
    require(record.contains("values") || record.contains("state"), "existing archive values");
}

dynsyslab::core::NpzRecord solution(const std::vector<double>& values) {
    dynsyslab::core::NpzRecord result;
    result.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({2}, {0.0, 1.0});
    result.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles({2, 2}, values);
    result.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"schema":"dynsyslab.solution/v1","system":"duffing","coordinate_system":"dimensionless","time_span":[0,1]})");
    return result;
}

void auditDuplicates() {
    const auto root = std::filesystem::temp_directory_path()
        / "dynsyslab-native-archive-audit-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    const auto directory = root / "solutions" / "duffing";
    std::filesystem::create_directories(directory);
    dynsyslab::core::writeNpzAtomic(directory / "one.npz", solution({1.0, 0.0, 0.5, -0.5}));
    dynsyslab::core::writeNpzAtomic(directory / "two.npz", solution({1.0, 0.0, 0.5, -0.5}));
    dynsyslab::core::writeNpzAtomic(directory / "three.npz", solution({1.0, 0.0, 0.25, -0.25}));
    dynsyslab::core::ArchiveAuditScope scope; scope.modifiedWithinDays = -1.0;
    const auto report = dynsyslab::core::auditArchives(root, scope);
    require(report.scannedFiles == 3, "audit file count");
    require(report.exactDuplicateGroups.size() == 1, "exact duplicate group");
    require(report.exactDuplicateGroups.front().files.size() == 2, "exact duplicate members");
    require(report.samePlanGroups.size() == 1, "same plan group");
    require(report.samePlanGroups.front().files.size() == 3, "same plan members");
    require(dynsyslab::core::archiveAuditJson(report, scope).find("read_only") != std::string::npos,
            "audit JSON contract");
    scope.sections = {".."};
    bool rejected = false;
    try { (void)dynsyslab::core::auditArchives(root, scope); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "audit section traversal rejection");
    std::filesystem::remove_all(root, ignored);
}

void deduplicateRecoverably() {
    const auto root = std::filesystem::temp_directory_path()
        / "dynsyslab-native-archive-deduplication-test";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    const auto directory = root / "solutions" / "duffing";
    std::filesystem::create_directories(directory);
    const auto one = directory / "one.npz";
    const auto two = directory / "two.npz";
    const auto three = directory / "three.npz";
    const auto values = solution({1.0, 0.0, 0.5, -0.5});
    dynsyslab::core::writeNpzAtomic(one, values);
    dynsyslab::core::writeNpzAtomic(two, values);
    dynsyslab::core::writeNpzAtomic(three, values);
    const auto now = std::filesystem::file_time_type::clock::now();
    std::filesystem::last_write_time(one, now - std::chrono::seconds(3));
    std::filesystem::last_write_time(two, now - std::chrono::seconds(2));
    std::filesystem::last_write_time(three, now - std::chrono::seconds(1));
    std::filesystem::create_directories(root / "hierarchies");
    std::ofstream(root / "hierarchies" / "reference.json")
        << "{\"solution\":\"" << two.generic_u8string() << "\"}";
    dynsyslab::core::ArchiveAuditScope scope; scope.modifiedWithinDays = -1.0;
    const auto dry = dynsyslab::core::deduplicateArchives(root, scope);
    require(dry.dryRun, "deduplication defaults to dry-run");
    require(dry.canonicalFiles == std::vector<std::filesystem::path>{three},
            "deduplication canonical selection");
    require(dry.candidates == std::vector<std::filesystem::path>{one},
            "deduplication unreferenced candidate");
    require(dry.skipped.size() == 1 && dry.skipped.front().candidate == two,
            "deduplication referenced skip");
    require(std::filesystem::exists(one), "dry-run leaves candidate intact");
    const auto applied = dynsyslab::core::deduplicateArchives(root, scope, true);
    require(!applied.dryRun && applied.movedFiles.size() == 1,
            "deduplication apply result");
    require(!std::filesystem::exists(one) && std::filesystem::exists(two)
            && std::filesystem::exists(three), "only unreferenced duplicate moved");
    require(std::filesystem::exists(applied.movedFiles.front()),
            "duplicate retained in recovery trash");
    require(std::filesystem::exists(applied.trashRoot / "deduplication.json"),
            "deduplication recovery manifest");
    require(dynsyslab::core::archiveDeduplicationJson(applied).find("recoverable")
            != std::string::npos, "deduplication JSON contract");
    std::filesystem::remove_all(root, ignored);
}

void createOnlyPublication() {
    const auto root=std::filesystem::temp_directory_path()/("dynsyslab-create-only-"
        +std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    const auto path=root/"result.npz";
    const auto original=solution({1.,2.,3.,4.});
    dynsyslab::core::writeNpzAtomicNew(path,original);
    bool rejected=false;
    try { dynsyslab::core::writeNpzAtomicNew(path,solution({5.,6.,7.,8.})); }
    catch(const std::exception&) { rejected=true; }
    require(rejected,"create-only publication replaced an existing result");
    require(dynsyslab::core::readNpz(path).at("values").doubles()
        ==original.at("values").doubles(),"create-only publication changed source bytes");
    for(const auto& entry:std::filesystem::directory_iterator(root))
        require(entry.path()==path,"failed publication left a partial archive");
    std::filesystem::remove_all(root);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        roundTrip();
        longPublicationPath();
        createOnlyPublication();
        auditDuplicates();
        deduplicateRecoverably();
        if (argc > 1) readExisting(argv[1]);
        std::cout << "dynsyslab_archive_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_archive_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
