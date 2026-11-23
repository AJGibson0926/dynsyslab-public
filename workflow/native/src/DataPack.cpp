#include <dynsyslab/workflow/DataPack.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/BuildInfo.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Version.hpp>
#include <dynsyslab/workflow/ArchiveDigest.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>

#include <algorithm>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace dynsyslab::workflow {
namespace {

constexpr const char* kReleaseUrl = "https://github.com/AJGibson0926/dynsyslab-public/releases/tag/";

std::string readText(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to read: " + path.u8string());
    std::ostringstream text; text << stream.rdbuf();
    return text.str();
}

void writeText(const std::filesystem::path& path, const std::string& text) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to write: " + path.u8string());
    stream << text;
}

// Every file of the pack other than the manifest and its own checksum, sorted.
std::vector<std::string> packFiles(const std::filesystem::path& folder) {
    std::vector<std::string> names;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(folder)) {
        if (!entry.is_regular_file()) continue;
        const auto name = std::filesystem::relative(entry.path(), folder).generic_u8string();
        if (name != "pack.json" && name != "pack.sha256") names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

json::Value fileTable(const std::filesystem::path& folder) {
    json::Value::Object table;
    for (const auto& name : packFiles(folder)) table[name] = json::Value(core::sha256File(folder / std::filesystem::u8path(name)));
    return json::Value(std::move(table));
}

void portabilityProblems(const json::Value& value, const std::string& file, json::Value::Array& problems) {
    if (const auto* text = std::get_if<std::string>(&value.storage())) {
        static const std::regex absolute(R"(^[A-Za-z]:[\\/]|^/[A-Za-z]|^\\\\)");
        static const std::regex stamped(R"(__20[0-9]{6}T[0-9]{6}Z)");
        if (std::regex_search(*text, absolute)) problems.emplace_back(file + ": absolute path: " + *text);
        if (std::regex_search(*text, stamped)) problems.emplace_back(file + ": wall-clock archive name: " + *text);
        if (text->find(".npz") != std::string::npos && text->find('*') != std::string::npos)
            problems.emplace_back(file + ": archive glob: " + *text);
    } else if (const auto* array = std::get_if<json::Value::Array>(&value.storage())) {
        for (const auto& item : *array) portabilityProblems(item, file, problems);
    } else if (const auto* object = std::get_if<json::Value::Object>(&value.storage())) {
        for (const auto& [key, item] : *object) { (void)key; portabilityProblems(item, file, problems); }
    }
}

}  // namespace

json::Value readDataPack(const std::filesystem::path& folder) {
    const auto path = folder / L"pack.json";
    if (!std::filesystem::is_regular_file(path)) throw std::invalid_argument("no pack.json in " + folder.u8string());
    auto manifest = json::parse(readText(path));
    if (manifest.at("schema").string().rfind("dynsyslab.data-pack/v", 0) != 0)
        throw std::invalid_argument("not a DynSysLab data pack: " + path.u8string());
    return manifest;
}

std::filesystem::path dataPackRootJob(const std::filesystem::path& folder) {
    const auto manifest = readDataPack(folder);
    return folder / std::filesystem::u8path(manifest.at("root_job").string("jobs/00-reproduce-all.job.json"));
}

json::Value inspectDataPack(const std::filesystem::path& folder) {
    const auto manifest = readDataPack(folder);
    const bool understood = manifest.at("schema").string() == kDataPackSchema;
    const auto pinned = manifest.at("release").at("version").string();
    return json::Value(json::Value::Object{
        {"schema", manifest.at("schema")}, {"schema_supported", json::Value(understood)},
        {"pack_id", manifest.at("pack_id")}, {"title", manifest.at("title")}, {"date", manifest.at("date")},
        {"release", manifest.at("release")}, {"build", manifest.at("build")}, {"sealed", json::Value(manifest.at("sealed").boolean(false))},
        {"running_version", json::Value(std::string(core::version))},
        {"runnable_here", json::Value(understood && pinned == core::version)},
        {"root_job", manifest.at("root_job")}, {"figures", manifest.at("figures")},
        {"file_count", json::Value(double(manifest.at("files").object().size()))}});
}

json::Value checkDataPackPin(const json::Value& job) {
    const auto reference = job.at("data_pack").string();
    if (reference.empty()) return json::Value();
    const auto file = job.at("job_file").string();
    if (file.empty()) throw std::invalid_argument("a data-pack job must be run from its file");
    const auto folder = std::filesystem::weakly_canonical(
        std::filesystem::u8path(file).parent_path() / std::filesystem::u8path(reference)).parent_path();
    const auto manifest = readDataPack(folder);
    if (manifest.at("schema").string() != kDataPackSchema)
        throw std::runtime_error("this release cannot run pack schema " + manifest.at("schema").string()
                                 + "; 'datapack inspect' still reads it");
    const auto pinned = manifest.at("release").at("version").string();
    const auto tag = manifest.at("release").at("tag").string("v" + pinned);
    const bool match = pinned == core::version;
    const bool allowed = job.at("runtime").at("allow_version_mismatch").boolean(false);
    if (!match && !allowed)
        throw std::runtime_error("data pack '" + manifest.at("pack_id").string() + "' is pinned to DynSysLab " + pinned
            + " but this is " + std::string(core::version) + ". Install the pinned release: " + kReleaseUrl + tag
            + " (--allow-version-mismatch runs it anyway, labelled non-reproducing)");
    return json::Value(json::Value::Object{
        {"pack_id", manifest.at("pack_id")}, {"pinned_version", json::Value(pinned)},
        {"running_version", json::Value(std::string(core::version))}, {"reproducing", json::Value(match)},
        {"label", json::Value(match ? "pinned-release" : "NON-REPRODUCING: version mismatch")},
        {"build", json::parse(core::buildIdentityJson())}});
}

void refreshDataPackHashes(const std::filesystem::path& folder) {
    auto manifest = readDataPack(folder);
    manifest.at("files") = fileTable(folder);
    const auto text = json::dump(manifest) + "\n";
    writeText(folder / L"pack.json", text);
    writeText(folder / L"pack.sha256", core::sha256Text(text) + "  pack.json\n");
}

json::Value validateDataPack(const std::filesystem::path& folder) {
    json::Value::Array problems;
    const auto manifest = readDataPack(folder);
    if (manifest.at("schema").string() != kDataPackSchema) problems.emplace_back("unsupported schema: " + manifest.at("schema").string());
    for (const auto* field : {"pack_id", "title", "date"})
        if (manifest.at(field).string().empty()) problems.emplace_back(std::string("manifest lacks ") + field);
    if (manifest.at("release").at("version").string().empty()) problems.emplace_back("manifest lacks release.version");
    const auto sidecar = folder / L"pack.sha256";
    if (!std::filesystem::is_regular_file(sidecar)) problems.emplace_back("pack.sha256 is missing");
    else if (readText(sidecar).substr(0, 64) != core::sha256File(folder / L"pack.json"))
        problems.emplace_back("pack.json does not match pack.sha256");
    const auto actual = fileTable(folder);
    for (const auto& [name, hash] : manifest.at("files").object()) {
        if (actual.at(name).isNull()) problems.emplace_back("missing file: " + name);
        else if (actual.at(name).string() != hash.string()) problems.emplace_back("file changed: " + name);
    }
    for (const auto& [name, hash] : actual.object()) {
        (void)hash;
        if (manifest.at("files").at(name).isNull()) problems.emplace_back("file not in the manifest: " + name);
    }
    double jobCount = 0.0;
    try {
        for (const auto& name : packFiles(folder))
            if (name.size() > 5 && name.compare(name.size() - 5, 5, ".json") == 0)
                portabilityProblems(json::parse(readText(folder / std::filesystem::u8path(name))), name, problems);
        const auto file = dataPackRootJob(folder);
        const auto root = bindJobFile(json::parse(readText(file)), file);
        if (root.at("data_pack").string().empty()) problems.emplace_back("the root job does not name its pack");
        jobCount = double(loadJobGroup(root).leaves.size());
    } catch (const std::exception& error) { problems.emplace_back(std::string("job graph: ") + error.what()); }
    const bool valid = problems.empty();
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.data-pack-validation/v1")},
        {"pack_id", manifest.at("pack_id")}, {"valid", json::Value(valid)}, {"jobs", json::Value(jobCount)},
        {"problems", json::Value(std::move(problems))}});
}

json::Value verifyDataPackRun(const std::filesystem::path& folder, const std::filesystem::path& ledgerPath) {
    const auto manifest = readDataPack(folder);
    const auto ledger = json::parse(readText(ledgerPath));
    json::Value::Object jobs;
    double exact = 0.0, equivalent = 0.0, failed = 0.0;
    const auto count = [&](const std::string& status) {
        if (status == "exact") exact += 1.0;
        else if (status == "numerically_equivalent") equivalent += 1.0;
        else failed += 1.0;
    };
    const auto expectedPath = folder / L"expected.json";
    const bool sealed = std::filesystem::is_regular_file(expectedPath);
    if (sealed) {
        const auto expected = json::parse(readText(expectedPath));
        for (const auto& [id, want] : expected.at("jobs").object()) {
            const auto& got = ledger.at("jobs").at(id);
            json::Value::Array results;
            std::string worst = "exact";
            if (got.at("status").string() != "complete" || got.at("digests").array().size() != want.array().size()) worst = "failed";
            else for (std::size_t k = 0; k < want.array().size(); ++k) {
                auto compared = compareDigests(got.at("digests").array()[k], want.array()[k]);
                const auto status = compared.at("status").string();
                if (status == "failed" || (status == "numerically_equivalent" && worst == "exact")) worst = status;
                if (status != "exact") results.push_back(std::move(compared));
            }
            count(worst);
            jobs[id] = json::Value(json::Value::Object{{"status", json::Value(worst)}, {"differences", json::Value(std::move(results))}});
        }
    }
    const bool reproducing = ledger.at("data_pack").at("reproducing").boolean(false);
    std::string verdict = !sealed ? "unsealed" : !reproducing ? "non_reproducing"
        : failed > 0.0 || ledger.at("status").string() != "complete" ? "failed"
        : equivalent > 0.0 ? "numerically_equivalent" : "verified_exact";
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.data-pack-verification/v1")},
        {"pack_id", manifest.at("pack_id")}, {"verdict", json::Value(verdict)}, {"reproducing", json::Value(reproducing)},
        {"run_status", ledger.at("status")}, {"exact", json::Value(exact)}, {"numerically_equivalent", json::Value(equivalent)},
        {"failed", json::Value(failed)}, {"jobs", json::Value(std::move(jobs))}});
}

json::Value sealDataPack(const std::filesystem::path& folder, const std::filesystem::path& ledgerPath,
                         const std::filesystem::path& output) {
    const auto ledger = json::parse(readText(ledgerPath));
    if (ledger.at("status").string() != "complete") throw std::runtime_error("seal needs a complete run");
    if (!ledger.at("data_pack").at("reproducing").boolean(false)) throw std::runtime_error("seal refuses a non-reproducing run");
    const auto identity = core::buildIdentity();
    if (identity.sourceModified || identity.sourceCommit == "unknown")
        throw std::runtime_error("seal needs a release build of a clean, committed source tree");
    auto manifest = readDataPack(folder);
    if (std::filesystem::weakly_canonical(folder) == std::filesystem::weakly_canonical(output))
        throw std::invalid_argument("seal writes a separate copy; choose another --output");
    std::filesystem::create_directories(output);
    for (const auto& name : packFiles(folder)) {
        if (name == "expected.json") continue;
        const auto target = output / std::filesystem::u8path(name);
        std::filesystem::create_directories(target.parent_path());
        std::filesystem::copy_file(folder / std::filesystem::u8path(name), target, std::filesystem::copy_options::overwrite_existing);
    }
    json::Value::Object jobs;
    for (const auto& [id, entry] : ledger.at("jobs").object()) jobs[id] = entry.at("digests");
    writeText(output / L"expected.json", json::dump(json::Value(json::Value::Object{
        {"schema", json::Value("dynsyslab.data-pack-expected/v1")}, {"jobs", json::Value(std::move(jobs))}})) + "\n");
    manifest.at("sealed") = json::Value(true);
    manifest.at("release").at("source_commit") = json::Value(identity.sourceCommit);
    manifest.at("build") = json::parse(core::buildIdentityJson());
    writeText(output / L"pack.json", json::dump(manifest) + "\n");
    refreshDataPackHashes(output);
    return validateDataPack(output);
}

json::Value exportDataPack(const std::filesystem::path& folder, const std::filesystem::path& directory) {
    const auto report = validateDataPack(folder);
    if (!report.at("valid").boolean()) throw std::runtime_error("export refuses an invalid pack: " + json::dump(report.at("problems")));
    const auto manifest = readDataPack(folder);
    const auto id = manifest.at("pack_id").string();
    const auto date = manifest.at("date").string();
    auto name = manifest.at("name").string(id);
    const auto suffix = "-" + date;
    if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
        name.resize(name.size() - suffix.size());
    const auto file = "DynSysLab-datapack-" + name + "-" + date + "-v" + manifest.at("release").at("version").string() + ".zip";
    std::vector<std::pair<std::string, std::string>> entries;
    auto names = packFiles(folder);
    names.push_back("pack.json"); names.push_back("pack.sha256");
    std::sort(names.begin(), names.end());
    for (const auto& item : names) entries.emplace_back(id + "/" + item, readText(folder / std::filesystem::u8path(item)));
    const auto path = directory / std::filesystem::u8path(file);
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    core::writeTextZip(path, entries);
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.data-pack-export/v1")},
        {"file", json::Value(path.generic_u8string())}, {"sha256", json::Value(core::sha256File(path))},
        {"sealed", json::Value(manifest.at("sealed").boolean(false))}});
}

}  // namespace dynsyslab::workflow
