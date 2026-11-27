// Composite jobs, symbolic links, cleanup, resume, version pin, digests and the generic K-ROM and
// figure-table stages, on a tiny self-contained pack.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/PortableNormal.hpp>
#include <dynsyslab/core/Version.hpp>
#include <dynsyslab/workflow/ArchiveDigest.hpp>
#include <dynsyslab/workflow/DataPack.hpp>
#include <dynsyslab/workflow/JobGroup.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

using namespace dynsyslab;
namespace fs = std::filesystem;
namespace json = core::json;

namespace {

void check(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

void write(const fs::path& path, const std::string& text) {
    fs::create_directories(path.parent_path());
    std::ofstream(path, std::ios::binary) << text;
}

std::string read(const fs::path& path) {
    std::ifstream stream(path, std::ios::binary); std::ostringstream out; out << stream.rdbuf(); return out.str();
}

template <class F> bool throwsWith(F&& action, const std::string& text) {
    try { action(); } catch (const std::exception& error) { return std::string(error.what()).find(text) != std::string::npos; }
    return false;
}

const char* kSolve = R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"tiny-solve",
 "outputs":{"trajectory":"archive"},
 "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
 "time":{"start":0,"stop":"${params.stop}","output_step":0.01,"coordinate_system":"automatic"},
 "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,"values":[[-1.0,0.1]],"ranges":[]},
 "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[-1.0,0.0]}},
 "solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},"archive":{"enabled":true,"data_role":"free_training"}})";

const char* kTable = R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"tiny-table",
 "sources":["${inputs.run}"],"methods":[{"name":"figure_dataset","settings":{"inputs":{"run":"${inputs.run}"},
 "arrays":[{"name":"time","from":"run","array":"time"},{"name":"x","from":"run","array":"values","column":0,"scale":2.0}]}}],
 "archive":{"enabled":true}})";

const char* kFigure = R"({"schema":"dynsyslab.compute-job/v1","kind":"viewer_figure","label":"tiny-figure",
 "source":"${inputs.table}","output":"tiny.png",
 "figure":{"width":320,"height":240,"title":"tiny","axes":["x"],"style":"${params.style}",
           "curves":[{"series":"x","label":"x","color":"#202020","line_width":1.0}]}})";

const char* kRoot = R"({"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"tiny-pack","data_pack":"../pack.json",
 "cleanup":{"enabled":true,"min_free_gb":0},
 "parameters":{"stop":0.5,"style":{"font_scale":1.2,"line_width_scale":2.0,"default_line_width":2.0}},
 "jobs":[{"id":"figure","job":"20/figure.job.json","inputs":{"table":"${table.archive}"}},
         {"id":"table","job":"10/table.job.json","inputs":{"run":"${solve.outputs.trajectory}"}},
         {"id":"solve","job":"10/solve.job.json"}]})";

json::Value loadRoot(const fs::path& pack, json::Value::Object runtime = {}) {
    const auto file = pack / "jobs/00-reproduce-all.job.json";
    auto job = workflow::bindJobFile(json::parse(read(file)), file);
    if (!runtime.empty()) job.at("runtime") = json::Value(std::move(runtime));
    return job;
}

void writePack(const fs::path& pack, const std::string& version) {
    write(pack / "jobs/00-reproduce-all.job.json", kRoot);
    write(pack / "jobs/10/solve.job.json", kSolve);
    write(pack / "jobs/10/table.job.json", kTable);
    write(pack / "jobs/20/figure.job.json", kFigure);
    write(pack / "pack.json", json::dump(json::Value(json::Value::Object{
        {"schema", json::Value(workflow::kDataPackSchema)}, {"pack_id", json::Value("tiny-pack")}, {"title", json::Value("tiny")},
        {"date", json::Value("2026-11-23")}, {"root_job", json::Value("jobs/00-reproduce-all.job.json")},
        {"release", json::Value(json::Value::Object{{"version", json::Value(version)}, {"tag", json::Value("v" + version)}})}})));
    workflow::refreshDataPackHashes(pack);
}

void testPortableNormal() {
    std::mt19937_64 generator(20261104);
    // Reference values: the same words through Python's statistics.NormalDist (also AS 241).
    for (const double expected : {0.15130098367977934, 1.4357790859171153, 0.42559972185461376, -0.9956739334924842})
        check(std::abs(core::portableNormal(generator) - expected) < 1e-14, "normal stream changed");
    check(std::abs(core::portableNormalQuantile(0.975) - 1.959963984540054) < 1e-14, "central quantile");
    check(std::abs(core::portableNormalQuantile(1e-10) + 6.361340902404056) < 1e-12, "tail quantile");
}

void testGroupRun(const fs::path& base) {
    const auto pack = base / "pack";
    writePack(pack, std::string(core::version));
    check(workflow::validateDataPack(pack).at("valid").boolean(), "pack should validate");
    const auto plan = workflow::planJob(loadRoot(pack));
    check(plan.at("job_count").number() == 3.0, "three leaf jobs");
    check(plan.at("jobs").array()[0].at("id").string() == "solve", "dependency order");

    const auto root = base / "run";
    const auto first = workflow::executeJob(loadRoot(pack), root);
    const auto ledgerPath = first.manifest;
    auto ledger = json::parse(read(ledgerPath));
    check(ledger.at("status").string() == "complete", "ledger complete");
    check(fs::is_regular_file(ledger.at("jobs").at("figure").at("archives").array()[0].string()), "figure written");
    for (const auto* id : {"solve", "table"}) {
        check(ledger.at("jobs").at(id).at("deleted").boolean(), std::string("intermediate kept: ") + id);
        check(!fs::exists(ledger.at("jobs").at(id).at("folder").string()), "intermediate folder remains");
    }
    check(ledger.at("jobs").at("table").at("inputs").object().size() == 1, "resolution recorded");

    (void)workflow::executeJob(loadRoot(pack), root);  // resume: nothing recomputes
    ledger = json::parse(read(ledgerPath));
    for (const auto* id : {"solve", "table", "figure"}) check(ledger.at("jobs").at(id).at("reused").boolean(), "resume reused");

    // Verification categories against sealed digests.
    json::Value::Object digests;
    for (const auto& [id, entry] : ledger.at("jobs").object()) digests[id] = entry.at("digests");
    const auto sealed = base / "sealed";
    fs::copy(pack, sealed, fs::copy_options::recursive);
    write(sealed / "expected.json", json::dump(json::Value(json::Value::Object{{"jobs", json::Value(digests)}})));
    workflow::refreshDataPackHashes(sealed);
    check(workflow::verifyDataPackRun(sealed, ledgerPath).at("verdict").string() == "verified_exact", "exact verdict");
    digests["figure"].array()[0].at("sha256") = json::Value(std::string(64, '0'));
    write(sealed / "expected.json", json::dump(json::Value(json::Value::Object{{"jobs", json::Value(digests)}})));
    check(workflow::verifyDataPackRun(sealed, ledgerPath).at("verdict").string() == "failed", "changed figure must fail");

    // Tampering with any file is detected.
    write(pack / "jobs/10/solve.job.json", std::string(kSolve) + " ");
    check(!workflow::validateDataPack(pack).at("valid").boolean(), "tamper detection");
    writePack(pack, std::string(core::version));
    write(pack / "pack.json", read(pack / "pack.json") + " ");
    check(!workflow::validateDataPack(pack).at("valid").boolean(), "manifest tamper detection");
    writePack(pack, std::string(core::version));

    check(throwsWith([&] { (void)workflow::executeJob(loadRoot(pack, {{"min_free_gb", json::Value(1e9)}}), base / "floor"); },
                     "disk preflight"), "disk floor refuses");
}

void testPin(const fs::path& base) {
    const auto pack = base / "old";
    writePack(pack, "0.0.1");
    check(throwsWith([&] { (void)workflow::executeJob(loadRoot(pack), base / "old-run"); }, "releases/tag/v0.0.1"), "version refusal");
    check(workflow::inspectDataPack(pack).at("pack_id").string() == "tiny-pack", "old pack stays inspectable");
    const auto run = workflow::executeJob(loadRoot(pack, {{"allow_version_mismatch", json::Value(true)}}), base / "old-run");
    const auto ledger = json::parse(read(run.manifest));
    check(!ledger.at("data_pack").at("reproducing").boolean(), "mismatch labelled");
    json::Value::Object digests;
    for (const auto& [id, entry] : ledger.at("jobs").object()) digests[id] = entry.at("digests");
    write(pack / "expected.json", json::dump(json::Value(json::Value::Object{{"jobs", json::Value(digests)}})));
    check(workflow::verifyDataPackRun(pack, run.manifest).at("verdict").string() == "non_reproducing", "never verified");
}

void testGraphErrors(const fs::path& base) {
    const auto dir = base / "graphs";
    write(dir / "a.job.json", R"({"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"a","jobs":[{"id":"b","job":"b.job.json"}]})");
    write(dir / "b.job.json", R"({"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"b","jobs":[{"id":"a","job":"a.job.json"}]})");
    auto job = workflow::bindJobFile(json::parse(read(dir / "a.job.json")), dir / "a.job.json");
    check(throwsWith([&] { (void)workflow::loadJobGroup(job); }, "includes itself"), "file cycle");
    write(dir / "c.job.json", R"({"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"c",
      "jobs":[{"id":"x","job":{"kind":"solve"},"depends_on":["y"]},{"id":"y","job":{"kind":"solve"},"depends_on":["x"]}]})");
    job = workflow::bindJobFile(json::parse(read(dir / "c.job.json")), dir / "c.job.json");
    check(throwsWith([&] { (void)workflow::loadJobGroup(job); }, "dependency cycle"), "dependency cycle");
    std::string nested = R"({"kind":"solve"})";
    for (int k = 0; k < 20; ++k)
        nested = R"({"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"n","jobs":[{"id":"n","job":)" + nested + "}]}";
    write(dir / "deep.job.json", nested);
    job = workflow::bindJobFile(json::parse(nested), dir / "deep.job.json");
    check(throwsWith([&] { (void)workflow::loadJobGroup(job); }, "nesting deeper"), "depth limit");
}

void testComparison(const fs::path& base) {
    core::NpzRecord a;
    a.arrays["x"] = core::NpyArray::fromDoubles({3}, std::vector<double>{1.0, 2.0, 3.0});
    a.arrays["metadata_json"] = core::NpyArray::fromString(R"({"created_utc":"1"})");
    auto b = a, c = a;
    b.arrays["x"] = core::NpyArray::fromDoubles({3}, std::vector<double>{1.0, 2.0 + 4e-12, 3.0});
    c.arrays["x"] = core::NpyArray::fromDoubles({3}, std::vector<double>{1.0, 2.5, 3.0});
    c.arrays["metadata_json"] = core::NpyArray::fromString(R"({"created_utc":"2"})");
    const auto pa = base / "a.npz", pb = base / "b.npz", pc = base / "c.npz";
    core::writeNpzAtomic(pa, a); core::writeNpzAtomic(pb, b); core::writeNpzAtomic(pc, c);
    check(workflow::compareResults(pa, pa).at("status").string() == "exact", "exact");
    const auto near = workflow::compareResults(pb, pa);
    check(near.at("status").string() == "numerically_equivalent", "equivalent");
    check(near.at("differing_arrays").at("x").at("flat_index").number() == 1.0, "location reported");
    check(workflow::compareResults(pc, pa).at("status").string() == "failed", "failed");
    check(workflow::resultDigest(pa).at("sha256").string() != workflow::resultDigest(pc).at("sha256").string(), "digest differs");
}

void testKromEditing(const fs::path& base) {
    using C = std::complex<double>;
    const std::vector<C> coefficients{C(1, 1), C(2, 0), C(0, 1), C(1, -1), C(2, 0), C(0, -1)};
    core::PolynomialDictionary monomials{2, {1, 0, 0, 1, 2, 0}, {0.5, 0.0}, {2.0, 1.0}};
    core::FunctionDomain plane;
    plane.dimension = 2;
    core::KoopmanModel krom(core::FunctionExpansion(std::move(monomials), coefficients, plane),
                            core::Spectrum::fromContinuousRates({C(-1, 2), C(-1, -2)}));
    core::ReconstructionResult readout{{C(1, 0), C(1, 0)}, 1, 2, 0.0, 0.0, 0.0};
    readout.offset = {0.25};
    krom.setReconstruction("y", std::move(readout));
    const auto source = base / "krom.npz";
    core::saveKoopmanModel(source, krom);
    const auto run = [&](const char* method, const std::string& settings) {
        return workflow::executeJob(json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"edit","sources":[")"
            + source.generic_u8string() + R"("],"methods":[{"name":")" + method + R"(","settings":)" + settings + "}]}"), base / "edit").archives.at(0);
    };
    const auto psi = [](const std::vector<C>& k, double cx, double x, double y, std::size_t m) {
        const double z0 = (x - cx) / 2.0, z1 = y;
        return k[m * 3] * z0 + k[m * 3 + 1] * z1 + k[m * 3 + 2] * (z0 * z0);
    };
    const auto movedPath = run("krom_recenter", R"({"target_state":[1.0,0.5]})");
    const auto moved = core::readNpz(movedPath);
    const auto& before = coefficients;
    const auto after = core::loadKoopmanModel(movedPath).expansion().coefficients();
    for (std::size_t m = 0; m < 2; ++m) {
        const auto shift = psi(before, 0.5, 1.0, 0.5, m);
        const auto expected = psi(before, 0.5, 1.7, -0.3, m) - shift;
        check(std::abs(psi(after, 1.0, 1.7, -0.3 - 0.5, m) - expected) < 1e-12, "recentred eigenfunction");
    }
    const double intercept = 0.25 + (psi(before, 0.5, 1.0, 0.5, 0) + psi(before, 0.5, 1.0, 0.5, 1)).real();
    check(std::abs(moved.at("reconstruction_offset__y").doubles()[0] - intercept) < 1e-12, "readout intercept");
    const auto split = core::readNpz(run("krom_readout", R"({"readouts":[{"name":"lifted_split","kind":"lifted_split"}]})"));
    const auto rows = split.at("reconstruction__lifted_split").complexDoubles();
    check(split.at("reconstruction__lifted_split").shape == std::vector<std::size_t>({4, 2}), "split shape");
    check(rows[0] == C(0.5, 0) && rows[1] == C(0.5, 0) && rows[2] == C(0, -0.5) && rows[3] == C(0, 0.5), "split pairing");
}

}  // namespace

int main() {
    const auto base = fs::temp_directory_path() / ("dynsyslab-job-group-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    int status = 0;
    try {
        fs::create_directories(base);
        testPortableNormal();
        testComparison(base);
        testKromEditing(base);
        testGraphErrors(base);
        testGroupRun(base);
        testPin(base);
        std::cout << "job group tests passed\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
    std::error_code ignored;
    fs::remove_all(base, ignored);
    return status;
}
