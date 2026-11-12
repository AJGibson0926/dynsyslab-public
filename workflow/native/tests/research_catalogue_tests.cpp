#include "../src/MetadataPaths.hpp"
#include <chrono>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <filesystem>
#include <iostream>
#include <stdexcept>
namespace {
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F action, const char *message) {
    bool rejected = false;
    try {
        action();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, message);
}
} // namespace
int main() {
    namespace j = dynsyslab::core::json;
    using namespace dynsyslab::workflow;
    const auto root =
        std::filesystem::temp_directory_path() /
        ("dynsyslab-research-test-" +
         std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    try {
        for (const auto &descriptor : dynsyslab::core::availableSystems()) {
            auto catalogue =
                j::parse(R"({"schema":"dynsyslab.research-query/v1","action":"list","system":{}})");
            catalogue.at("system").object()["name"] = descriptor.identifier;
            if (evaluateResearchCatalogueQuery(catalogue, root).at("papers").array().empty())
                throw std::runtime_error("Native system is missing its researched bibliography: " +
                                         descriptor.identifier);
        }
        auto query = j::parse(
            R"({"schema":"dynsyslab.research-query/v1","action":"list","system":{"name":"rpe"}})");
        auto listed = evaluateResearchCatalogueQuery(query, root);
        const auto count = listed.at("papers").array().size();
        require(count >= 6, "Built-in RPE bibliography absent");
        require(!std::filesystem::exists(root), "Read-only catalogue created files");
        bool chang = false;
        for (const auto &paper : listed.at("papers").array())
            if (paper.at("doi").string() == "10.1063/1.865830")
                chang = paper.at("year").number() == 1986;
        require(chang, "Chang-Chen reference missing or incorrectly dated");
        query.object()["action"] = "add";
        query.object()["paper"] =
            j::parse(R"({"id":"test-extension","title":"Local reference fixture","year":2026,
            "doi":"https://doi.org/10.1234/TEST", "authors":["Test Author"], "notes":"fixture only"})");
        auto added = evaluateResearchCatalogueQuery(query, root);
        require(added.at("status").string() == "added", "Paper was not added");
        const auto file = std::filesystem::u8path(added.at("added_path").string());
        const auto before = dynsyslab::core::sha256File(detail::metadataIoPath(file));
        auto repeated = evaluateResearchCatalogueQuery(query, root);
        require(repeated.at("status").string() == "already_present",
                "Repeated addition was not idempotent");
        require(dynsyslab::core::sha256File(detail::metadataIoPath(file)) == before,
                "Repeated add changed immutable evidence");
        query.at("paper").object()["id"] = "alternate-id";
        require(evaluateResearchCatalogueQuery(query, root).at("papers").array().size() ==
                    count + 1,
                "Same DOI under another id duplicated a paper");
        query.at("paper").object()["id"] = "test-extension";
        query.at("paper").object()["doi"] = "10.1234/different";
        rejects([&] { evaluateResearchCatalogueQuery(query, root); },
                "Conflicting id was accepted");
        query.at("paper").object()["id"] = "../invalid";
        rejects([&] { evaluateResearchCatalogueQuery(query, root); }, "Unsafe paper id accepted");
        dynsyslab::core::SystemConfiguration config;
        config.identifier = "rpe";
        config.parameters["mu_pa_s"] = 0;
        auto system = dynsyslab::core::createSystem(config);
        require(loadSystemKnowledge(*system, root).papers.size() == count + 1,
                "Bibliography did not survive a new system instance and parameter change");
        auto phase = j::parse(
            R"({"schema":"dynsyslab.phase-space-query/v1","mode":"lookup", "system":{"name":"duffing"},
            "options":{"lower":[-2,-2],"upper":[2,2],"manifolds":false,"periodic_orbits":false,"vector_field":false}})");
        auto miss = evaluatePhaseSpaceQuery(phase, root);
        require(miss.at("cache").at("status").string() == "miss",
                "Fresh phase query was not a miss");
        phase.object()["mode"] = "compute";
        auto computed = evaluatePhaseSpaceQuery(phase, root);
        require(computed.at("cache").at("status").string() == "computed",
                "Phase result not cached");
        require(!computed.at("knowledge_checks").array().empty(), "System guides were not checked");
        phase.object()["mode"] = "lookup";
        require(evaluatePhaseSpaceQuery(phase, root).at("cache").at("status").string() == "hit",
                "Cache not reused");
        phase.at("options").object()["lower"] = j::parse("[2,2]");
        rejects([&] { evaluatePhaseSpaceQuery(phase, root); },
                "Invalid box accepted during cache lookup");
        for (const auto *id : {"heat_equation_1d", "mackey_glass", "ising", "ising_3d", "potts"}) {
            auto familyQuery =
                j::parse(R"({"schema":"dynsyslab.research-query/v1","action":"list","system":{}})");
            familyQuery.at("system").object()["name"] = id;
            require(evaluateResearchCatalogueQuery(familyQuery, root).at("system").string() == id,
                    "Non-ODE system has no bibliography interface");
        }
        auto pdeQuery = j::parse(R"({"schema":"dynsyslab.research-query/v1","action":"add",
            "system":{"name":"heat_equation_1d"},"paper":{"id":"pde-fixture","title":"PDE catalogue fixture"}})");
        pdeQuery.object()["action"] = "list";
        const auto pdeBuiltinCount =
            evaluateResearchCatalogueQuery(pdeQuery, root).at("papers").array().size();
        pdeQuery.object()["action"] = "add";
        require(evaluateResearchCatalogueQuery(pdeQuery, root).at("status").string() == "added",
                "PDE catalogue could not retain an addition");
        pdeQuery.object()["action"] = "list";
        require(evaluateResearchCatalogueQuery(pdeQuery, root).at("papers").array().size() ==
                    pdeBuiltinCount + 1,
                "PDE bibliography addition did not persist");
        // Nested hash-based metadata must remain readable beyond Win32 MAX_PATH.
        const auto deepRoot =
            root / std::string(80, 'a') / std::string(80, 'b') / std::string(80, 'c');
        pdeQuery.object()["action"] = "add";
        const auto deepAdded = evaluateResearchCatalogueQuery(pdeQuery, deepRoot);
        require(deepAdded.at("status").string() == "added", "Deep catalogue write failed");
        require(deepAdded.at("added_path").string().size() > 260,
                "Regression fixture did not exceed the legacy path limit");
        pdeQuery.object()["action"] = "list";
        require(evaluateResearchCatalogueQuery(pdeQuery, deepRoot).at("papers").array().size() ==
                    pdeBuiltinCount + 1,
                "Deep catalogue could not be reloaded");
        pdeQuery.object()["action"] = "add";
        require(evaluateResearchCatalogueQuery(pdeQuery, deepRoot).at("status").string() ==
                    "already_present",
                "Deep catalogue duplicate was not detected");
        phase.at("options").object()["lower"] = j::parse("[-2,-2]");
        phase.object()["mode"] = "compute";
        require(evaluatePhaseSpaceQuery(phase, deepRoot).at("cache").at("status").string() ==
                    "computed",
                "Deep phase cache write failed");
        phase.object()["mode"] = "lookup";
        require(evaluatePhaseSpaceQuery(phase, deepRoot).at("cache").at("status").string() == "hit",
                "Deep phase cache was not reused");
        const auto solveJob = j::parse(
            R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"bibliography-ownership",
            "system":{"name":"harmonic_oscillator"},
            "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
            "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,"values":[],"ranges":[]},
            "controller":{"kind":"none"},"solver":{"method":"RK45","rtol":1e-8,"atol":1e-10}})");
        const auto solved = executeSolve(solveJob, root);
        require(solved.archives.size() == 1, "Solution archive missing");
        const auto archive = dynsyslab::core::readNpz(solved.archives.front());
        const auto metadata = j::parse(archive.metadataJson());
        require(metadata.at("dynamical_knowledge").isNull(),
                "Bibliography copied into solution metadata");
        require(metadata.at("analytic_solution").at("solution").at("references").isNull(),
                "Analytic solution archive contains paper references");
        require(archive.contains("analytic_values"),
                "Raw analytic trajectory missing from solution archive");
        require(std::filesystem::weakly_canonical(root).parent_path() ==
                        std::filesystem::weakly_canonical(std::filesystem::temp_directory_path()) &&
                    root.filename().string().rfind("dynsyslab-research-test-", 0) == 0,
                "Refusing cleanup outside the test's generated temporary directory");
        std::filesystem::remove_all(detail::metadataIoPath(root));
        std::cout
            << "Research catalogue persistence, duplicates, ownership and cache checks passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << "\nEvidence retained in " << root << '\n';
        return 1;
    }
}
