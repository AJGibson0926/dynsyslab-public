#include "dynsyslab_viewer/ArchiveFocus.hpp"
#include "dynsyslab_viewer/RuntimePaths.hpp"
#include "dynsyslab_viewer/StartupDiagnostics.hpp"

#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>

namespace {

bool expect(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / L"dynsyslab-focus-test";
    std::filesystem::create_directories(root / L"archives");
    const auto manifest = root / L"focus.json";
    std::ofstream output(manifest);
    output << R"({"schema":"dynsyslab.viewer.archive-focus/v1","active":true,)"
              R"("label":"focused test","allowed_sources":["archives/kept.npz"]})";
    output.close();

    dsl::Catalog catalog;
    dsl::CatalogItem batch; batch.id = "batch"; batch.source = "hierarchy.json";
    dsl::CatalogItem system; system.id = "system"; system.parentId = "batch";
    system.source = "hierarchy.json";
    dsl::CatalogItem kept; kept.id = "kept"; kept.parentId = "system";
    kept.source = "archives/kept.npz";
    kept.relations.push_back({"hidden", "derived from", "hidden"});
    dsl::CatalogItem hidden; hidden.id = "hidden"; hidden.source = "archives/hidden.npz";
    catalog.items = {batch, system, kept, hidden};

    dsl::applyArchiveFocus(catalog, manifest, root);
    bool ok = expect(catalog.items.size() == 3, "focus did not retain leaf ancestry only");
    ok &= expect(catalog.indexOf("batch") >= 0 && catalog.indexOf("system") >= 0
                 && catalog.indexOf("kept") >= 0 && catalog.indexOf("hidden") < 0,
                 "focus retained the wrong catalog objects");
    ok &= expect(catalog.items[static_cast<std::size_t>(catalog.indexOf("kept"))]
                     .relations.empty(), "focus retained a dangling relation");
    ok &= expect(catalog.focusLabel == "focused test" && catalog.hiddenItemCount == 1,
                 "focus diagnostics were not recorded");

    dsl::Catalog noMatch;
    noMatch.items = {hidden};
    dsl::applyArchiveFocus(noMatch, manifest, root);
    ok &= expect(noMatch.items.empty() && noMatch.focusLabel == "focused test"
        && noMatch.hiddenItemCount == 1, "an unmatched Mask silently exposed the full tree");

    const auto invalid = root / "invalid.json";
    for (const auto& document : {"{", "{}",
         "{\"schema\":\"dynsyslab.viewer.archive-focus/v1\",\"active\":true}",
         "{\"schema\":\"dynsyslab.viewer.archive-focus/v1\",\"active\":true,\"allowed_sources\":[12]}"}) {
        std::ofstream(invalid) << document;
        dsl::Catalog untouched; untouched.items = {hidden};
        bool rejected = false;
        try { dsl::applyArchiveFocus(untouched, invalid, root); }
        catch (const std::exception&) { rejected = true; }
        ok &= expect(rejected && untouched.items.size() == 1,
            "an invalid Mask must fail before the catalog can be displayed");
    }
    std::ofstream(invalid) << R"({"schema":"dynsyslab.viewer.archive-focus/v1","active":true,"label":"empty","allowed_sources":[]})";
    dsl::Catalog emptyMask; emptyMask.items = {hidden};
    dsl::applyArchiveFocus(emptyMask, invalid, root);
    ok &= expect(emptyMask.items.empty() && emptyMask.focusLabel == "empty",
        "an empty active Mask silently exposed the full tree");
    std::ofstream(invalid) << R"({"schema":"dynsyslab.viewer.archive-focus/v1","active":false})";
    dsl::Catalog disabled; disabled.items = {hidden};
    dsl::applyArchiveFocus(disabled, invalid, root);
    ok &= expect(disabled.items.size() == 1 && disabled.focusLabel.empty(),
        "an explicitly disabled Mask must permit the full tree");

    const auto preferences = root / "user" / "Viewer";
    std::filesystem::create_directories(preferences);
    std::filesystem::copy_file(manifest, preferences / "archive-focus.json",
        std::filesystem::copy_options::overwrite_existing);
    for (const auto& cache : {root / "custom-cache", root / "Compute" / "viewer"}) {
        const auto paths = dsl::resolveRuntimePaths(L"Viewer", root / "bin" / "Viewer.exe",
            root / "user", root, root / "archives", cache);
        dsl::Catalog launched;
        launched.items = {batch, system, kept, hidden};
        const auto focus = paths.preferencesRoot / "archive-focus.json";
        dsl::applyArchiveFocus(launched, focus, paths.projectRoot);
        ok &= expect(launched.items.size() == 3 && launched.indexOf("hidden") < 0,
            "a custom or Compute cache dropped the saved Mask");
        dsl::writeStartupDiagnostics("Viewer", paths.cacheRoot, paths.projectRoot,
            paths.archiveRoot, {}, launched.items.size(), 3, focus,
            launched.focusLabel, launched.hiddenItemCount);
        const auto diagnostics = dynsyslab::core::json::load(
            (cache / "startup-diagnostics.json").wstring());
        ok &= expect(diagnostics.at("archive_focus_applied").boolean()
            && diagnostics.at("archive_focus_hidden_items").number() == 1
            && diagnostics.at("archive_focus_path").string() == focus.generic_u8string(),
            "startup diagnostics lost the applied Mask");
    }
    std::filesystem::remove_all(root);
    return ok ? 0 : 1;
}
