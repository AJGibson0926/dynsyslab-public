// User-supplied output paths: bare file names resolve against the working directory, missing
// parent folders are created, and failures name the option and the resolved path.
#include <dynsyslab/core/OutputPath.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {
namespace fs = std::filesystem;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
}  // namespace

int main() {
    const auto previous = fs::current_path();
    const auto root = fs::temp_directory_path() / ("dynsyslab-output-path-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        fs::create_directories(root);
        fs::current_path(root);
        using dynsyslab::core::prepareOutputFile;

        const auto bare = prepareOutputFile("progress.jsonl", "--progress");  // the reported bug
        require(bare.is_absolute() && fs::equivalent(bare.parent_path(), root), "bare file name resolves to the cwd");
        const auto nested = prepareOutputFile(fs::path("sub") / "deeper" / "progress.jsonl", "--progress");
        require(fs::is_directory(root / "sub" / "deeper") && nested.is_absolute(), "relative sub-folder created");
        const auto absolute = prepareOutputFile(root / "abs" / "out.npz", "--output");
        require(fs::is_directory(root / "abs") && absolute == (root / "abs" / "out.npz").lexically_normal(), "absolute path");
        require(dynsyslab::core::prepareOutputDirectory("archives", "--archive-root") == (root / "archives").lexically_normal()
                && fs::is_directory(root / "archives"), "relative archive folder created");

        std::ofstream(root / "blocker") << "a file, not a folder";
        bool rejected = false;
        try { (void)prepareOutputFile(fs::path("blocker") / "progress.jsonl", "--progress"); }
        catch (const dynsyslab::core::OutputPathError& error) {
            const std::string text = error.what();
            rejected = error.option() == "--progress" && text.rfind("--progress: cannot create parent directory of", 0) == 0
                && text.find("blocker") != std::string::npos && error.path().is_absolute();
        }
        require(rejected, "an uncreatable parent is a clear error naming the option and resolved path");
        bool empty = false;
        try { (void)prepareOutputFile("", "--progress"); } catch (const dynsyslab::core::OutputPathError&) { empty = true; }
        require(empty, "empty output path rejected");

        fs::current_path(previous);
        fs::remove_all(root);
        std::cout << "dynsyslab_output_path_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored; fs::current_path(previous, ignored); fs::remove_all(root, ignored);
        std::cerr << "dynsyslab_output_path_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
