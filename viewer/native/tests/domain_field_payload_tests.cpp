// Viewer figures for grid-field archives: complex scalars on 2-D domains default to domain coloring
// and switch in place; push-forward status is kept apart from the values; frames become pages; a
// 3-D grid gives staggered slices and a translucent 3-D slice view.
#include "../src/ComplexHeatmap.hpp"
#include "../src/AnimationPayload.hpp"
#include "../src/DomainFieldPayload.hpp"
#include <dynsyslab/core/GridFieldArchive.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
c::GridField field(std::size_t dims, bool complex) {
    std::vector<double> lower(dims, -1), upper(dims, 1);
    auto domain = c::FunctionDomain::physicalState(dims);
    auto grid = c::DomainGrid::box(domain, lower, upper, 9);
    auto f = c::evaluateOnGrid([complex](const double* x) {
        return std::vector<std::complex<double>>{complex ? std::complex<double>(x[0], x[1]) : std::complex<double>(x[0] * x[1])};
    }, 1, grid);
    f.complexValued = complex;
    f.outputNames = {"g"};
    return f;
}
dsl::Payload payload(std::vector<c::GridField> frames, const std::string& operation = "evaluate") {
    c::GridFieldSeries series;
    for (std::size_t k = 0; k < frames.size(); ++k) series.frameTime.push_back(0.5 * double(k));
    series.frames = std::move(frames);
    series.metadata = c::json::Value::Object{{"label", "test"},
        {"domain_field_settings", c::json::Value::Object{{"operation", operation}}}};
    return dsl::buildDomainFieldPayload("test.npz", c::gridFieldArchive(series));
}
}  // namespace

int main() {
    try {
        // Complex scalar on a 2-D domain: domain coloring by default, switchable in place.
        auto complexMap = payload({field(2, true)});
        require(complexMap.complexDisplay == "domain-coloring", "domain coloring is the default");
        require(complexMap.figures.size() == 2 && complexMap.figures[0].traces.front().encoding == "domain-color",
                "domain-colored panel and magnitude panel");
        for (auto& f : complexMap.figures) dsl::encodeComplexHeatmap(f, "real-imaginary");
        require(complexMap.figures[1].traces.front().encoding == "scalar"
                    && std::abs(complexMap.figures[1].traces.front().z.back() - 1.0) < 1e-12, "imaginary panel after switching");

        // Real scalar: one heat map.
        require(payload({field(2, false)}).figures.size() == 1, "real scalar is one heat map");

        // Push-forward status stays apart from the values: conflicts are blank, status has its own panel.
        auto pushed = field(2, false);
        pushed.status[3] = c::GridStatus::conflict;
        pushed.spread.assign(pushed.status.size(), 0.0);
        pushed.counts.assign(pushed.status.size(), 1);
        const auto statusPayload = payload({pushed}, "push_forward");
        require(std::isnan(statusPayload.figures.front().traces.front().z[3]), "conflict is blank in the value panel");
        bool statusPanel = false;
        for (const auto& f : statusPayload.figures) statusPanel = statusPanel || f.title.rfind("Status", 0) == 0;
        require(statusPanel && statusPayload.figures.size() == 3, "value, spread and status panels");

        // Frames are pages.
        const auto frames = payload({field(2, false), field(2, false), field(2, false)}, "koopman_action");
        require(frames.pageLabels.size() == 3 && frames.figures.back().dashboardPage == 2, "one page per frame");

        // 3-D grid: five slices plus the translucent stacked view.
        const auto cube = payload({field(3, false)});
        std::size_t slicesSeen = 0;
        bool stacked = false;
        for (const auto& f : cube.figures) {
            if (f.kind == "domain_field") ++slicesSeen;
            if (f.kind == "three_d") stacked = f.traces.size() == 5 && f.traces.front().encoding == "slice";
        }
        require(slicesSeen == 5 && stacked, "3-D slices");

        // Frame series are animations; an animation manifest strings panels into frames: a
        // three-frame panel and a one-page panel give three frames that repeat the static panel.
        require(frames.animation && frames.animation->autoplay, "frame series autoplay");
        const auto dir = std::filesystem::temp_directory_path() / "dynsyslab-animation-payload-test";
        std::filesystem::create_directories(dir);
        const auto write = [&](const std::string& name, std::vector<c::GridField> fields) {
            c::GridFieldSeries series;
            for (std::size_t k = 0; k < fields.size(); ++k) series.frameTime.push_back(double(k));
            series.frames = std::move(fields);
            series.metadata = c::json::Value::Object{{"label", name}};
            const auto path = dir / (name + ".npz");
            c::writeNpzAtomic(path, c::gridFieldArchive(series));
            return path.generic_u8string();
        };
        const auto moving = write("moving", {field(2, false), field(2, false), field(2, false)});
        const auto still = write("still", {field(2, false)});
        c::NpzRecord manifest;
        manifest.arrays["metadata_json"] = c::NpyArray::fromString(c::json::dump(c::json::Value(c::json::Value::Object{
            {"format", "dynsyslab.animation/v1"}, {"title", "demo"}, {"frame_seconds", 0.1},
            {"panels", c::json::Value::Array{c::json::Value::Object{{"source", moving}, {"view", "domain_field"}},
                                             c::json::Value::Object{{"source", still}, {"view", "domain_field"}}}}})));
        const auto movie = dsl::buildAnimationPayload("demo.npz", manifest);
        require(movie.pageLabels.size() == 3 && movie.figures.size() == 6, "two panels in each of three frames");
        require(movie.animation && std::abs(movie.animation->frameSeconds - 0.1) < 1e-12, "animation timing");
        require(movie.figures[5].dashboardPage == 2, "static panel repeats in the last frame");
        std::filesystem::remove_all(dir);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "domain field payload tests passed\n";
    return 0;
}
