// Viewer payload for dynsyslab.animation/v1 manifests (the `animation` job): any views of any
// archives strung together frame by frame. Each panel is one archive drawn in one view; frame k
// shows page k of every panel. A panel has either one page (shown in every frame) or the full frame
// count, so panels stay in place and consistent between frames. Figures that belong to no page
// (dashboardPage < 0) appear in every frame.
#include "AnimationPayload.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <stdexcept>

namespace dsl {
namespace c = dynsyslab::core;

namespace {
int pages(const Payload& p) {
    int n = 1;
    for (const auto& f : p.figures) n = std::max(n, f.dashboardPage + 1);
    return std::max<int>(n, static_cast<int>(p.pageLabels.size()));
}
}  // namespace

Payload buildAnimationPayload(const std::filesystem::path& source, const c::NpzRecord& record) {
    const auto meta = c::json::parse(record.metadataJson());
    std::vector<Payload> panels;
    for (const auto& spec : meta.at("panels").array()) {
        const auto path = std::filesystem::u8path(spec.at("source").string());
        const auto view = spec.at("view").string();
        if (view == "animation") throw std::invalid_argument("an animation panel cannot itself be an animation");
        const auto archive = c::readNpz(path);
        panels.push_back(buildNativePayload(path, view, archive, c::json::parse(archive.metadataJson())));
        const auto label = spec.at("label").string();
        if (!label.empty())
            for (auto& f : panels.back().figures) f.title = label + " — " + f.title;
    }
    if (panels.empty()) throw std::invalid_argument("an animation needs at least one panel");
    int frames = 1;
    for (const auto& p : panels) frames = std::max(frames, pages(p));
    for (const auto& p : panels)
        if (pages(p) != 1 && pages(p) != frames)
            throw std::invalid_argument("animation panels must have one page or " + std::to_string(frames) + " frames");
    Payload out{source.u8string(), "animation", meta.at("title").string(meta.at("label").string(source.stem().u8string())), 1};
    int columns = 1;
    for (const auto& p : panels) columns = std::max(columns, p.columns);
    out.columns = static_cast<int>(meta.at("columns").number(double(columns)));
    const Payload* timeline = nullptr;
    for (const auto& p : panels) if (pages(p) == frames && !timeline) timeline = &p;
    for (int k = 0; k < frames; ++k) {
        for (const auto& p : panels) {
            const int page = pages(p) == 1 ? 0 : k;
            for (const auto& f : p.figures) {
                if (f.dashboardPage >= 0 && f.dashboardPage != page) continue;
                auto copy = f;
                copy.dashboardPage = k;
                out.figures.push_back(std::move(copy));
            }
        }
        out.pageLabels.push_back(timeline && k < static_cast<int>(timeline->pageLabels.size())
            ? timeline->pageLabels[k] : "Frame " + std::to_string(k + 1));
    }
    for (const auto& p : panels) {
        if (out.complexDisplay.empty()) out.complexDisplay = p.complexDisplay;
        out.warnings.insert(out.warnings.end(), p.warnings.begin(), p.warnings.end());
    }
    Payload::Animation animation;
    animation.frameSeconds = meta.at("frame_seconds").number(animation.frameSeconds);
    animation.loop = meta.at("loop").boolean(true);
    animation.autoplay = meta.at("autoplay").boolean(true);
    out.animation = animation;
    return out;
}

}  // namespace dsl
