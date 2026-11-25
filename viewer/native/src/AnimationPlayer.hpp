#pragma once
// A self-contained HTML animation player: frames rendered by the Viewer canvas (PNG), played with
// the same controls as the Viewer toolbar (first, back, play/pause, forward, stop) plus a frame
// slider and speed choice. Written by `DynSysLabViewerCli animate`; the HTML opens in any browser
// and can be embedded in a chat or a report.
#include <string>
#include <vector>

namespace dsl {

struct AnimationPlayerFrame {
    std::string png;    // raw PNG bytes
    std::string label;  // e.g. "t = 0.5"
};

[[nodiscard]] std::string base64(const std::string& bytes);
[[nodiscard]] std::string animationPlayerHtml(const std::string& title, const std::vector<AnimationPlayerFrame>& frames,
                                              double frameSeconds, bool loop, bool autoplay);

}  // namespace dsl
