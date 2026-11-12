#include "dynsyslab_viewer/OrnamentalScrollBar.hpp"

#include <iostream>

namespace {

bool expect(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

}  // namespace

int main() {
    dsl::OrnamentalScrollBar bar;
    bool ok = true;
    bar.configure(1000, 200);
    ok &= expect(bar.visible(), "overflow should display the custom scrollbar");
    ok &= expect(bar.scrollBy(125) && bar.position() == 125,
        "scrolling should advance through content coordinates");
    ok &= expect(bar.scrollBy(10000) && bar.position() == 800,
        "scrolling should clamp at the content boundary");
    ok &= expect(!bar.scrollBy(1) && bar.position() == 800,
        "scrolling at the boundary should report no change");
    bar.configure(100, 200);
    ok &= expect(!bar.visible() && bar.position() == 0,
        "a fitting document should hide the scrollbar and reset its position");
    return ok ? 0 : 1;
}
