#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <string>
#include <vector>

namespace dsl::compute {

// Persistent, job-wide numerical execution controls. Unlike the System page,
// this pane remains visible while inspecting any proposed hierarchy node.
class SolverExecutionPane {
public:
    bool create(HWND parent, const dsl::theme::Fonts& fonts);
    void show(bool visible, bool enabled);
    void layout(int x, int y, int width, int height);
    void resetForFamily(const std::string& family);
    void load(const std::string& solverJson);
    [[nodiscard]] std::string settings() const;
    void focus();

private:
    HWND parent_{};
    std::vector<HWND> controls_;
    std::vector<HWND> labels_;
    HWND section_{};
    HWND method_{};
    HWND internalStep_{};
    HWND relativeTolerance_{};
    HWND absoluteTolerance_{};
    HWND maximumStep_{};
    HWND maximumInternalSteps_{};
    HWND backend_{};
    std::string baseSettings_{"{}"};

};

}  // namespace dsl::compute
