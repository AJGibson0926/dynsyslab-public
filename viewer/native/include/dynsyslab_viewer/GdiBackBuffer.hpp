#pragma once

#include <windows.h>

namespace dsl {

// Persistent off-screen surface for small Win32 owner-drawn controls.  A
// complete frame is drawn away from the display and copied to the window in
// one operation, avoiding intermediate clear/layout frames during scrolling.
class GdiBackBuffer {
public:
    GdiBackBuffer() = default;
    ~GdiBackBuffer();
    GdiBackBuffer(const GdiBackBuffer&) = delete;
    GdiBackBuffer& operator=(const GdiBackBuffer&) = delete;

    HDC target(HDC compatibleWith, int width, int height);
    void present(HDC target, const RECT& region) const;

private:
    void release();

    HDC memory_ = nullptr;
    HBITMAP bitmap_ = nullptr;
    HGDIOBJ previousBitmap_ = nullptr;
    int width_ = 0;
    int height_ = 0;
};

}  // namespace dsl
