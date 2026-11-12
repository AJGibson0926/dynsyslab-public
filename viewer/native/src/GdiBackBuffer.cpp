#include "dynsyslab_viewer/GdiBackBuffer.hpp"

#include <algorithm>

namespace dsl {

GdiBackBuffer::~GdiBackBuffer() {
    release();
}

void GdiBackBuffer::release() {
    if (memory_ && previousBitmap_) SelectObject(memory_, previousBitmap_);
    if (bitmap_) DeleteObject(bitmap_);
    if (memory_) DeleteDC(memory_);
    memory_ = nullptr;
    bitmap_ = nullptr;
    previousBitmap_ = nullptr;
    width_ = 0;
    height_ = 0;
}

HDC GdiBackBuffer::target(HDC compatibleWith, int width, int height) {
    width = std::max(1, width);
    height = std::max(1, height);
    if (memory_ && bitmap_ && width_ == width && height_ == height)
        return memory_;

    release();
    memory_ = CreateCompatibleDC(compatibleWith);
    if (!memory_) return compatibleWith;

    BITMAPINFO description{};
    description.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    description.bmiHeader.biWidth = width;
    description.bmiHeader.biHeight = -height;
    description.bmiHeader.biPlanes = 1;
    description.bmiHeader.biBitCount = 32;
    description.bmiHeader.biCompression = BI_RGB;
    void* pixels = nullptr;
    bitmap_ = CreateDIBSection(compatibleWith, &description, DIB_RGB_COLORS,
                               &pixels, nullptr, 0);
    if (!bitmap_) {
        DeleteDC(memory_);
        memory_ = nullptr;
        return compatibleWith;
    }
    previousBitmap_ = SelectObject(memory_, bitmap_);
    width_ = width;
    height_ = height;
    return memory_;
}

void GdiBackBuffer::present(HDC target, const RECT& region) const {
    if (!memory_ || !bitmap_ || region.right <= region.left
        || region.bottom <= region.top) return;
    BitBlt(target, region.left, region.top, region.right - region.left,
           region.bottom - region.top, memory_, region.left, region.top,
           SRCCOPY);
}

}  // namespace dsl
