#include "PngCompare.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
#include <objbase.h>
#include <wincodec.h>
#include <windows.h>
#endif

namespace dynsyslab::workflow::detail {

#ifdef _WIN32
namespace {

struct Image {
    UINT width{}, height{};
    std::vector<std::uint8_t> pixels;
};

template <class T>
struct Released {
    T* item{};
    ~Released() { if (item) item->Release(); }
};

void check(HRESULT result, const char* what) {
    if (FAILED(result)) throw std::runtime_error(std::string("PNG decode failed: ") + what);
}

Image decode(IWICImagingFactory* factory, const std::filesystem::path& path) {
    Released<IWICBitmapDecoder> decoder;
    check(factory->CreateDecoderFromFilename(path.c_str(), nullptr, GENERIC_READ,
                                             WICDecodeMetadataCacheOnDemand, &decoder.item), "open");
    Released<IWICBitmapFrameDecode> frame;
    check(decoder.item->GetFrame(0, &frame.item), "frame");
    Released<IWICFormatConverter> converter;
    check(factory->CreateFormatConverter(&converter.item), "converter");
    check(converter.item->Initialize(frame.item, GUID_WICPixelFormat32bppBGRA, WICBitmapDitherTypeNone,
                                     nullptr, 0.0, WICBitmapPaletteTypeCustom), "convert");
    Image image;
    check(converter.item->GetSize(&image.width, &image.height), "size");
    image.pixels.resize(std::size_t(image.width) * image.height * 4);
    check(converter.item->CopyPixels(nullptr, image.width * 4, static_cast<UINT>(image.pixels.size()),
                                     image.pixels.data()), "pixels");
    return image;
}

}  // namespace

PixelComparison comparePngPixels(const std::filesystem::path& first, const std::filesystem::path& second) {
    const HRESULT initialized = CoInitializeEx(nullptr, COINIT_MULTITHREADED);
    const bool balance = SUCCEEDED(initialized);
    PixelComparison result;
    {
        Released<IWICImagingFactory> factory;
        check(CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER,
                               IID_PPV_ARGS(&factory.item)), "factory");
        const auto a = decode(factory.item, first);
        const auto b = decode(factory.item, second);
        result.sameSize = a.width == b.width && a.height == b.height;
        if (result.sameSize)
            for (std::size_t k = 0; k + 3 < a.pixels.size(); k += 4) {
                int largest = 0;
                for (std::size_t c = 0; c < 4; ++c)
                    largest = std::max(largest, std::abs(int(a.pixels[k + c]) - int(b.pixels[k + c])));
                if (largest) { ++result.differing; result.maximumChannelDifference = std::max(result.maximumChannelDifference, largest); }
            }
    }
    if (balance) CoUninitialize();
    return result;
}
#else
PixelComparison comparePngPixels(const std::filesystem::path&, const std::filesystem::path&) {
    throw std::runtime_error("PNG pixel comparison needs the Windows Imaging Component");
}
#endif

}  // namespace dynsyslab::workflow::detail
