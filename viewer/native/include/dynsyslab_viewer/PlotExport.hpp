#pragma once

#include <windows.h>

#include <filesystem>
#include <string>

namespace dsl {

struct PlotExportResult {
    bool succeeded = false;
    bool cancelled = false;
    std::wstring detail;
};

PlotExportResult copyPlotToClipboard(HWND owner, HWND plotWindow);
PlotExportResult savePlotAsPng(HWND owner, HWND plotWindow,
                               const std::wstring& suggestedStem);
PlotExportResult savePlotToDisk(HWND owner, HWND plotWindow,
                                const std::wstring& suggestedStem);
bool saveBitmapAsPdf(const std::filesystem::path& path,HBITMAP bitmap,
                     int width,int height);

}  // namespace dsl
