#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>
#include <string>

namespace dsl::compute {

HWND makeWaveformPreview(HWND parent, const dsl::theme::Fonts& fonts);
void setWaveformPreview(HWND preview, const std::string& settings,
                        double startTime, double stopTime);

}  // namespace dsl::compute
