#pragma once
#include "dynsyslab_viewer/Theme.hpp"
#include <windows.h>
#include <optional>
#include <string>
namespace dsl::compute {
std::optional<std::string> showSpectralKoopmanDialog(HWND,const std::string&,const dsl::theme::Fonts&);
}
