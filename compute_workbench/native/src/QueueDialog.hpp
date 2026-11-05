#pragma once

#include "dynsyslab_compute/Backend.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <string>

namespace dsl::compute {

// Opens the durable native job queue. Returns true when the user asks the
// workbench to claim and execute the next queued job asynchronously.
bool showQueueDialog(HWND owner, Backend& backend, const std::string& activeJob,
                     bool taskActive, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
