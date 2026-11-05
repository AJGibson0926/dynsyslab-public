#pragma once
#include "dynsyslab_viewer/Json.hpp"
namespace dsl {
// Read-only inspection of actual running Viewer instances. Does not activate
// windows or substitute a newly loaded catalog for the displayed tree.
json::Value runningViewerWindows();
// Requests the existing Viewer's normal catalog refresh without showing or
// activating its window. The live-window query reports the resulting tree.
bool reloadRunningViewer(unsigned long processId);
}
