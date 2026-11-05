#pragma once

#include "dynsyslab_viewer/Theme.hpp"
#include <dynsyslab/core/ArchiveAudit.hpp>

#include <windows.h>

#include <optional>

namespace dsl::compute {

struct ArchiveToolRequest {
    enum class Action { audit, deduplicatePreview, deduplicateApply };
    Action action{Action::audit};
    dynsyslab::core::ArchiveAuditScope scope;
};

std::optional<ArchiveToolRequest> showArchiveToolsDialog(
    HWND owner, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
