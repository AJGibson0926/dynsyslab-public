#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

struct KromControlSearchRequest {
    enum class Action { plan,run,queue,save };
    std::string jobJson;
    Action action{Action::run};
};

[[nodiscard]] std::optional<KromControlSearchRequest> showKromControlSearchDialog(
    HWND owner,const std::filesystem::path& selected,const Registry& registry,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
