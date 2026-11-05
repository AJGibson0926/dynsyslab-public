#include "KromPayload.hpp"
#include "dynsyslab_viewer/KromWorkspace.hpp"
namespace dsl {
Payload buildKromDashboardPayload(const std::filesystem::path& source,
    const dynsyslab::core::NpzRecord&, const std::string&) {
    return buildKromWorkspace(source);
}
}
