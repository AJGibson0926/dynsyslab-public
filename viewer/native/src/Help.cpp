#include "dynsyslab_viewer/Help.hpp"
#include "dynsyslab_viewer/Process.hpp"

#include <sstream>

namespace dsl {

void showHelpAbout(HWND parent, const std::wstring& application,
                   const std::filesystem::path& projectRoot,
                   const std::filesystem::path& archiveRoot) {
    const auto guide = projectRoot / L"docs" / L"START_HERE.md";
    std::wostringstream text;
    text << application << L"\n\n"
         << L"Native UI support: " << kUiSupportVersion << L"\n"
         << L"Viewer catalog schema: dynsyslab.viewer.catalog/v1\n"
         << L"Compute job schema: dynsyslab.compute-job/v1\n\n"
         << L"Project: " << projectRoot.wstring() << L"\n"
         << L"Archives: " << archiveRoot.wstring() << L"\n\n";
    if (std::filesystem::is_regular_file(guide))
        text << L"Select Yes to open the Start Here guide.";
    else
        text << L"The source-checkout Start Here guide is not available.";
    const int choice = MessageBoxW(parent, text.str().c_str(), L"DynSysLab Help / About",
        MB_YESNO | MB_ICONINFORMATION);
    if (choice == IDYES && std::filesystem::is_regular_file(guide))
        openDocument(guide.wstring());
}

}  // namespace dsl
