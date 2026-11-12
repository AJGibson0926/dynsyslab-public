#pragma once

#include <string>
#include <vector>

namespace dsl {

struct ProcessResult {
    unsigned long exitCode = 0;
    std::wstring error;
};

ProcessResult runHidden(const std::vector<std::wstring>& arguments,
                        const std::wstring& workingDirectory);
ProcessResult launchDetached(const std::vector<std::wstring>& arguments,
                             const std::wstring& workingDirectory);
bool openDocument(const std::wstring& path);

}  // namespace dsl
