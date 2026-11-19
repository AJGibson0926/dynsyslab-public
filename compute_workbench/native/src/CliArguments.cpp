#include "CliArguments.hpp"
#include "CliRuntimePaths.hpp"

#include <dynsyslab/core/Text.hpp>

#include <stdexcept>

namespace dynsyslab::cli {
Arguments arguments(int count, wchar_t** values) {
    Arguments result;
    for (int index = 1; index < count; ++index) {
        const std::wstring value = values[index];
        if ((value == L"--archive-root" || value == L"--progress") && index + 1 < count) {
            const auto path = std::filesystem::path(values[++index]);
            if (value == L"--archive-root") result.archiveRoot = path;
            else result.progress = path;
        } else if (value == L"--days" && index + 1 < count) {
            result.auditDays = std::stod(values[++index]);
        } else if (value == L"--help" || value == L"-h" || value == L"/?") {
            result.help = true;
        } else if (value == L"--format" && index + 1 < count) {
            result.format = dynsyslab::core::narrow(values[++index]);
        } else if (value == L"--all-ages") {
            result.auditDays = -1.0;
        } else if ((value == L"--section" || value == L"--system") && index + 1 < count) {
            auto name = dynsyslab::core::narrow(values[++index]);
            if (value == L"--section") result.auditSections.push_back(std::move(name));
            else result.auditSystems.push_back(std::move(name));
        } else if (value == L"--include-derived") {
            result.auditDerived = true;
        } else if (value == L"--apply") {
            result.apply = true;
        } else if ((value == L"--map-system" || value == L"--preset"
                    || value == L"--coordinates" || value == L"--label"
                    || value == L"--parameters-json") && index + 1 < count) {
            const auto item=dynsyslab::core::narrow(values[++index]);
            if(value==L"--map-system")result.importSystem=item;
            else if(value==L"--preset")result.importPreset=item;
            else if(value==L"--coordinates")result.importCoordinates=item;
            else if(value==L"--label")result.importLabel=item;
            else result.importParameters=item;
        } else if (value == L"--fresh") {
            result.sets.emplace_back("__fresh", "true");
        } else if (value == L"--refresh" && index + 1 < count) {
            result.sets.emplace_back("__refresh", dynsyslab::core::narrow(values[++index]));
        } else if (value == L"--until" && index + 1 < count) {
            result.sets.emplace_back("__until", dynsyslab::core::narrow(values[++index]));
        } else if ((value == L"--param" || value == L"--set") && index + 1 < count) {
            const auto item = dynsyslab::core::narrow(values[++index]);
            const auto equals = item.find('=');
            if (equals == std::string::npos || equals == 0)
                throw std::invalid_argument("--param/--set expects NAME=VALUE");
            (value == L"--param" ? result.params : result.sets)
                .emplace_back(item.substr(0, equals), item.substr(equals + 1));
        } else if (value == L"--input-scale" && index + 1 < count) {
            result.inputScale=std::stod(values[++index]);
        } else if (value == L"--allow-unknown-records") {
            result.allowUnknownRecords=true;
        } else if (value == L"--component" && index + 1 < count) {
            result.fourierComponent=static_cast<std::size_t>(std::stoull(values[++index]));
        } else if (value == L"--window" && index + 1 < count) {
            result.fourierWindow=dynsyslab::core::narrow(values[++index]);
        } else if (value == L"--keep-mean") {
            result.fourierRemoveMean=false;
        } else if (value == L"--two-sided") {
            result.fourierOneSided=false;
        } else result.positional.push_back(value);
    }
    if (result.archiveRoot.empty()) {
        auto choice = defaultArchiveRoot();
        result.archiveRoot = std::move(choice.path);
        result.archiveRootSource = std::move(choice.source);
    }
    return result;
}
}  // namespace dynsyslab::cli
