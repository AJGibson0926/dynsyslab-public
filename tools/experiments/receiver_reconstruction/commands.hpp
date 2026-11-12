#pragma once
#include <filesystem>

void receiverSignals(const std::filesystem::path& source,const std::filesystem::path& geometry,
                     const std::filesystem::path& output);
void ringErrors(const std::filesystem::path& reference,const std::filesystem::path& predictions,
                const std::filesystem::path& output);
