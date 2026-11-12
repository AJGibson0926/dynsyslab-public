#pragma once
#include <filesystem>
void fixedRollout(const std::filesystem::path& model, const std::filesystem::path& query,
                  double step, const std::filesystem::path& destination);
