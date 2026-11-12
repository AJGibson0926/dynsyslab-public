#pragma once
#include "dynsyslab_viewer/Model.hpp"
#include <atomic>
namespace dsl {
json::Value callPhaseSidecar(const json::Value &request, const std::filesystem::path &archives,
                             const std::atomic_bool *cancelled = nullptr);
void appendPhaseSurfaces(Figure &, const json::Value &, std::size_t, double, double);
void appendPhaseContext(Figure &figure, const json::Value &metadata, const json::Value &result);
} // namespace dsl
