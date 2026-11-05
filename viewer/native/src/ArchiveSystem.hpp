#pragma once
#include "dynsyslab_viewer/Model.hpp"
#include <dynsyslab/core/System.hpp>
namespace dsl {
dynsyslab::core::SystemConfiguration archiveSystemConfiguration(const json::Value &metadata);
void appendVerifiedEquilibria(Figure &figure, const json::Value &metadata);
} // namespace dsl
