#pragma once

#include <dynsyslab/core/Json.hpp>

// Compatibility alias for the desktop presentation layer. New workflow code
// uses dynsyslab::core::json directly and has no Viewer dependency.
namespace dsl {
namespace json = dynsyslab::core::json;
}
