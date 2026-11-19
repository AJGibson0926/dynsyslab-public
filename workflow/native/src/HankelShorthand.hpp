#pragma once
#include "PostprocessInternal.hpp"

namespace dynsyslab::workflow::postprocess_detail {
// Rewrites Hankel methods that name physics observables into an explicit
// coordinate_transform {"kind":"channels"} step plus the same Hankel method
// reading that step's columns. Jobs already in explicit form are unchanged.
json::Value expandHankelCoordinateShorthand(const json::Value& postprocessJob);
} // namespace dynsyslab::workflow::postprocess_detail
