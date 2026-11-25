#pragma once
// Selecting arrays from an archive record and sampling them for display (NativePayload.cpp and the
// payload builders beside it). Display sampling never changes scientific data.
#include <dynsyslab/core/Archive.hpp>

#include <cstddef>
#include <initializer_list>
#include <vector>

namespace dsl {

inline constexpr std::size_t kDisplaySamples = 2400;

// The first of `names` the record holds (with `dimensions` axes when nonzero), or null.
[[nodiscard]] const dynsyslab::core::NpyArray* select(const dynsyslab::core::NpzRecord& record,
    std::initializer_list<const char*> names, std::size_t dimensions = 0);
// At most `maximum` evenly spaced indices of `count` samples.
[[nodiscard]] std::vector<std::size_t> displayIndices(std::size_t count, std::size_t maximum = kDisplaySamples);
// Error-bounded display subset of a long trajectory (core PathSimplification): fast excursions keep
// their vertices, smooth stretches are thinned.
[[nodiscard]] std::vector<std::size_t> trajectoryDisplayIndices(const dynsyslab::core::NpyArray& states, std::size_t maximum);
// The array's values at `indices` (sample indices when the array is missing or short).
[[nodiscard]] std::vector<double> sampledVector(const dynsyslab::core::NpyArray* array, std::size_t count,
                                                const std::vector<std::size_t>& indices);
[[nodiscard]] std::vector<double> sampledColumn(const dynsyslab::core::NpyArray& array, std::size_t column,
                                                const std::vector<std::size_t>& indices);

}  // namespace dsl
