#include "PayloadSampling.hpp"
#include <dynsyslab/core/PathSimplification.hpp>

namespace dsl {
using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;

const NpyArray* select(const NpzRecord& record,
                       std::initializer_list<const char*> names,
                       std::size_t dimensions) {
    for (const char* name : names) {
        const auto found = record.arrays.find(name);
        if (found != record.arrays.end() && (!dimensions || found->second.shape.size() == dimensions))
            return &found->second;
    }
    return nullptr;
}

std::vector<std::size_t> displayIndices(std::size_t count, std::size_t maximum) {
    if (count <= maximum) {
        std::vector<std::size_t> result(count);
        for (std::size_t index = 0; index < count; ++index) result[index] = index;
        return result;
    }
    std::vector<std::size_t> result;
    result.reserve(maximum);
    for (std::size_t index = 0; index < maximum; ++index)
        result.push_back(index * (count - 1) / (maximum - 1));
    return result;
}

// Error-bounded display subset of a long trajectory (core PathSimplification,
// after the trajectory-resolution side-car): fast excursions keep their
// vertices, smooth stretches are thinned; scientific data is untouched.
std::vector<std::size_t> trajectoryDisplayIndices(const NpyArray& states, std::size_t maximum) {
    const std::size_t count = states.shape[0], width = states.shape.size() > 1 ? states.shape[1] : 1;
    if (count <= maximum || width == 0) return displayIndices(count, maximum);
    return dynsyslab::core::simplifiedPathIndices(states.doubles(), count, width, maximum);
}

std::vector<double> sampledVector(const NpyArray* array, std::size_t count,
                                  const std::vector<std::size_t>& indices) {
    std::vector<double> values;
    if (array) values = array->doubles();
    if (values.size() < count) {
        values.resize(count);
        for (std::size_t index = 0; index < count; ++index) values[index] = static_cast<double>(index);
    }
    std::vector<double> result;
    result.reserve(indices.size());
    for (const auto index : indices) result.push_back(values[index]);
    return result;
}

std::vector<double> sampledColumn(const NpyArray& array, std::size_t column,
                                  const std::vector<std::size_t>& indices) {
    if (array.shape.size() != 2 || column >= array.shape[1]) return {};
    const auto values = array.doubles();
    const std::size_t columns = array.shape[1];
    std::vector<double> result;
    result.reserve(indices.size());
    for (const auto row : indices) result.push_back(values[row * columns + column]);
    return result;
}

}  // namespace dsl
