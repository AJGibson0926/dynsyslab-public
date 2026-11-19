#include "dynsyslab_viewer/RecordWindow.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace dsl {

using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;

NpzRecord windowRecord(const NpzRecord& record, const TimeWindow& window) {
    if (!window.active()) return record;
    if ((window.start && !std::isfinite(*window.start)) || (window.stop && !std::isfinite(*window.stop)))
        throw std::invalid_argument("Expected finite time limits");
    if (window.start && window.stop && *window.stop <= *window.start)
        throw std::invalid_argument("The time window must have stop > start");
    // Leading length -> retained row range. The first clock of a given length wins.
    std::map<std::size_t, std::pair<std::size_t, std::size_t>> ranges;
    for (const char* name : {"time", "controlled_time", "control_time", "free_time", "estimator_time"}) {
        const auto found = record.arrays.find(name);
        if (found == record.arrays.end() || found->second.shape.size() != 1) continue;
        const auto clock = found->second.doubles();
        if (clock.size() < 2 || ranges.count(clock.size())) continue;
        const auto first = std::lower_bound(clock.begin(), clock.end(), window.start.value_or(clock.front()));
        const auto last = std::upper_bound(clock.begin(), clock.end(), window.stop.value_or(clock.back()));
        if (last - first < 2)
            throw std::invalid_argument(std::string("The time window keeps fewer than two samples of ") + name);
        ranges[clock.size()] = {static_cast<std::size_t>(first - clock.begin()),
                                static_cast<std::size_t>(last - clock.begin())};
    }
    if (ranges.empty()) throw std::invalid_argument("This archive has no clock to window");
    NpzRecord result;
    for (const auto& [name, array] : record.arrays) {
        const auto range = array.shape.empty() ? ranges.end() : ranges.find(array.shape[0]);
        if (range == ranges.end() || name == "metadata_json") { result.arrays.emplace(name, array); continue; }
        const auto [first, last] = range->second;
        const std::size_t rowBytes = array.data.size() / array.shape[0];
        NpyArray sliced;
        sliced.type = array.type; sliced.descriptor = array.descriptor; sliced.shape = array.shape;
        sliced.shape[0] = last - first;
        sliced.data.assign(array.data.begin() + static_cast<std::ptrdiff_t>(first * rowBytes),
                           array.data.begin() + static_cast<std::ptrdiff_t>(last * rowBytes));
        result.arrays.emplace(name, std::move(sliced));
    }
    return result;
}

}  // namespace dsl
