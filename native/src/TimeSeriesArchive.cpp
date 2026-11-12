#include <dynsyslab/core/TimeSeriesArchive.hpp>

#include <stdexcept>

namespace dynsyslab::core {
namespace {

const NpyArray* select(const NpzRecord& archive,
                       const std::vector<std::string>& names) {
    for (const auto& name : names) {
        const auto found = archive.arrays.find(name);
        if (found != archive.arrays.end()) return &found->second;
    }
    return nullptr;
}

}  // namespace

TimeSeries loadTimeSeries(const NpzRecord& archive,
                          const TimeSeriesArchiveSelection& selection) {
    const auto* time = select(archive, selection.timeArrays);
    const auto* values = select(archive, selection.valueArrays);
    if (!time || time->shape.size() != 1)
        throw std::invalid_argument("TimeSeries archive requires a one-dimensional time array");
    if (!values || values->shape.empty() || values->shape.size() > 2)
        throw std::invalid_argument("TimeSeries archive values must be a vector or matrix");
    const std::size_t samples = time->elementCount();
    const std::size_t components = values->shape.size() == 1 ? 1 : values->shape[1];
    if (values->shape[0] != samples)
        throw std::invalid_argument("TimeSeries archive time and value sample counts differ");
    return TimeSeries(time->doubles(), values->doubles(), components);
}

void storeTimeSeries(NpzRecord& archive, const TimeSeries& series,
                     const std::string& timeArray, const std::string& valueArray) {
    if (valueArray.empty()) throw std::invalid_argument("TimeSeries value array name is empty");
    if (!timeArray.empty())
        archive.arrays[timeArray] = NpyArray::fromDoubles(
            {series.sampleCount()}, series.timeData(), series.sampleCount());
    archive.arrays[valueArray] = NpyArray::fromDoubles(
        {series.sampleCount(), series.componentCount()}, series.valuesData(),
        series.sampleCount() * series.componentCount());
}

}  // namespace dynsyslab::core
