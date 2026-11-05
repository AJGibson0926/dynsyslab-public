#include <dynsyslab/core/Archive.hpp>

#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

template <typename T>
std::vector<double> numeric(const NpyArray& array) {
    if (array.data.size() != array.elementCount() * sizeof(T))
        throw std::runtime_error("NPY array byte count does not match its shape and dtype");
    std::vector<double> result(array.elementCount());
    const auto* values = reinterpret_cast<const T*>(array.data.data());
    std::transform(values, values + result.size(), result.begin(),
                   [](T value) { return static_cast<double>(value); });
    return result;
}

}  // namespace

std::size_t NpyArray::elementCount() const {
    return shape.empty() ? 1 : std::accumulate(shape.begin(), shape.end(), std::size_t{1},
        [](std::size_t left, std::size_t right) { return left * right; });
}

std::vector<double> NpyArray::doubles() const {
    switch (type) {
    case ArrayType::float64: return numeric<double>(*this);
    case ArrayType::float32: return numeric<float>(*this);
    case ArrayType::int64: return numeric<std::int64_t>(*this);
    case ArrayType::int32: return numeric<std::int32_t>(*this);
    case ArrayType::uint64: return numeric<std::uint64_t>(*this);
    case ArrayType::uint32: return numeric<std::uint32_t>(*this);
    case ArrayType::uint8: return numeric<std::uint8_t>(*this);
    case ArrayType::boolean: return numeric<std::uint8_t>(*this);
    default: throw std::runtime_error("NPY array is not real numeric data");
    }
}

std::vector<std::complex<double>> NpyArray::complexDoubles() const {
    if (type != ArrayType::complex128 || data.size() != elementCount() * sizeof(std::complex<double>))
        throw std::runtime_error("NPY array is not complex128 data");
    std::vector<std::complex<double>> result(elementCount());
    std::memcpy(result.data(), data.data(), data.size());
    return result;
}

std::string NpyArray::string() const {
    if (type != ArrayType::unicode) throw std::runtime_error("NPY array is not Unicode data");
    const std::size_t characters = descriptor.size() > 2
        ? static_cast<std::size_t>(std::stoull(descriptor.substr(2))) : 0;
    if (data.size() < characters * 4) throw std::runtime_error("truncated NPY Unicode value");
    std::string result;
    for (std::size_t index = 0; index < characters; ++index) {
        std::uint32_t code = 0;
        std::memcpy(&code, data.data() + index * 4, 4);
        if (!code) break;
        if (code <= 0x7f) result.push_back(static_cast<char>(code));
        else if (code <= 0x7ff) {
            result.push_back(static_cast<char>(0xc0 | (code >> 6)));
            result.push_back(static_cast<char>(0x80 | (code & 0x3f)));
        } else if (code <= 0xffff) {
            result.push_back(static_cast<char>(0xe0 | (code >> 12)));
            result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3f)));
            result.push_back(static_cast<char>(0x80 | (code & 0x3f)));
        } else {
            result.push_back(static_cast<char>(0xf0 | (code >> 18)));
            result.push_back(static_cast<char>(0x80 | ((code >> 12) & 0x3f)));
            result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3f)));
            result.push_back(static_cast<char>(0x80 | (code & 0x3f)));
        }
    }
    return result;
}

NpyArray NpyArray::fromDoubles(std::vector<std::size_t> shape,
                               const std::vector<double>& values) {
    return fromDoubles(std::move(shape), values.data(), values.size());
}

NpyArray NpyArray::fromDoubles(std::vector<std::size_t> shape,
                               const double* values, std::size_t count) {
    NpyArray result{ArrayType::float64, "<f8", std::move(shape), {}};
    if (result.elementCount() != count || (count && !values))
        throw std::invalid_argument("array shape mismatch");
    result.data.resize(count * sizeof(double));
    if (count) std::memcpy(result.data.data(), values, result.data.size());
    return result;
}

NpyArray NpyArray::fromInt64s(std::vector<std::size_t> shape,
                              const std::vector<std::int64_t>& values) {
    NpyArray result;
    result.type = ArrayType::int64;
    result.descriptor = "<i8";
    result.shape = std::move(shape);
    if (result.elementCount() != values.size()) throw std::invalid_argument("array shape mismatch");
    result.data.resize(values.size() * sizeof(std::int64_t));
    std::memcpy(result.data.data(), values.data(), result.data.size());
    return result;
}

NpyArray NpyArray::fromComplex(std::vector<std::size_t> shape,
                               const std::vector<std::complex<double>>& values) {
    NpyArray result{ArrayType::complex128, "<c16", std::move(shape), {}};
    if (result.elementCount() != values.size()) throw std::invalid_argument("array shape mismatch");
    result.data.resize(values.size() * sizeof(std::complex<double>));
    std::memcpy(result.data.data(), values.data(), result.data.size());
    return result;
}

NpyArray NpyArray::fromString(const std::string& value) {
    NpyArray result{ArrayType::unicode, "<U" + std::to_string(value.size()), {}, {}};
    result.data.resize(value.size() * 4, 0);
    for (std::size_t index = 0; index < value.size(); ++index) {
        const std::uint32_t code = static_cast<unsigned char>(value[index]);
        std::memcpy(result.data.data() + index * 4, &code, 4);
    }
    return result;
}

bool NpzRecord::contains(const std::string& name) const { return arrays.count(name) != 0; }

const NpyArray& NpzRecord::at(const std::string& name) const {
    const auto found = arrays.find(name);
    if (found == arrays.end()) throw std::out_of_range("NPZ array is missing: " + name);
    return found->second;
}

std::string NpzRecord::metadataJson() const {
    return contains("metadata_json") ? at("metadata_json").string() : std::string{};
}

}  // namespace dynsyslab::core
