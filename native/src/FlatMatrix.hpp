#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dynsyslab::detail {

template<class Scalar>
struct FlatMatrix {
    FlatMatrix()=default;
    FlatMatrix(std::size_t rowCount,std::size_t columnCount)
        : values(rowCount*columnCount),rows(rowCount),columns(columnCount) {}
    FlatMatrix(std::vector<Scalar> data,std::size_t rowCount,
               std::size_t columnCount)
        : values(std::move(data)),rows(rowCount),columns(columnCount) {
        if(!rows||!columns||values.size()!=rows*columns)
            throw std::invalid_argument("flat matrix shape does not match its data");
    }

    Scalar& operator()(std::size_t row,std::size_t column) {
        return values[row*columns+column];
    }
    const Scalar& operator()(std::size_t row,std::size_t column) const {
        return values[row*columns+column];
    }

    std::vector<Scalar> values;
    std::size_t rows{};
    std::size_t columns{};
};

template<class Scalar>
double squaredNorm(const FlatMatrix<Scalar>& matrix) {
    double result=0.0;
    for(const auto& value:matrix.values) {
        const double magnitude=std::abs(value);result+=magnitude*magnitude;
    }
    return result;
}

}  // namespace dynsyslab::detail
