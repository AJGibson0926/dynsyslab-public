#pragma once

#include <dynsyslab/core/EigenfunctionValidation.hpp>

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {


EigenfunctionValidationResult validateKoopmanSemigroup(
    const KoopmanModel& model,
    const std::vector<double>& time,
    const std::vector<double>& states,
    std::size_t samples,
    std::size_t stateDimension,
    std::vector<double> horizons = {},
    const std::vector<double>& trainingMinimum = {},
    const std::vector<double>& trainingMaximum = {});

}  // namespace dynsyslab::core
