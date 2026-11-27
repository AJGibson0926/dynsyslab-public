#pragma once
#include "DisplayCoordinates.hpp"
#include <dynsyslab/core/KoopmanModel.hpp>
namespace dsl {
struct KromTraining {
    dynsyslab::core::NpyArray states;
    std::vector<double> time;
    std::vector<std::size_t> starts;
    std::vector<std::string> labels, warnings;
    json::Value metadata;
    DisplayCoordinateTransform coordinates;
    bool sampleIndex = false;
};
KromTraining readKromTraining(const std::filesystem::path&, const dynsyslab::core::NpzRecord&,
    const json::Value&, const dynsyslab::core::KoopmanModel&);
dynsyslab::core::KoopmanModel loadKromEigenfunctions(const std::filesystem::path&,
    const dynsyslab::core::NpzRecord&, const json::Value&);
double kromSampleTime(const json::Value&);
// A krom as shown: a child (one that references a parent krom) as the whole K-ROM it describes, with
// its ancestors' method arrays (training samples, diagnostics) where it has none; others unchanged.
dynsyslab::core::NpzRecord composedKrom(const std::filesystem::path& source, dynsyslab::core::NpzRecord record);
}
