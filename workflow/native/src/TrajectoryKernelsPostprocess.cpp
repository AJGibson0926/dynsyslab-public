// Post-processing methods over one trajectory (docs/architecture/TRAJECTORY_KERNELS.md):
//
//   affine_dynamics  {"model": "<dynamical_model.npz>", "scale": 1}
//       A, c of an affine fitted vector field (core/AffineDynamics.hpp) and its exact propagation
//       from the source's first sample on the source's clock, with the error against the source.
//   prehistory       {"duration": 2.0, "step": 0 (source step), "maximum_step": 0 (step / 10),
//                     "roundtrip_tolerance": 1e-8}
//       The source preceded by a backward-integrated prehistory of its own (unforced) system,
//       verified by integrating forward again (core/TrajectoryAudit.hpp).
//   clock_audit      {"reference": "<dense.npz>", "window": [a, b], "array": "values", "reference_array": "values"}
//       The source's samples against a dense reference record of the same quantity.
#include "PostprocessInternal.hpp"
#include "TrajectoryKernelsPostprocess.hpp"

#include <dynsyslab/core/AffineDynamics.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/TrajectoryAudit.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::filesystem::path runAffineDynamics(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                        const std::string& label) {
    const auto source = loadSource(path);
    const auto modelPath = std::filesystem::u8path(settings.at("model").string());
    const auto modelArchive = c::readNpz(modelPath);
    if (!modelArchive.contains("dynamics_model_json")) throw std::invalid_argument("affine_dynamics needs a dynamical_model with dynamics_model_json");
    const auto model = c::SindyDynamicsModel::fromDefinition(json::parse(modelArchive.at("dynamics_model_json").string()));
    const auto equations = c::affineEquations(model.vectorField(), settings.at("scale").number(1.0));
    const auto n = equations.dimension;
    if (source.dimension() != n) throw std::invalid_argument("affine_dynamics: the source and the model differ in dimension");
    const std::vector<double> x0(source.values().begin(), source.values().begin() + std::ptrdiff_t(n));
    const auto rollout = c::affinePropagation(equations, x0, source.time());
    J::Array error;
    for (std::size_t j = 0; j < n; ++j) {
        double e = 0, r = 0;
        for (std::size_t k = 0; k < source.samples(); ++k) {
            e += std::pow(rollout[k * n + j] - source.values()[k * n + j], 2);
            r += std::pow(source.values()[k * n + j], 2);
        }
        error.push_back(r > 0 ? J(std::sqrt(e / r)) : J());
    }
    c::NpzRecord archive;
    archive.arrays["time"] = c::NpyArray::fromDoubles({source.samples()}, source.time());
    archive.arrays["values"] = c::NpyArray::fromDoubles({source.samples(), n}, rollout);
    archive.arrays["A"] = c::NpyArray::fromDoubles({n, n}, equations.A);
    archive.arrays["c"] = c::NpyArray::fromDoubles({n}, equations.c);
    auto meta = json::parse(metadata(source, "affine-dynamics", label, "")).object();
    meta["model_archive"] = std::filesystem::absolute(modelPath).generic_u8string();
    meta["equations"] = "x' = A x + c";
    meta["relative_rms_error"] = error;
    meta["references"] = J::Array{"Van Loan 1978, IEEE Trans. Autom. Control 23:395-404"};
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::timeSeries, "affine-dynamics", std::move(archive));
}

std::filesystem::path runPrehistory(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                    const std::string& label) {
    const auto source = loadSource(path);
    const auto configuration = systemConfiguration(source);
    if (configuration.driver.amplitude != 0) throw std::invalid_argument("prehistory integrates the unforced system; the source is driven");
    const auto system = c::createSystem(configuration);
    const auto n = source.dimension();
    const double step = settings.at("step").number(0) > 0 ? settings.at("step").number() : uniformDeltaTime(source);
    const double duration = settings.at("duration").number(0);
    if (!(duration > 0)) throw std::invalid_argument("prehistory needs a positive duration");
    const auto count = std::size_t(std::floor(duration / step + 1e-9));
    const double t0 = source.time().front();
    std::vector<double> times;
    for (std::size_t k = count; k >= 1; --k) times.push_back(t0 - double(k) * step);
    const double maximum = settings.at("maximum_step").number(0) > 0 ? settings.at("maximum_step").number() : step / 10;
    const std::vector<double> x0(source.values().begin(), source.values().begin() + std::ptrdiff_t(n));
    const auto history = c::backwardHistory(*system, x0, t0, times, maximum);
    const double tolerance = settings.at("roundtrip_tolerance").number(1e-8);
    if (history.roundtripError > tolerance)
        throw std::runtime_error("prehistory: backward and forward integration differ by " + std::to_string(history.roundtripError));
    auto time = times;
    time.insert(time.end(), source.time().begin(), source.time().end());
    auto values = history.states;
    values.insert(values.end(), source.values().begin(), source.values().end());
    c::NpzRecord archive;
    archive.arrays["time"] = c::NpyArray::fromDoubles({time.size()}, time);
    archive.arrays["values"] = c::NpyArray::fromDoubles({time.size(), n}, values);
    auto meta = json::parse(metadata(source, "prehistory", label, "")).object();
    for (const char* key : {"forcing_type", "archive_coordinate_system", "integration_form"})
        if (!json::parse(source.metadataJson).at(key).isNull()) meta[key] = json::parse(source.metadataJson).at(key);
    meta["prehistory_samples"] = double(count);
    meta["roundtrip_error"] = history.roundtripError;
    meta["prehistory_note"] = "simulated: generated by the source's own model, not measured";
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::timeSeries, "prehistory", std::move(archive));
}

std::filesystem::path runClockAudit(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                    const std::string& label) {
    const auto source = loadSource(path);
    const auto sampled = c::readNpz(path), dense = c::readNpz(std::filesystem::u8path(settings.at("reference").string()));
    const auto array = settings.at("array").string("values"), referenceArray = settings.at("reference_array").string("values");
    const auto& s = sampled.at(array);
    const auto& r = dense.at(referenceArray);
    const auto width = s.shape.size() == 2 ? s.shape[1] : 1;
    if ((r.shape.size() == 2 ? r.shape[1] : 1) != width) throw std::invalid_argument("clock_audit: the records differ in width");
    double lo = -1e300, hi = 1e300;
    if (!settings.at("window").isNull()) { lo = settings.at("window").array().at(0).number(); hi = settings.at("window").array().at(1).number(); }
    const auto result = c::compareClocks(dense.at("time").doubles(), r.doubles(), sampled.at("time").doubles(), s.doubles(), width, lo, hi);
    const auto list = [](const std::vector<double>& v) { J::Array a; for (const auto x : v) a.push_back(x); return J(a); };
    auto meta = json::parse(metadata(source, "clock-audit", label, "")).object();
    meta["clock_audit"] = J::Object{{"compared_time", J::Array{result.first, result.last}}, {"reference_count", double(result.referenceCount)},
        {"sample_count", double(result.sampleCount)}, {"max_abs", list(result.maximumAbsolute)}, {"rms", list(result.rms)},
        {"worst_time", list(result.worstTime)}, {"reference", settings.at("reference")}};
    c::NpzRecord archive;
    archive.arrays["max_abs"] = c::NpyArray::fromDoubles({width}, result.maximumAbsolute);
    archive.arrays["rms"] = c::NpyArray::fromDoubles({width}, result.rms);
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::report, "clock-audit", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
