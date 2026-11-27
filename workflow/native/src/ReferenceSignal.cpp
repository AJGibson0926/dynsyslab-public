// Post-processing method "reference_signal": samples a controller target and the input that makes
// it an exact solution of a control-affine system, x' = f(t, x) + g(t, x) u.
//
// Settings: {"system":{...} or taken from the first source's metadata, "target":<waveform or
// constant target, optional ramp>, "time":{"start","stop","step"} or "clock":{"array":"time"}}.
// The target's time derivative is analytic. The input is the least-squares solution
// u = g^T (x' - f) / (g^T g), and matching_residual reports what g cannot reach.
// Output arrays: time, target_values, target_derivative, inverse_drive, matching_residual.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number());
    return out;
}

struct TargetSample { double position, velocity, acceleration; };

// Position, velocity and acceleration of a waveform target with an optional quintic ramp.
class TargetJet {
public:
    explicit TargetJet(const json::Value& target) {
        const auto kind = target.at("kind").string();
        if (kind == "constant") {
            const auto state = numbers(target.at("state"));
            offset_ = state.at(0); constantVelocity_ = state.size() > 1 ? state[1] : 0.0;
            return;
        }
        if (kind != "waveform") throw std::invalid_argument("reference_signal supports constant and waveform targets");
        offset_ = numbers(target.at("offset")).at(0);
        amplitudes_ = numbers(target.at("amplitudes"));
        frequencies_ = numbers(target.at("angular_frequencies"));
        if (frequencies_.empty()) {
            frequencies_ = numbers(target.at("frequencies"));
            for (auto& value : frequencies_) value *= 6.283185307179586476925286766559;
        }
        phases_ = numbers(target.at("phases_rad"));
        if (phases_.empty()) phases_.assign(amplitudes_.size(), 0.0);
        if (amplitudes_.empty() || amplitudes_.size() != frequencies_.size() || amplitudes_.size() != phases_.size())
            throw std::invalid_argument("reference_signal: waveform amplitudes, frequencies and phases must match");
        const auto& ramp = target.at("ramp");
        if (!ramp.isNull() && ramp.at("enabled").boolean(true)) {
            duration_ = ramp.at("duration").number();
            from_ = numbers(ramp.at("from")).at(0);
        }
    }

    TargetSample operator()(double t) const {
        double sine = 0.0, cosine = 0.0, curvature = 0.0;
        for (std::size_t k = 0; k < amplitudes_.size(); ++k) {
            const double angle = t * frequencies_[k] + phases_[k];
            const double S = std::sin(angle), C = std::cos(angle);
            sine += S * amplitudes_[k];
            cosine += C * (amplitudes_[k] * frequencies_[k]);
            curvature += -S * (amplitudes_[k] * (frequencies_[k] * frequencies_[k]));
        }
        const double r = offset_ + sine, v = amplitudes_.empty() ? constantVelocity_ : cosine;
        if (!(duration_ > 0.0)) return {r, v, curvature};
        const double x = std::min(1.0, std::max(0.0, t / duration_));
        const double s = std::pow(x, 3.0) * (10.0 - 15.0 * x + 6.0 * x * x);
        const bool inside = t > 0.0 && t < duration_;
        const double ds = inside ? 30.0 * x * x * ((1.0 - x) * (1.0 - x)) / duration_ : 0.0;
        const double dds = inside ? 60.0 * x * (1.0 - x) * (1.0 - 2.0 * x) / (duration_ * duration_) : 0.0;
        const double lead = r - from_;
        return {from_ + s * lead, s * v + ds * lead, s * curvature + 2.0 * ds * v + dds * lead};
    }

private:
    double offset_{}, constantVelocity_{}, duration_{}, from_{};
    std::vector<double> amplitudes_, frequencies_, phases_;
};

c::SystemConfiguration systemFrom(const json::Value& settings, const std::vector<std::filesystem::path>& sources) {
    c::SystemConfiguration out;
    if (!settings.at("system").isNull()) {
        out.identifier = settings.at("system").at("name").string();
        out.preset = settings.at("system").at("preset").string();
        for (const auto& [key, value] : settings.at("system").at("parameters").object()) out.parameters[key] = value.number();
    } else {
        if (sources.empty()) throw std::invalid_argument("reference_signal needs a system or a source archive");
        const auto metadata = json::parse(c::readNpzIndex(sources.front()).metadataJson);
        out.identifier = metadata.at("system").string();
        out.preset = metadata.at("system_preset").string();
        for (const auto& [key, value] : metadata.at("resolved_system_parameters").object()) out.parameters[key] = value.number();
    }
    if (out.identifier.empty()) throw std::invalid_argument("reference_signal could not determine the system");
    return out;
}

}  // namespace

std::filesystem::path runReferenceSignal(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                         const std::filesystem::path& root, const std::string& label) {
    const auto system = c::createSystem(systemFrom(settings, sources));
    const auto dimension = system->initialState().size();
    if (dimension != 2) throw std::invalid_argument("reference_signal currently supports position-velocity targets");
    const auto& targetSpec = settings.at("target");
    (void)resolveControlTarget(targetSpec, dimension);  // the controller's own validation of the spec
    const TargetJet jet(targetSpec);
    std::vector<double> clock;
    if (!settings.at("clock").isNull()) {
        if (sources.empty()) throw std::invalid_argument("reference_signal clock needs a source archive");
        clock = c::readNpz(sources.front()).at(settings.at("clock").at("array").string("time")).doubles();
    } else {
        const auto& time = settings.at("time");
        const double start = time.at("start").number(0.0), stop = time.at("stop").number(NAN), step = time.at("step").number(NAN);
        if (!(stop > start) || !(step > 0.0)) throw std::invalid_argument("reference_signal time needs start < stop and step > 0");
        const auto count = static_cast<std::size_t>(std::floor((stop - start) / step + 1e-9)) + 1;
        for (std::size_t k = 0; k < count; ++k) clock.push_back(start + static_cast<double>(k) * step);
    }
    const auto n = clock.size();
    std::vector<double> states(n * dimension), derivatives(n * dimension), drive(n), residual(n);
    std::vector<double> free(dimension), unit(dimension);
    for (std::size_t k = 0; k < n; ++k) {
        const auto sample = jet(clock[k]);
        const std::vector<double> x{sample.position, sample.velocity}, xdot{sample.velocity, sample.acceleration};
        system->rhsWithInput(clock[k], x.data(), 0.0, free.data());
        system->rhsWithInput(clock[k], x.data(), 1.0, unit.data());
        double gg = 0.0, gr = 0.0;
        for (std::size_t i = 0; i < dimension; ++i) {
            const double g = unit[i] - free[i];
            gg += g * g; gr += g * (xdot[i] - free[i]);
        }
        if (!(gg > 0.0)) throw std::runtime_error("reference_signal: the input has no effect at this target state");
        // Matched input (only the last state equation is driven): the plain quotient keeps its
        // rounding; a general g uses the least-squares solution.
        const double u = unit[0] == free[0] ? (xdot[1] - free[1]) / (unit[1] - free[1]) : gr / gg;
        double leftover = 0.0;
        for (std::size_t i = 0; i < dimension; ++i) {
            const double r = xdot[i] - free[i] - u * (unit[i] - free[i]);
            leftover += r * r;
        }
        drive[k] = u; residual[k] = std::sqrt(leftover);
        for (std::size_t i = 0; i < dimension; ++i) { states[k * dimension + i] = x[i]; derivatives[k * dimension + i] = xdot[i]; }
    }
    c::NpzRecord record;
    record.arrays["time"] = c::NpyArray::fromDoubles({n}, clock);
    record.arrays["target_values"] = c::NpyArray::fromDoubles({n, dimension}, states);
    record.arrays["target_derivative"] = c::NpyArray::fromDoubles({n, dimension}, derivatives);
    record.arrays["inverse_drive"] = c::NpyArray::fromDoubles({n, 1}, drive);
    record.arrays["matching_residual"] = c::NpyArray::fromDoubles({n}, residual);
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(json::Value::Object{
        {"schema", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("reference-signal")},
        {"label", json::Value(label)}, {"system", json::Value(system->descriptor().identifier)},
        {"target", targetSpec}, {"data_role", json::Value("reference_signal")}})));
    const auto directory = kindFolder(root, dynsyslab::core::object_kind::timeSeries);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    stampDataObject(record, dynsyslab::core::object_kind::timeSeries, "reference-signal", {});
    c::writeNpzAtomicNew(path, record);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
