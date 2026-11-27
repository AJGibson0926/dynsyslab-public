#include "PrescribedControl.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>

// Designed open-loop inputs (docs/architecture/NATIVE_SOLVER.md):
//   tabulated            piecewise-linear (or held) input from times/values, ends held.
//   command_tape         exact interval-held commands, inline or from an .npy/.npz file,
//                        trajectory- or transition-major (MATLAB Ubig), with a scale factor
//                        and a SHA-256 digest of the tape.
//   telegraph_multisine  persistently exciting tape: a balanced multilevel random telegraph
//                        blended with a random-phase multisine, zero mean, peak clipped
//                        (Ljung, System Identification, 2nd ed. 1999, Sec. 13.3; Pintelon &
//                        Schoukens, System Identification, 2012, Ch. 2).
//   symmetric_pulse      +/- pulse pairs: trajectory 2k is -magnitude, 2k+1 is +magnitude, so
//                        paired first responses isolate the input direction (pair_id, polarity).
//   pulse_branches       held pulses of several amplitudes and widths from shared anchors
//                        (anchor, branch; +/- branches pair as pair_id, polarity).
// Command tapes for Duffing follow Korda & Mezic, Automatica 93 (2018) 149-160.

namespace dynsyslab::workflow::detail {
namespace {

std::vector<double> numbers(const json::Value& source) {
    std::vector<double> result;
    for (const auto& item : source.array()) {
        const double value = item.number(std::numeric_limits<double>::quiet_NaN());
        if (!std::isfinite(value)) throw std::invalid_argument("prescribed input values must be finite");
        result.push_back(value);
    }
    return result;
}

// Rows of a 1-D (shared) or 2-D (per trajectory) JSON array.
std::vector<std::vector<double>> rows(const json::Value& source) {
    const auto& items = source.array();
    if (!items.empty() && !items.front().array().empty()) {
        std::vector<std::vector<double>> result;
        for (const auto& row : items) result.push_back(numbers(row));
        return result;
    }
    return {numbers(source)};
}

std::uint64_t seedOf(const json::Value& source) {
    const double value = source.number();
    if (!std::isfinite(value) || value < 0.0 || std::floor(value) != value || value > 9007199254740991.0)
        throw std::invalid_argument("input seed must be a nonnegative integer");
    return static_cast<std::uint64_t>(value);
}

bool physicalUnits(const json::Value& settings, const char* key) {
    const auto units = settings.at(key).string("model");
    if (units != "model" && units != "physical")
        throw std::invalid_argument(std::string(key) + " must be model or physical");
    return units == "physical";
}

std::string tapeDigest(const std::vector<std::vector<double>>& tape) {
    std::string bytes;
    for (const auto& row : tape)
        for (const double value : row) {
            char raw[sizeof(double)];
            std::memcpy(raw, &value, sizeof(double));
            bytes.append(raw, sizeof(double));
        }
    return core::sha256Text(bytes);
}

const std::vector<double>& rowFor(const std::vector<std::vector<double>>& table, std::size_t trajectory) {
    if (table.size() == 1) return table.front();
    if (trajectory >= table.size())
        throw std::out_of_range("prescribed input has no row for trajectory " + std::to_string(trajectory));
    return table[trajectory];
}

// Interval-held tape law: command k on [start + k dt, start + (k + 1) dt), last held after.
core::ScalarControlLaw heldTape(std::shared_ptr<const std::vector<std::vector<double>>> tape,
                                double start, double interval) {
    return [tape = std::move(tape), start, interval](std::size_t trajectory, double time, const double*) {
        const auto& row = rowFor(*tape, trajectory);
        const double position = (time - start) / interval;
        // Times within roundoff of a boundary belong to the interval that starts there.
        const double index = std::floor(position + 1.0e-9);
        const auto k = index < 0.0 ? std::size_t{0}
            : std::min(static_cast<std::size_t>(index), row.size() - 1);
        return row[k];
    };
}

PrescribedControl tabulated(const json::Value& settings, const core::CoordinateScaling& scaling) {
    auto times = numbers(settings.at("times"));
    auto table = rows(settings.at("values"));
    if (times.size() < 2) throw std::invalid_argument("tabulated input needs at least two times");
    for (std::size_t index = 1; index < times.size(); ++index)
        if (!(times[index] > times[index - 1])) throw std::invalid_argument("tabulated input times must increase");
    for (const auto& row : table)
        if (row.size() != times.size()) throw std::invalid_argument("tabulated input values must match its times");
    const bool physicalTime = physicalUnits(settings, "time_units");
    const bool physicalInput = physicalUnits(settings, "input_units");
    if (physicalTime) for (auto& time : times) time = scaling.modelTime(time);
    if (physicalInput) for (auto& row : table) for (auto& value : row) value = scaling.modelInput(value);
    const auto mode = settings.at("interpolation").string("linear");
    if (mode != "linear" && mode != "previous")
        throw std::invalid_argument("tabulated input interpolation must be linear or previous");
    PrescribedControl result;
    result.law = [times = std::move(times), table = std::move(table), held = mode == "previous"]
        (std::size_t trajectory, double time, const double*) {
        const auto& row = rowFor(table, trajectory);
        if (time <= times.front()) return row.front();
        if (time >= times.back()) return row.back();
        const auto upper = static_cast<std::size_t>(std::upper_bound(times.begin(), times.end(), time) - times.begin());
        const double t0 = times[upper - 1], t1 = times[upper];
        if (held) return row[upper - 1];
        return row[upper - 1] + (time - t0) / (t1 - t0) * (row[upper] - row[upper - 1]);
    };
    result.settingsCoordinates = physicalTime == physicalInput ? (physicalTime ? "physical" : "model") : "mixed";
    return result;
}

std::vector<std::vector<double>> tapeFromFile(const json::Value& settings) {
    const auto path = std::filesystem::u8path(settings.at("file").string());
    const auto record = core::readNpz(path);
    const auto name = settings.at("array").string(record.arrays.size() == 1 ? record.arrays.begin()->first : "commands");
    const auto& array = record.at(name);
    const auto values = array.doubles();
    if (array.shape.empty() || array.shape.size() > 3)
        throw std::invalid_argument("command tape array must be 1-D, 2-D or (trajectory, transition, 1)");
    const std::size_t first = array.shape[0], second = array.shape.size() > 1 ? array.shape[1] : 1;
    if (array.shape.size() == 1) return {values};
    std::vector<std::vector<double>> result(first, std::vector<double>(second));
    for (std::size_t r = 0; r < first; ++r)
        for (std::size_t c = 0; c < second; ++c) result[r][c] = values[r * second + c];
    return result;
}

PrescribedControl commandTape(const json::Value& settings, const core::CoordinateScaling& scaling) {
    auto tape = settings.at("file").string().empty() ? rows(settings.at("values")) : tapeFromFile(settings);
    const auto layout = settings.at("layout").string("trajectory_major");
    if (layout == "transition_major") {
        // MATLAB Ubig layout: rows are transitions, columns trajectories.
        std::vector<std::vector<double>> transposed(tape.front().size(), std::vector<double>(tape.size()));
        for (std::size_t r = 0; r < tape.size(); ++r)
            for (std::size_t c = 0; c < tape[r].size(); ++c) transposed[c][r] = tape[r][c];
        tape = std::move(transposed);
    } else if (layout != "trajectory_major") {
        throw std::invalid_argument("command tape layout must be trajectory_major or transition_major");
    }
    if (tape.empty() || tape.front().empty()) throw std::invalid_argument("command tape is empty");
    for (const auto& row : tape)
        if (row.size() != tape.front().size()) throw std::invalid_argument("command tape rows differ in length");
    const auto digest = tapeDigest(tape);
    const double scale = settings.at("scale").number(1.0);
    const bool physicalInput = physicalUnits(settings, "input_units");
    for (auto& row : tape)
        for (auto& value : row) {
            value *= scale;
            if (physicalInput) value = scaling.modelInput(value);
            if (!std::isfinite(value)) throw std::invalid_argument("command tape values must be finite");
        }
    double interval = settings.at("interval").number(0.0);
    double start = settings.at("start").number(0.0);
    const bool physicalTime = physicalUnits(settings, "time_units");
    if (physicalTime) { interval = scaling.modelTime(interval); start = scaling.modelTime(start); }
    if (!(interval > 0.0)) throw std::invalid_argument("command tape needs a positive interval");
    PrescribedControl result;
    const auto shared = std::make_shared<const std::vector<std::vector<double>>>(std::move(tape));
    result.law = heldTape(shared, start, interval);
    result.trajectoryMetadata = [digest, scale, interval, shared](std::size_t trajectory) {
        const auto& row = rowFor(*shared, trajectory);
        double peak = 0.0;
        for (const double value : row) peak = std::max(peak, std::abs(value));
        return json::Value(json::Value::Object{{"kind", std::string("command_tape")},
            {"tape_sha256", digest}, {"scale", scale}, {"interval", interval},
            {"transitions", static_cast<double>(row.size())}, {"command_peak", peak},
            {"command_tape_exact", true}});
    };
    result.settingsCoordinates = physicalTime == physicalInput ? (physicalTime ? "physical" : "model") : "mixed";
    return result;
}

// One telegraph-plus-multisine tape for trajectory `member` (Python generator semantics with a
// per-trajectory std::mt19937_64 stream seeded by (seed, member)).
struct TelegraphDesign {
    std::uint64_t seed{};
    std::size_t transitions{};
    double interval{};
    std::vector<double> amplitudes, frequencies, levels;
    std::vector<std::size_t> holds;
    double fraction{};

    std::vector<double> tape(std::size_t member, std::size_t& hold, double& amplitude) const {
        std::seed_seq sequence{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32),
                               static_cast<std::uint32_t>(member), static_cast<std::uint32_t>(member >> 32)};
        std::mt19937_64 engine(sequence);
        amplitude = amplitudes[member % amplitudes.size()];
        hold = holds[(member / amplitudes.size()) % holds.size()];
        const std::size_t blocks = (transitions + hold - 1) / hold;
        std::vector<double> blockLevels(blocks);
        for (std::size_t block = 0; block < blocks; ++block) blockLevels[block] = levels[block % levels.size()];
        std::shuffle(blockLevels.begin(), blockLevels.end(), engine);
        std::uniform_real_distribution<double> phase(0.0, 6.283185307179586476925286766559);
        std::vector<double> phases(frequencies.size());
        for (auto& value : phases) value = phase(engine);
        std::vector<double> multisine(transitions, 0.0), command(transitions);
        double peakSine = 0.0;
        for (std::size_t k = 0; k < transitions; ++k) {
            const double time = (static_cast<double>(k) + 0.5) * interval;
            for (std::size_t f = 0; f < frequencies.size(); ++f)
                multisine[k] += std::sin(time * frequencies[f] + phases[f]);
            peakSine = std::max(peakSine, std::abs(multisine[k]));
        }
        peakSine = std::max(peakSine, std::numeric_limits<double>::min());
        double mean = 0.0;
        for (std::size_t k = 0; k < transitions; ++k) {
            command[k] = amplitude * ((1.0 - fraction) * blockLevels[k / hold] + fraction * multisine[k] / peakSine);
            mean += command[k];
        }
        mean /= static_cast<double>(transitions);
        double peak = 0.0;
        for (auto& value : command) { value -= mean; peak = std::max(peak, std::abs(value)); }
        if (peak > amplitude) for (auto& value : command) value *= amplitude / peak;
        return command;
    }
};

PrescribedControl telegraphMultisine(const json::Value& settings) {
    TelegraphDesign design;
    design.seed = seedOf(settings.at("seed"));
    const double transitions = settings.at("transition_count").number(0.0);
    if (!(transitions >= 1.0) || std::floor(transitions) != transitions)
        throw std::invalid_argument("telegraph_multisine needs a positive integer transition_count");
    design.transitions = static_cast<std::size_t>(transitions);
    design.interval = settings.at("interval").number(0.0);
    design.amplitudes = numbers(settings.at("amplitudes"));
    design.frequencies = numbers(settings.at("angular_frequencies"));
    design.levels = numbers(settings.at("levels"));
    if (design.levels.empty()) design.levels = {-1.0, -0.5, 0.5, 1.0};
    for (const double hold : numbers(settings.at("hold_steps"))) {
        if (!(hold >= 1.0) || std::floor(hold) != hold)
            throw std::invalid_argument("telegraph_multisine hold_steps must be positive integers");
        design.holds.push_back(static_cast<std::size_t>(hold));
    }
    design.fraction = settings.at("multisine_fraction").number(0.0);
    if (!(design.interval > 0.0) || design.amplitudes.empty() || design.holds.empty()
        || std::any_of(design.amplitudes.begin(), design.amplitudes.end(), [](double a) { return !(a > 0.0); })
        || !(design.fraction >= 0.0 && design.fraction <= 1.0)
        || (design.fraction > 0.0 && design.frequencies.empty()))
        throw std::invalid_argument("telegraph_multisine needs interval > 0, positive amplitudes, hold_steps, "
                                    "multisine_fraction in [0, 1] and angular_frequencies for the multisine");
    const double start = settings.at("start").number(0.0);
    auto cache = std::make_shared<std::map<std::size_t, std::vector<double>>>();
    auto shared = std::make_shared<const TelegraphDesign>(std::move(design));
    const auto tapeFor = [cache, shared](std::size_t member) -> const std::vector<double>& {
        auto found = cache->find(member);
        if (found == cache->end()) {
            std::size_t hold{}; double amplitude{};
            found = cache->emplace(member, shared->tape(member, hold, amplitude)).first;
        }
        return found->second;
    };
    PrescribedControl result;
    result.law = [tapeFor, shared, start](std::size_t trajectory, double time, const double*) {
        const auto& row = tapeFor(trajectory);
        const double index = std::floor((time - start) / shared->interval + 1.0e-9);
        const auto k = index < 0.0 ? std::size_t{0} : std::min(static_cast<std::size_t>(index), row.size() - 1);
        return row[k];
    };
    result.trajectoryMetadata = [tapeFor, shared](std::size_t trajectory) {
        std::size_t hold{}; double amplitude{};
        (void)shared->tape(trajectory, hold, amplitude);
        const auto& row = tapeFor(trajectory);
        double peak = 0.0;
        for (const double value : row) peak = std::max(peak, std::abs(value));
        return json::Value(json::Value::Object{{"kind", std::string("telegraph_multisine")},
            {"tape_sha256", tapeDigest({row})}, {"hold_steps", static_cast<double>(hold)},
            {"amplitude", amplitude}, {"command_peak", peak}, {"command_tape_exact", true}});
    };
    return result;
}

PrescribedControl symmetricPulse(const json::Value& settings, const core::CoordinateScaling& scaling) {
    const bool physicalInput = physicalUnits(settings, "input_units");
    double magnitude = settings.at("magnitude").number(1.0);
    if (physicalInput) magnitude = scaling.modelInput(magnitude);
    if (!(magnitude > 0.0) || !std::isfinite(magnitude))
        throw std::invalid_argument("symmetric_pulse magnitude must be positive");
    PrescribedControl result;
    result.law = [magnitude](std::size_t trajectory, double, const double*) {
        return trajectory % 2 == 0 ? -magnitude : magnitude;
    };
    result.trajectoryMetadata = [magnitude](std::size_t trajectory) {
        return json::Value(json::Value::Object{{"kind", std::string("symmetric_pulse")},
            {"pair_id", static_cast<double>(trajectory / 2)},
            {"polarity", trajectory % 2 == 0 ? -1.0 : 1.0},
            {"command", trajectory % 2 == 0 ? -magnitude : magnitude},
            {"command_tape_exact", true}});
    };
    result.settingsCoordinates = physicalInput ? "physical" : "model";
    return result;
}

// Pulse branches from shared anchors: trajectory k is branch k % B of anchor k / B and holds
// amplitudes[b] for t - start < widths[b], then zero; solve with solver.control_evaluation
// "sampled" so each pulse is held exactly over whole samples (widths multiples of the sample time). Branches +a and -a of equal width form a
// symmetric pair (pair_id, polarity), so input_map fit "symmetric_pulse" reads their paired
// first responses and pulse_response their causal sensitivities.
PrescribedControl pulseBranches(const json::Value& settings, const core::CoordinateScaling& scaling) {
    const bool physicalInput = physicalUnits(settings, "input_units");
    auto amplitudes = numbers(settings.at("amplitudes"));
    if (amplitudes.empty()) throw std::invalid_argument("pulse_branches needs amplitudes");
    if (physicalInput) for (auto& a : amplitudes) a = scaling.modelInput(a);
    const auto B = amplitudes.size();
    std::vector<double> widths = std::holds_alternative<double>(settings.at("widths").storage())
        ? std::vector<double>(B, settings.at("widths").number()) : numbers(settings.at("widths"));
    if (widths.size() != B) throw std::invalid_argument("pulse_branches widths: one value or one per branch");
    for (const double w : widths) if (!(w > 0) || !std::isfinite(w)) throw std::invalid_argument("pulse_branches widths must be positive");
    const double start = settings.at("start").number(0.0);
    std::vector<long long> partner(B, -1);
    for (std::size_t b = 0; b < B; ++b)
        for (std::size_t c = 0; c < B && amplitudes[b] != 0; ++c)
            if (c != b && amplitudes[c] == -amplitudes[b] && widths[c] == widths[b]) { partner[b] = (long long)c; break; }
    PrescribedControl result;
    result.law = [amplitudes, widths, start, B](std::size_t trajectory, double t, const double*) {
        const auto b = trajectory % B;
        return t - start < widths[b] - 1e-12 * std::max(1.0, std::abs(widths[b])) ? amplitudes[b] : 0.0;
    };
    result.trajectoryMetadata = [amplitudes, widths, partner, B](std::size_t trajectory) {
        const auto b = trajectory % B, anchor = trajectory / B;
        json::Value::Object meta{{"kind", std::string("pulse_branches")}, {"anchor", double(anchor)}, {"branch", double(b)},
            {"command", amplitudes[b]}, {"width", widths[b]}, {"branches", double(B)}, {"command_tape_exact", true}};
        if (partner[b] >= 0) {
            meta["pair_id"] = double(anchor * B + std::min<std::size_t>(b, std::size_t(partner[b])));
            meta["polarity"] = amplitudes[b] > 0 ? 1.0 : -1.0;
        }
        return json::Value(std::move(meta));
    };
    result.settingsCoordinates = physicalInput ? "physical" : "model";
    return result;
}

}  // namespace

PrescribedControl buildDesignedInput(const std::string& kind, const json::Value& settings,
                                     const core::CoordinateScaling& scaling) {
    if (kind == "tabulated") return tabulated(settings, scaling);
    if (kind == "command_tape") return commandTape(settings, scaling);
    if (kind == "telegraph_multisine") return telegraphMultisine(settings);
    if (kind == "symmetric_pulse") return symmetricPulse(settings, scaling);
    if (kind == "pulse_branches") return pulseBranches(settings, scaling);
    throw std::invalid_argument("unknown prescribed input: " + kind);
}

}  // namespace dynsyslab::workflow::detail
