#pragma once
// Held-input pairs (psi_k, u_k) -> psi_{k+1} from forced records, for every coordinate
// kind: state/observable KROMs, measured coordinates, and causal-history KROMs (through
// the frozen transform). Shared by input_map and bilinear_certificate.
//
// A pair is usable only when u really was constant across the interval:
//   control_interval_values            interval-aligned, always held;
//   sampled control, clock = N * dt    held across each output interval;
//   sampled control, other clock       refused when strict;
//   other sample-aligned records       strict keeps only pairs with u_{k+1} == u_k.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <complex>
#include <cstddef>
#include <filesystem>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

struct HeldInputRow {
    std::vector<std::complex<double>> psi, next;
    double u{};
};

struct HeldInputRecords {
    std::vector<HeldInputRow> rows;
    double sampleTime{};
    std::size_t droppedAtInputChange{};  // strict mode: pairs whose input was not held
};

// `transform` is required for causal-history KROMs and ignored otherwise. Reads
// fit_start, measurement_observable and (not for history) state_window from settings.
HeldInputRecords heldInputRecords(const dynsyslab::core::KoopmanModel& model,
    const dynsyslab::core::FrozenHistoryTransform* transform,
    const std::vector<std::filesystem::path>& sources, const json::Value& settings, bool strictHold);

}  // namespace dynsyslab::workflow::postprocess_detail
