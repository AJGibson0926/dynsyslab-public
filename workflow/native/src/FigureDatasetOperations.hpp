#pragma once
// Operations of a figure_dataset array spec ("operation"), applied to the columns the spec selects:
//   "column"      the columns as selected (default)
//   "norm"        Euclidean norm across the selected columns, one value per row
//   "difference"  these columns minus those of spec "minus" (same row count; one column broadcasts)
//   "diff"        successive differences v[k+1] - v[k] of each column (n - 1 rows); with
//                 "per_time": true divided by the matching differences of "clock" (default "time")
//   "psi"         K-ROM eigenfunctions ("krom": PATH, "modes": [k,...] default all) evaluated on
//                 every row of the source array ("array"), one column per mode; "part": "real"
//                 (default), "imag", "abs" or "phase". Rows are physical states or the K-ROM's own
//                 coordinates; a causal-history K-ROM (psi of xi) takes coordinates only, e.g. xi
//                 lifted by control_target. "domain": "state" | "coordinates" overrides the choice.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

using FigureSeries = std::vector<std::vector<double>>;  // columns

[[nodiscard]] FigureSeries psiSeries(const core::json::Value& spec, const core::NpzRecord& source);
void applyOperation(const std::string& operation, const core::json::Value& spec, FigureSeries& series,
                    const core::NpzRecord& source,
                    const std::function<FigureSeries(const core::json::Value&)>& seriesOf);

}  // namespace dynsyslab::workflow::postprocess_detail
