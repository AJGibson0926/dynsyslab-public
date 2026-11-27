#pragma once

#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

// {"kind":"viewer_figure","source":"${inputs.table}","output":"name.png",
//  "axis_limits":{"reference":{"source":"${inputs.other}","axis":0,"figure":{...}}},
//  "figure":{"view":"time-series","width":1400,"height":800,"columns":1,
//            "style":{"font_scale":1.8,"line_width_scale":1.8,"default_line_width":3.6},
//            "start":0,"stop":20,"x_range":[a,b],"y_range":[a,b],"title":"...","time_label":"...",
//            "x_label":"...","value_label":"...","axes":["..."],
//            "curves":[{"series":"r,v","label":"...","color":"#202020","line_width":5.4,
//                       "line_style":"solid"}]}}
[[nodiscard]] std::vector<std::string> viewerFigureArguments(const json::Value& job,
                                                             const std::filesystem::path& output);
// Curves with their own "source" and "clock", and "figure.clock" (ViewerFigureSources.cpp):
//   "curves":[{"series":"values","source":"${inputs.other}","clock":"estimator_time",
//               "interpolation":"linear"}], "clock":{"source":"${inputs.table}","array":"time"}
// Every curve is resampled onto the figure clock into one scratch archive before rendering.
[[nodiscard]] bool hasCurveSources(const json::Value& job);
[[nodiscard]] json::Value assembleCurveSources(const json::Value& job, const std::filesystem::path& scratch);
// "figure.fit_y_to":[curve indices]: the ordinate range of a model-scale page (ViewerFigureSources.cpp).
[[nodiscard]] json::Value fitOrdinate(const json::Value& job);
[[nodiscard]] json::Value planViewerFigure(const json::Value& job);
[[nodiscard]] RunResult executeViewerFigure(const json::Value& job, const std::filesystem::path& archiveRoot);

}  // namespace dynsyslab::workflow::detail
