#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>

#include <optional>

namespace dsl {

// Draws the legend at `origin` (top-left corner), or at the default top-right corner.
void drawPlotLegend(Gdiplus::Graphics& graphics, const Figure& figure,
                    const Gdiplus::RectF& plot, const theme::Fonts& fonts,
                    const Gdiplus::PointF* origin = nullptr);
Gdiplus::SizeF plotLegendSize(Gdiplus::Graphics& graphics, const Figure& figure,
                              const Gdiplus::RectF& plot, const theme::Fonts& fonts);
// legendPlacement "top": the entries side by side in rows across `width`, centred. The height is
// what plotArea reserves under the title; zero when the figure has no labelled traces.
float plotTopLegendHeight(Gdiplus::Graphics& graphics, const Figure& figure, float width,
                          const theme::Fonts& fonts);
// Titles the card does not draw because the panel already says them.
bool redundantPlotTitle(const Figure& figure);
// The band under the title that holds the legend as one row (std::nullopt: the legend goes inside
// the plot). `plot` is the plot laid out without a band. Figure::legendPlacement "top" always gets
// the band; the automatic placement of time and spectrum plots gets it only when no corner of the
// plot is clear of data (clearLegendCorner).
std::optional<Gdiplus::RectF> topLegendBand(const Figure& figure, const Gdiplus::RectF& card,
                                            const Gdiplus::RectF& plot, const theme::Fonts& fonts);
// Time and spectrum plots: the first corner of `plot` (top right, top left, bottom right, bottom
// left) where the overlay legend covers no sample of any trace; std::nullopt when every corner
// covers some, or for other kinds of plot.
std::optional<Gdiplus::RectF> clearLegendCorner(const Figure& figure, const Gdiplus::RectF& plot,
                                                const theme::Fonts& fonts);
void drawPlotTopLegend(Gdiplus::Graphics& graphics, const Figure& figure,
                       const Gdiplus::RectF& band, const theme::Fonts& fonts);

}  // namespace dsl
