#pragma once
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Theme.hpp"
#include <gdiplus.h>
namespace dsl {
inline constexpr float curveLineWidthPoints = 3.6f;
inline constexpr float phaseLineWidthPoints = 2.5f;
inline constexpr float targetLineWidthPoints = 3.0f;
inline constexpr Gdiplus::REAL targetDashPattern[]{4.f, 6.f};
float plotSeriesWidthPoints(const Trace&);
void stylePlotSeriesPen(Gdiplus::Pen&, const Trace&);
Gdiplus::Color applyPlotSeriesOpacity(Gdiplus::Color, const Trace&);
Gdiplus::Color plotSeriesColor(const Trace&, std::size_t traceIndex);
Gdiplus::Color trajectoryTimeColor(double fraction);
std::pair<double,double> trajectoryTimeRange(const Figure&);
Gdiplus::RectF plotReferenceLabelBounds(const Gdiplus::PointF&,float markerRadius,
    const Gdiplus::SizeF&,const Gdiplus::RectF& plotArea);
void drawPlotSeries(Gdiplus::Graphics&,const Trace&,std::size_t,const PlotTransform&,
    std::optional<std::pair<double,double>> timeRange = std::nullopt);
void drawPlotReferences(Gdiplus::Graphics&,const Figure&,const PlotTransform&,const theme::Fonts&);
void drawTimeColorbar(Gdiplus::Graphics&,const Figure&,const Gdiplus::RectF&,const theme::Fonts&);
}
