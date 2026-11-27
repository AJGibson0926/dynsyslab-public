#include "PlotHeatmap.hpp"
#include "MathTextDrawing.hpp"
#include "DomainColor.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dsl {
namespace {
std::pair<double,double> colorRange(const Trace& trace, const ValueScale& scale = {}) {
    if (trace.encoding == "domain-color") return {-3.141592653589793,3.141592653589793};
    if (trace.colorLimits && trace.colorLimits->second > trace.colorLimits->first) return *trace.colorLimits;
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (const auto value : trace.z) if (std::isfinite(value)) {
        low = std::min(low,value); high = std::max(high,value);
    }
    if (!std::isfinite(low)) return {-1,1};
    if (low == high) return {scale.inverse(scale.forward(low)-.5),scale.inverse(scale.forward(high)+.5)};
    return {low,high};
}
}
void drawScalarHeatmap(Gdiplus::Graphics& graphics, const Trace& trace, const PlotTransform& t, const ValueScale& scale) {
    if (trace.zShape.size() < 2 || trace.x.size() < 2 || trace.y.size() < 2) return;
    const auto rows = trace.zShape[trace.zShape.size()-2], columns = trace.zShape.back();
    if (rows < 2 || columns < 2 || trace.z.size() < static_cast<std::size_t>(rows*columns)) return;
    const auto [low,high] = colorRange(trace,scale);
    const auto lower = scale.forward(low), upper = scale.forward(high);
    Gdiplus::Bitmap bitmap(columns,rows,PixelFormat32bppARGB);
    for (int row=0;row<rows;++row) for (int column=0;column<columns;++column)
        bitmap.SetPixel(column,rows-1-row,trace.encoding == "domain-color"
            ? domainColor(trace.z[row*columns+column],static_cast<std::size_t>(row*columns+column) < trace.magnitudes.size() ? trace.magnitudes[row*columns+column] : 1.)
            : !std::isfinite(trace.z[row*columns+column]) ? theme::color(theme::plotBackground)
            : theme::heatColor((scale.forward(trace.z[row*columns+column])-lower)/(upper-lower)));
    const auto state = graphics.Save();
    const Gdiplus::RectF area(static_cast<float>(t.left),static_cast<float>(t.top),
        static_cast<float>(t.width),static_cast<float>(t.height));
    graphics.SetClip(area,Gdiplus::CombineModeIntersect);
    graphics.SetInterpolationMode(Gdiplus::InterpolationModeHighQualityBilinear);
    const float left = static_cast<float>(t.left+(trace.x.front()-t.xMin)*t.width/(t.xMax-t.xMin));
    const float top = static_cast<float>(t.top+(t.yMax-trace.y.back())*t.height/(t.yMax-t.yMin));
    const float width = static_cast<float>((trace.x.back()-trace.x.front())*t.width/(t.xMax-t.xMin));
    const float height = static_cast<float>((trace.y.back()-trace.y.front())*t.height/(t.yMax-t.yMin));
    Gdiplus::ImageAttributes attributes; attributes.SetWrapMode(Gdiplus::WrapModeTileFlipXY);
    graphics.DrawImage(&bitmap,Gdiplus::RectF(left,top,width,height),0.f,0.f,static_cast<float>(columns),
        static_cast<float>(rows),Gdiplus::UnitPixel,&attributes);
    graphics.Restore(state);
}
std::pair<double,double> heatmapColorRange(const Trace& trace, const ValueScale& scale) {
    return colorRange(trace,scale);
}
}
