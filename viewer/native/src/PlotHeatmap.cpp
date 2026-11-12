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
void drawScalarColorbar(Gdiplus::Graphics& graphics, const Figure& figure,
    const Gdiplus::RectF& plot, const theme::Fonts& fonts) {
    if (figure.traces.empty() || (figure.kind != "observable_heatmap" && figure.kind != "krom_heatmap"
        && figure.traces.front().encoding != "domain-color")) return;
    const bool domain = figure.traces.front().encoding == "domain-color";
    const auto [low,high] = colorRange(figure.traces.front(),resolveValueScale(figure));
    const auto valueScale = resolveValueScale(figure);
    const float scale = scientificTextScale(fonts.body);
    const float width = std::min(180.f*scale,plot.Width*.45f), x = plot.GetRight()-width-64*scale, y = plot.Y-20*scale;
    for (int i=0;i<128;++i) {
        Gdiplus::SolidBrush brush(domain ? domainColor(low+(high-low)*i/127.) : theme::heatColor(i/127.));
        graphics.FillRectangle(&brush,x+width*i/128,y,width/128+.5f,6.f*scale);
    }
    if (!domain && valueScale.kind != ValueScaleKind::linear) {
        const auto label = valueScale.kind == ValueScaleKind::symmetric_log ? L"Symmetric log" : L"Logarithmic";
        drawMathText(graphics,label,{x,y+5*scale,width,19*scale},fonts.caption,theme::plotMuted,Gdiplus::StringAlignmentCenter);
    }
    wchar_t lower[32]{},upper[32]{};
    swprintf_s(lower,L"%.3g",low); swprintf_s(upper,L"%.3g",high);
    drawMathText(graphics,lower,{x-68*scale,y-9*scale,64*scale,22*scale},fonts.caption,theme::plotMuted,Gdiplus::StringAlignmentFar);
    drawMathText(graphics,upper,{x+width+4*scale,y-9*scale,64*scale,22*scale},fonts.caption,theme::plotMuted);
}
}
