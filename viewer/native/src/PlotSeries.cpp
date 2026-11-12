#include "PlotSeries.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "MathTextDrawing.hpp"
#include "DomainColor.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dsl {
namespace {
Gdiplus::PointF project(double x,double y,const PlotTransform& t) {
    return {static_cast<float>(t.left+(x-t.xMin)/(t.xMax-t.xMin)*t.width),
        static_cast<float>(plotPixelY(t,y))};
}
std::pair<double,double> timeRange(const std::vector<double>& values) {
    double low=std::numeric_limits<double>::infinity(),high=-low;
    for(const auto value:values) if(std::isfinite(value)) { low=std::min(low,value);high=std::max(high,value); }
    return std::isfinite(low)?std::make_pair(low,high):std::make_pair(0.,1.);
}
}
float plotSeriesWidthPoints(const Trace& trace) {
    if (std::isfinite(trace.lineWidthPoints) && trace.lineWidthPoints>0.)
        return static_cast<float>(std::min(32.,trace.lineWidthPoints));
    if (trace.role == "phase_vector_field") return .65f;
    if (trace.role.rfind("phase_",0)==0) return 1.05f;
    return trace.role == "control_target" ? targetLineWidthPoints
        : trace.encoding == "time" ? phaseLineWidthPoints : curveLineWidthPoints;
}
void stylePlotSeriesPen(Gdiplus::Pen& pen, const Trace& trace) {
    pen.SetWidth(plotSeriesWidthPoints(trace)*96.f/72.f);
    pen.SetLineJoin(Gdiplus::LineJoinRound);
    pen.SetStartCap(Gdiplus::LineCapRound); pen.SetEndCap(Gdiplus::LineCapRound);
    if (trace.lineStyle == "dash" || trace.lineStyle == "time-dash") {
        if (trace.role == "control_target") pen.SetDashPattern(targetDashPattern,2);
        else pen.SetDashStyle(Gdiplus::DashStyleDash);
    }
}
Gdiplus::Color applyPlotSeriesOpacity(Gdiplus::Color color, const Trace& trace) {
    const auto alpha = static_cast<BYTE>(std::lround(
        static_cast<double>(color.GetA()) * std::clamp(trace.opacity, 0.0, 1.0)));
    return {alpha, color.GetR(), color.GetG(), color.GetB()};
}
Gdiplus::Color plotSeriesColor(const Trace& trace, std::size_t traceIndex) {
    return applyPlotSeriesOpacity(theme::plotColor(trace.color, traceIndex), trace);
}
Gdiplus::Color trajectoryTimeColor(double fraction) {
    return theme::heatColor(fraction);
}
std::pair<double,double> trajectoryTimeRange(const Figure& figure) {
    double low=std::numeric_limits<double>::infinity(),high=-low;
    for(const auto& trace:figure.traces) if(trace.encoding=="time") {
        const auto limits=timeRange(trace.colorValues);
        low=std::min(low,limits.first);high=std::max(high,limits.second);
    }
    return std::isfinite(low)?std::make_pair(low,high):std::make_pair(0.,1.);
}
Gdiplus::RectF plotReferenceLabelBounds(const Gdiplus::PointF& point,float markerRadius,
    const Gdiplus::SizeF& measured,const Gdiplus::RectF& area) {
    const float width=std::min(area.Width,std::max(1.f,std::ceil(measured.Width)+2.f));
    const float height=std::min(area.Height,std::max(1.f,std::ceil(measured.Height)));
    constexpr float gap=3.f;
    float x=point.X+markerRadius+gap;
    if(x+width>area.GetRight()) x=point.X-markerRadius-gap-width;
    x=std::clamp(x,area.X,std::max(area.X,area.GetRight()-width));
    const float y=std::clamp(point.Y-height*.5f,area.Y,
        std::max(area.Y,area.GetBottom()-height));
    return {x,y,width,height};
}
void drawPlotSeries(Gdiplus::Graphics& g,const Trace& trace,std::size_t traceIndex,const PlotTransform& t,
    std::optional<std::pair<double,double>> sharedTimeRange) {
    const auto count=std::min(trace.x.size(),trace.y.size()); if(!count) return;
    const auto state=g.Save();
    g.SetClip(Gdiplus::RectF(static_cast<float>(t.left),static_cast<float>(t.top),
        static_cast<float>(t.width),static_cast<float>(t.height)),Gdiplus::CombineModeIntersect);
    Gdiplus::Pen pen(plotSeriesColor(trace,traceIndex));
    stylePlotSeriesPen(pen,trace);
    if (trace.lineStyle=="time-dash") {
        // Dash along the time axis so a dense curve leaves the underlying trace visible.
        pen.SetDashStyle(Gdiplus::DashStyleSolid);
        const float dash=8.f*g.GetDpiX()/96.f;
        Gdiplus::Region stripes;
        stripes.MakeEmpty();
        for (float x=0; x<t.width; x+=2*dash)
            stripes.Union(Gdiplus::RectF(static_cast<float>(t.left)+x,static_cast<float>(t.top),dash,static_cast<float>(t.height)));
        g.SetClip(&stripes,Gdiplus::CombineModeIntersect);
    }
    const auto times=sharedTimeRange.value_or(timeRange(trace.colorValues));
    const double maxMagnitude=trace.magnitudes.empty()?1.:
        *std::max_element(trace.magnitudes.begin(),trace.magnitudes.end());
    std::vector<Gdiplus::PointF> points;
    auto flush=[&](){if(points.size()>1) g.DrawLines(&pen,points.data(),static_cast<INT>(points.size()));points.clear();};
    std::optional<Gdiplus::PointF> previous;
    for(std::size_t i=0;i<count;++i) {
        double x=trace.x[i],y=trace.y[i];
        if(i<trace.z.size()) {x-=.32*trace.z[i];y+=.20*trace.z[i];}
        if(!std::isfinite(x)||!std::isfinite(y)) {flush();previous.reset();continue;}
        const auto point=project(x,y,t);
        if(trace.colorValues.empty()) points.push_back(point);
        else if(previous) {
            const auto colorIndex=trace.encoding=="time"?i-1:i;
            const auto value=colorIndex<trace.colorValues.size()?trace.colorValues[colorIndex]:0.;
            if(trace.encoding=="time") pen.SetColor(applyPlotSeriesOpacity(trajectoryTimeColor(
                times.second>times.first?(value-times.first)/(times.second-times.first):0.),trace));
            else if(trace.encoding=="domain-color") pen.SetColor(applyPlotSeriesOpacity(
                domainColor(value,i<trace.magnitudes.size()?trace.magnitudes[i]:1.),trace));
            else pen.SetColor(applyPlotSeriesOpacity(theme::phaseColor(value,
                i<trace.magnitudes.size()&&maxMagnitude>0?trace.magnitudes[i]/maxMagnitude:1.),trace));
            g.DrawLine(&pen,*previous,point);
        }
        if(trace.marker||count==1) {
            Gdiplus::SolidBrush brush(plotSeriesColor(trace,traceIndex));
            if (trace.selectionId >= 0) {
                const float radius = trace.selected ? 5.f : 3.5f;
                if (trace.selected) g.FillEllipse(&brush,point.X-radius,point.Y-radius,radius*2,radius*2);
                Gdiplus::Pen ring(trace.selected ? theme::color(theme::plotText) : theme::color(theme::plotMuted),1.2f);
                g.DrawEllipse(&ring,point.X-radius,point.Y-radius,radius*2,radius*2);
            } else g.FillEllipse(&brush,point.X-2.2f,point.Y-2.2f,4.4f,4.4f);
        }
        previous=point;
    }
    flush();g.Restore(state);
}
void drawPlotReferences(Gdiplus::Graphics& g,const Figure& f,const PlotTransform& t,const theme::Fonts& fonts) {
    const auto state=g.Save();
    const Gdiplus::RectF area(static_cast<float>(t.left),static_cast<float>(t.top),
        static_cast<float>(t.width),static_cast<float>(t.height));
    g.SetClip(area,Gdiplus::CombineModeIntersect);
    Gdiplus::Pen pen(Gdiplus::Color(255, 255, 0, 0),1.5f);
    Gdiplus::SolidBrush red(Gdiplus::Color(255, 255, 0, 0));
    Gdiplus::SolidBrush background(theme::color(theme::plotBackground));
    for(const auto& ref:f.references) {
        if (!plotReferenceVisible(f,ref)) continue;
        const auto p=project(ref.x,ref.y,t);
        const bool coincident=std::any_of(f.references.begin(),f.references.end(),[&](const auto& other) {
            if(!plotReferenceVisible(f,other)||other.equilibrium==ref.equilibrium) return false;
            const auto q=project(other.x,other.y,t);
            return std::hypot(p.X-q.X,p.Y-q.Y)<8.f;
        });
        if (ref.equilibrium) {
            g.FillEllipse(&red,p.X-3,p.Y-3,6.f,6.f);
        } else {
            const float radius = coincident ? 6.f : 4.f;
            if (!coincident && f.kind != "observable_heatmap") g.FillEllipse(&background,p.X-radius,p.Y-radius,2*radius,2*radius);
            g.DrawEllipse(&pen,p.X-radius,p.Y-radius,2*radius,2*radius);
        }
        if(f.kind=="phase"&&(ref.equilibrium||!coincident)) {
            const auto label=wide(coincident?"IC / EC":ref.label);
            const float markerRadius=coincident?6.f:ref.equilibrium?3.f:4.f;
            const auto bounds=plotReferenceLabelBounds(p,markerRadius,
                measureMathText(g,label,fonts.caption),area);
            drawMathText(g,label,bounds,fonts.caption,theme::plotMuted);
        }
    }
    g.Restore(state);
}
void drawTimeColorbar(Gdiplus::Graphics& g,const Figure& f,const Gdiplus::RectF& plot,const theme::Fonts& fonts) {
    const Trace* colored=nullptr;
    for(const auto& trace:f.traces) if(trace.encoding=="time"&&!trace.colorValues.empty()) {colored=&trace;break;}
    if(!colored) return;
    const auto values=trajectoryTimeRange(f);
    const float scale=scientificTextScale(fonts.body);
    const float width=std::min(140.f*scale,plot.Width*.38f),x=plot.GetRight()-width-50*scale,y=plot.Y-19*scale;
    for(int i=0;i<64;++i) {
        Gdiplus::SolidBrush brush(trajectoryTimeColor(i/63.));
        g.FillRectangle(&brush,x+width*i/64,y,width/64+.5f,5.f*scale);
    }
    wchar_t low[32]{},high[32]{};
    swprintf_s(low,L"%.3g",values.first);swprintf_s(high,L"%.3g",values.second);
    drawMathText(g,low,{x-50*scale,y-9*scale,46*scale,20*scale},fonts.caption,theme::plotMuted,Gdiplus::StringAlignmentFar);
    drawMathText(g,high,{x+width+4*scale,y-9*scale,48*scale,20*scale},fonts.caption,theme::plotMuted);
    drawMathText(g,wide(f.colorLabel),{x-88*scale,y-9*scale,34*scale,20*scale},fonts.caption,theme::plotMuted);
}
}
