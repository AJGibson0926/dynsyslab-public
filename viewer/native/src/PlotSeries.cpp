#include "PlotSeries.hpp"
#include <cstdio>
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "MathTextDrawing.hpp"
#include "DomainColor.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace dsl {
namespace {
Gdiplus::PointF project(double x,double y,const PlotTransform& t) {
    return {static_cast<float>(t.left+(x-t.xMin)/(t.xMax-t.xMin)*t.width),
        static_cast<float>(plotPixelY(t,y))};
}
struct PixelPoint { double x, y; };
PixelPoint projectExact(double x,double y,const PlotTransform& t) {
    return {t.left+(x-t.xMin)/(t.xMax-t.xMin)*t.width, plotPixelY(t,y)};
}
// Liang-Barsky clip of segment a-b to a box. Points far outside the plot (a run that
// diverged to 1e7) overflow GDI+ coordinates and make it drop the whole polyline, so
// every segment is clipped in double precision before it reaches GDI+.
struct PixelBox { double left, top, right, bottom; };
bool clipSegment(PixelPoint& a,PixelPoint& b,const PixelBox& box) {
    double t0=0., t1=1.; const double dx=b.x-a.x, dy=b.y-a.y;
    const double p[4]={-dx,dx,-dy,dy}, q[4]={a.x-box.left,box.right-a.x,a.y-box.top,box.bottom-a.y};
    for(int k=0;k<4;++k) {
        if(p[k]==0.) { if(q[k]<0.) return false; continue; }
        const double r=q[k]/p[k];
        if(p[k]<0.) { if(r>t1) return false; t0=std::max(t0,r); }
        else { if(r<t0) return false; t1=std::min(t1,r); }
    }
    const PixelPoint start{a.x+t0*dx,a.y+t0*dy}, end{a.x+t1*dx,a.y+t1*dy};
    a=start; b=end; return true;
}
bool inside(const PixelPoint& p,const PixelBox& box) {
    return p.x>=box.left&&p.x<=box.right&&p.y>=box.top&&p.y<=box.bottom;
}
Gdiplus::PointF toPoint(const PixelPoint& p) { return {static_cast<float>(p.x),static_cast<float>(p.y)}; }
std::pair<double,double> timeRange(const std::vector<double>& values) {
    double low=std::numeric_limits<double>::infinity(),high=-low;
    for(const auto value:values) if(std::isfinite(value)) { low=std::min(low,value);high=std::max(high,value); }
    return std::isfinite(low)?std::make_pair(low,high):std::make_pair(0.,1.);
}
}
float plotSeriesWidthPoints(const Trace& trace) {
    if (std::isfinite(trace.lineWidthPoints) && trace.lineWidthPoints>0.)
        return static_cast<float>(std::min(32.,trace.lineWidthPoints));
    if (trace.role == "phase_vector_field") return 1.3f;
    if (trace.role.rfind("phase_",0)==0) return 2.1f;
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
    } else if (trace.lineStyle == "sparse-dash") {
        // Long dashes with wide gaps (in pen widths): an overlaid reference leaves the curve beneath visible.
        static const Gdiplus::REAL sparse[] = {5.0f, 4.0f};
        pen.SetDashPattern(sparse, 2);
    } else if (trace.lineStyle == "dot") {
        static const Gdiplus::REAL dots[] = {0.6f, 2.4f};
        pen.SetDashPattern(dots, 2);
    } else if (trace.lineStyle == "dash-dot") {
        static const Gdiplus::REAL dashDot[] = {5.0f, 2.5f, 0.8f, 2.5f};
        pen.SetDashPattern(dashDot, 4);
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
    const double maxMagnitude=trace.magnitudeScale>0?trace.magnitudeScale:trace.magnitudes.empty()?1.:
        *std::max_element(trace.magnitudes.begin(),trace.magnitudes.end());
    std::vector<Gdiplus::PointF> points;
    auto flush=[&](){if(points.size()>1) g.DrawLines(&pen,points.data(),static_cast<INT>(points.size()));points.clear();};
    // Guard box: the plot area plus one plot size on every side (the visible clip is the plot).
    const PixelBox box{t.left-t.width,t.top-t.height,t.left+2.*t.width,t.top+2.*t.height};
    std::optional<PixelPoint> previousExact;
    for(std::size_t i=0;i<count;++i) {
        double x=trace.x[i],y=trace.y[i];
        if(i<trace.z.size()) {x-=.32*trace.z[i];y+=.20*trace.z[i];}
        if(!std::isfinite(x)||!std::isfinite(y)) {flush();previousExact.reset();continue;}
        const auto exact=projectExact(x,y,t);
        const bool visible=inside(exact,box);
        const auto point=toPoint(exact);
        if(trace.colorValues.empty()) {
            // A "none" line style is a scatter: points only, never joined in storage order.
            if(trace.lineStyle!="none") {
                if(!previousExact) { if(visible) points.push_back(point); }
                else {
                    PixelPoint a=*previousExact,b=exact;
                    if(!clipSegment(a,b,box)) flush();
                    else {
                        if(!inside(*previousExact,box)) { flush(); points.push_back(toPoint(a)); }
                        else if(points.empty()) points.push_back(toPoint(a));
                        points.push_back(toPoint(b));
                        if(!visible) flush();
                    }
                }
            }
        } else if(previousExact) {
            PixelPoint a=*previousExact,b=exact;
            const bool drawSegment=clipSegment(a,b,box);
            const auto colorIndex=trace.encoding=="time"?i-1:i;
            const auto value=colorIndex<trace.colorValues.size()?trace.colorValues[colorIndex]:0.;
            if(trace.encoding=="time") pen.SetColor(applyPlotSeriesOpacity(trajectoryTimeColor(
                times.second>times.first?(value-times.first)/(times.second-times.first):0.),trace));
            else if(trace.encoding=="domain-color") pen.SetColor(applyPlotSeriesOpacity(
                domainColor(value,i<trace.magnitudes.size()?trace.magnitudes[i]:1.),trace));
            else pen.SetColor(applyPlotSeriesOpacity(theme::phaseColor(value,
                i<trace.magnitudes.size()&&maxMagnitude>0?trace.magnitudes[i]/maxMagnitude:1.),trace));
            if(drawSegment) g.DrawLine(&pen,toPoint(a),toPoint(b));
        }
        if((trace.marker||count==1)&&visible) {
            Gdiplus::SolidBrush brush(plotSeriesColor(trace,traceIndex));
            if (trace.selectionId >= 0) {
                const float radius = trace.selected ? 5.f : 3.5f;
                if (trace.selected) g.FillEllipse(&brush,point.X-radius,point.Y-radius,radius*2,radius*2);
                Gdiplus::Pen ring(trace.selected ? theme::color(theme::plotText) : theme::color(theme::plotMuted),1.2f);
                g.DrawEllipse(&ring,point.X-radius,point.Y-radius,radius*2,radius*2);
            } else if (trace.markerShape == "open-circle") {
                Gdiplus::Pen ring(plotSeriesColor(trace,traceIndex),1.3f);
                g.DrawEllipse(&ring,point.X-4.5f,point.Y-4.5f,9.f,9.f);
            } else if (trace.markerShape == "large-circle") {
                g.FillEllipse(&brush,point.X-9.f,point.Y-9.f,18.f,18.f);
                Gdiplus::Pen edge(theme::color(theme::plotText),1.6f);
                g.DrawEllipse(&edge,point.X-9.f,point.Y-9.f,18.f,18.f);
            } else if (trace.markerShape == "large-open-circle") {
                Gdiplus::Pen ring(plotSeriesColor(trace,traceIndex),4.f);
                g.DrawEllipse(&ring,point.X-9.f,point.Y-9.f,18.f,18.f);
            } else g.FillEllipse(&brush,point.X-2.2f,point.Y-2.2f,4.4f,4.4f);
        }
        previousExact=exact;
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
    std::vector<Gdiplus::RectF> placed;  // labels already drawn: later ones move to a clear side
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
            // An initial condition is labelled with its coordinates in the current display units.
            char coordinates[64]{};
            std::snprintf(coordinates,sizeof coordinates,"(%.3g, %.3g)",ref.x,ref.y);
            const auto label=wide(coincident?"IC / EC":(!ref.equilibrium&&ref.label=="IC")?std::string(coordinates):ref.label);
            const float markerRadius=coincident?6.f:ref.equilibrium?3.f:4.f;
            auto bounds=plotReferenceLabelBounds(p,markerRadius,measureMathText(g,label,fonts.caption),area);
            const auto clear=[&](const Gdiplus::RectF& r) {
                return std::none_of(placed.begin(),placed.end(),[&](const auto& q) { return r.IntersectsWith(q); });
            };
            if(!clear(bounds)) {  // one, then two label heights above or below, keeping the same side
                const float h=bounds.Height;
                for(const float rows:{-1.f,1.f,-2.f,2.f}) {
                    auto candidate=bounds; candidate.Y=std::clamp(bounds.Y+rows*h,area.Y,std::max(area.Y,area.GetBottom()-h));
                    if(clear(candidate)) { bounds=candidate; break; }
                }
            }
            placed.push_back(bounds);
            drawMathText(g,label,bounds,fonts.caption,theme::plotMuted);
        }
    }
    g.Restore(state);
}
}
