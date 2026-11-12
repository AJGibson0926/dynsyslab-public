#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "PlotFigureDrawing.hpp"
#include <algorithm>
#include <array>
#include <cmath>
namespace dsl {
namespace {
constexpr std::size_t budget = 64u * 1024u * 1024u;
bool sameLimits(const std::optional<AxisLimits>& a, const std::optional<AxisLimits>& b) {
    return a.has_value() == b.has_value() && (!a ||
        (a->xMin == b->xMin && a->xMax == b->xMax && a->yMin == b->yMin && a->yMax == b->yMax));
}
void quality(Gdiplus::Graphics& g) {
    g.SetSmoothingMode(Gdiplus::SmoothingModeAntiAlias);
    g.SetTextRenderingHint(Gdiplus::TextRenderingHintAntiAliasGridFit);
    g.SetPixelOffsetMode(Gdiplus::PixelOffsetModeHalf);
}
}
struct PlotCanvas::RasterCache {
    struct Panel {
        float width{}, height{};
        std::array<HFONT,4> fonts{};
        std::optional<AxisLimits> limits;
        std::optional<Plot3DCameraState> camera3D;
        Plot3DControl hovered3DControl = Plot3DControl::none;
        PlotTransform transform;
        std::shared_ptr<Gdiplus::Bitmap> image;
        std::size_t bytes{}, used{};
    };
    std::map<std::size_t,Panel> panels;
    std::map<std::size_t,Gdiplus::RectF> pending;
    std::size_t bytes{}, clock{};
};
PlotRenderStatistics PlotCanvas::renderStatistics() const {
    return {panelRenders_,cacheHits_,rasterCache_ ? rasterCache_->bytes : 0};
}
PlotTransform PlotCanvas::drawCachedPanel(Gdiplus::Graphics& g, std::size_t index,
    const Figure& figure, const Gdiplus::RectF& card, const std::optional<AxisLimits>& limits,
    const Plot3DCameraState* camera3D, Plot3DControl hovered3DControl) {
    if (!rasterCache_) rasterCache_ = std::make_shared<RasterCache>();
    auto& cache = *rasterCache_;
    cache.pending.erase(index);
    const std::array<HFONT,4> fonts{fonts_->body,fonts_->caption,fonts_->heading,fonts_->title};
    auto found = cache.panels.find(index);
    if (found != cache.panels.end() && (found->second.width != card.Width || found->second.height != card.Height
        || found->second.fonts != fonts || !sameLimits(found->second.limits,limits)
        || found->second.camera3D.has_value() != (camera3D != nullptr)
        || (camera3D && found->second.camera3D != *camera3D)
        || found->second.hovered3DControl != hovered3DControl)) {
        cache.bytes -= found->second.bytes; cache.panels.erase(found); found = cache.panels.end();
    }
    if (found == cache.panels.end()) {
        const int width = std::max(1,static_cast<int>(std::ceil(card.Width)));
        const int height = std::max(1,static_cast<int>(std::ceil(card.Height)));
        const auto bytes = static_cast<std::size_t>(width)*height*4;
        if (bytes > budget || width > 16384 || height > 16384) {
            const auto transform = renderTransform(figure,card,limits,*fonts_);
            const Plot3DCamera camera(camera3D ? *camera3D : Plot3DCameraState{});
            figureCard(g,figure,card,*fonts_,transform,
                figure.kind == "three_d" ? &camera : nullptr, hovered3DControl);
            ++panelRenders_; return transform;
        }
        while (cache.bytes + bytes > budget && !cache.panels.empty()) {
            const auto oldest = std::min_element(cache.panels.begin(),cache.panels.end(),
                [](const auto& a,const auto& b) { return a.second.used < b.second.used; });
            cache.bytes -= oldest->second.bytes; cache.panels.erase(oldest);
        }
        RasterCache::Panel panel;
        panel.width = card.Width; panel.height = card.Height; panel.fonts = fonts; panel.limits = limits;
        if (camera3D) panel.camera3D = *camera3D;
        panel.hovered3DControl = hovered3DControl;
        panel.bytes = bytes;
        const Gdiplus::RectF local(0,0,card.Width,card.Height);
        panel.transform = renderTransform(figure,local,limits,*fonts_);
        const int scale = bytes <= 8000000u*4u ? 2 : 1;
        Gdiplus::Bitmap detailed(width*scale,height*scale,PixelFormat32bppPARGB);
        {
            Gdiplus::Graphics detail(&detailed); quality(detail); detail.ScaleTransform(static_cast<float>(scale),static_cast<float>(scale));
            detail.Clear(theme::color(theme::background));
            const Plot3DCamera camera(camera3D ? *camera3D : Plot3DCameraState{});
            figureCard(detail,figure,local,*fonts_,panel.transform,
                figure.kind == "three_d" ? &camera : nullptr, hovered3DControl);
        }
        panel.image = std::make_shared<Gdiplus::Bitmap>(width,height,PixelFormat32bppPARGB);
        {
            Gdiplus::Graphics reduced(panel.image.get()); quality(reduced);
            reduced.SetInterpolationMode(Gdiplus::InterpolationModeHighQualityBicubic);
            reduced.DrawImage(&detailed,Gdiplus::Rect(0,0,width,height),0,0,width*scale,height*scale,Gdiplus::UnitPixel);
        }
        cache.bytes += bytes; ++panelRenders_;
        found = cache.panels.emplace(index,std::move(panel)).first;
    } else ++cacheHits_;
    auto& panel = found->second; panel.used = ++cache.clock;
    const auto saved = g.Save();
    g.SetInterpolationMode(Gdiplus::InterpolationModeNearestNeighbor);
    g.DrawImage(panel.image.get(),Gdiplus::RectF(card.X,card.Y,
        static_cast<float>(panel.image->GetWidth()),static_cast<float>(panel.image->GetHeight())),
        0,0,static_cast<float>(panel.image->GetWidth()),static_cast<float>(panel.image->GetHeight()),Gdiplus::UnitPixel);
    g.Restore(saved);
    auto transform = panel.transform; transform.left += card.X; transform.top += card.Y;
    return transform;
}
void PlotCanvas::beginCacheFrame() {
    if (!rasterCache_) rasterCache_ = std::make_shared<RasterCache>();
    rasterCache_->pending.clear();
}
void PlotCanvas::queueCachedPanel(std::size_t index, const Gdiplus::RectF& card) {
    auto& cache = *rasterCache_;
    const auto found = cache.panels.find(index);
    if (found == cache.panels.end() || found->second.width != card.Width || found->second.height != card.Height)
        cache.pending[index] = card;
}
void PlotCanvas::scheduleCacheWarmup() {
    if (window_ && rasterCache_ && !rasterCache_->pending.empty() && rasterCache_->bytes < budget)
        SetTimer(window_,cacheTimerId_,30,nullptr);
}
void PlotCanvas::warmPlotCache() {
    KillTimer(window_,cacheTimerId_);
    if (!payload_ || !rasterCache_ || scrollTick_ || dragMode_ != DragMode::none) return;
    auto& pending = rasterCache_->pending;
    if (pending.empty()) return;
    const auto closest = std::min_element(pending.begin(),pending.end(),[&](const auto& a,const auto& b) {
        const auto distance = [&](const auto& card) {
            return std::max(0.f,std::max(-card.GetBottom(),card.Y-scrollViewport_.bottom));
        };
        return distance(a.second) < distance(b.second);
    });
    const auto index = closest->first; const auto card = closest->second;
    const auto bytes = static_cast<std::size_t>(std::ceil(card.Width))*static_cast<std::size_t>(std::ceil(card.Height))*4;
    if (index >= payload_->figures.size() || rasterCache_->bytes+bytes > budget) { pending.clear(); return; }
    const auto found = axisLimits_.find(index);
    const auto limits = found == axisLimits_.end() ? std::nullopt : std::optional<AxisLimits>(found->second);
    Gdiplus::Bitmap scratch(1,1,PixelFormat32bppPARGB);Gdiplus::Graphics graphics(&scratch);
    const auto camera = cameras3D_.find(index);
    const auto hovered = hovered3DFigure_ && *hovered3DFigure_ == index
        ? hovered3DControl_ : Plot3DControl::none;
    drawCachedPanel(graphics,index,payload_->figures[index],{0,0,card.Width,card.Height},limits,
        camera == cameras3D_.end() ? nullptr : &camera->second, hovered);
    scheduleCacheWarmup();
}

}
