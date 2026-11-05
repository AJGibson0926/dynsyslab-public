#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PanelWorkspace.hpp"
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "PlotPanelChrome.hpp"
#include "PlotFigureDrawing.hpp"
#include "PlotLegend.hpp"
#include "PlotSeries.hpp"
#include "PlotHeatmap.hpp"
#include "MathTextDrawing.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "PlotProbeOverlay.hpp"

#include <gdiplus.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace dsl {
namespace {
constexpr wchar_t kClassName[] = L"DynSysLabPlotCanvas";
constexpr int kFiguresPerPage = 4;

CurveProbe curveProbe(const Plot3DProbe& probe) {
    CurveProbe result;
    result.traceIndex = probe.traceIndex;
    result.segmentIndex = probe.segmentIndex;
    result.segmentFraction = probe.segmentFraction;
    result.x = probe.value.x; result.y = probe.value.y; result.z = probe.value.z;
    result.pixelX = probe.projected.x; result.pixelY = probe.projected.y;
    result.distancePixels = probe.distancePixels;
    return result;
}

std::vector<std::size_t> figureIndices(const Payload& payload, int page,
                                       bool composed = false) {
    const bool dashboard = std::any_of(payload.figures.begin(), payload.figures.end(),
        [](const Figure& figure) { return figure.dashboardPage >= 0; });
    std::vector<std::size_t> result;
    if (payload.layout == "krom" || payload.view == "observable_workspace") {
        for (std::size_t i = 0; i < payload.figures.size(); ++i)
            if (payload.figures[i].kind != "krom_spectrum") result.push_back(i);
        return result;
    }
    if (composed) {
        result.resize(payload.figures.size());
        for (std::size_t index = 0; index < result.size(); ++index) result[index] = index;
        return result;
    }
    if (dashboard || payload.layout == "krom") {
        for (std::size_t index = 0; index < payload.figures.size(); ++index)
            if (payload.figures[index].dashboardPage == page) result.push_back(index);
    } else {
        const std::size_t start = static_cast<std::size_t>(page * kFiguresPerPage);
        const std::size_t end = std::min(payload.figures.size(), start + kFiguresPerPage);
        for (std::size_t index = start; index < end; ++index) result.push_back(index);
    }
    return result;
}

}  // namespace

bool PlotCanvas::registerClass(HINSTANCE instance) {
    WNDCLASSEXW type{sizeof(type)};
    type.style = CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS;
    type.lpfnWndProc = procedure;
    type.hInstance = instance;
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.lpszClassName = kClassName;
    type.hbrBackground = CreateSolidBrush(theme::background);
    if (RegisterClassExW(&type)) return true;
    const auto error = GetLastError();
    DeleteObject(type.hbrBackground);
    return error == ERROR_CLASS_ALREADY_EXISTS;
}

bool PlotCanvas::create(HWND parent, const theme::Fonts& fonts) {
    fonts_ = &fonts;
    window_ = CreateWindowExW(0, kClassName,
        L"Interactive scientific plots. Three-dimensional plots support rotation, movement, zoom, and keyboard control.",
        WS_CHILD | WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS,
        0, 0, 0, 0, parent, nullptr, GetModuleHandleW(nullptr), this);
    return window_ != nullptr;
}

void PlotCanvas::setPayload(std::shared_ptr<Payload> payload) {
    stopScroll(); rasterCache_.reset();
    payload_ = std::move(payload);
    scrollBar_.scrollBy(-scrollBar_.position());
    axisLimits_.clear(); cameras3D_.clear(); page_ = 0; composed_ = false; hitRegions_.clear();
    hoverProbe_.reset(); pinnedProbes_.clear(); active3DFigure_.reset();
    hovered3DFigure_.reset(); hovered3DControl_ = Plot3DControl::none;
    if (!preferredCoordinates_.empty() && canDisplayCoordinates(preferredCoordinates_))
        setDisplayCoordinates(preferredCoordinates_);
    if (!preferredValueScale_.empty()) setValueScaling(preferredValueScale_);
    InvalidateRect(window_, nullptr, FALSE);
}

void PlotCanvas::setMessage(std::wstring title, std::wstring detail) {
    scrollBar_.scrollBy(-scrollBar_.position());
    stopScroll(); rasterCache_.reset();
    payload_.reset(); axisLimits_.clear(); cameras3D_.clear(); page_ = 0; composed_ = false;
    hitRegions_.clear(); hoverProbe_.reset(); pinnedProbes_.clear();
    active3DFigure_.reset(); hovered3DFigure_.reset();
    hovered3DControl_ = Plot3DControl::none;
    messageTitle_ = std::move(title); messageDetail_ = std::move(detail);
    InvalidateRect(window_, nullptr, FALSE);
}

int PlotCanvas::pageCount() const {
    if (!payload_ || payload_->figures.empty()) return 1;
    if (composed_ || verticallyScrollable()) return 1;
    int greatest = -1;
    for (const auto& figure : payload_->figures)
        greatest = std::max(greatest, figure.dashboardPage);
    if (greatest >= 0) return greatest + 1;
    return static_cast<int>((payload_->figures.size() + kFiguresPerPage - 1) / kFiguresPerPage);
}

std::wstring PlotCanvas::pageLabel(int page) const {
    if (payload_ && page >= 0 && static_cast<std::size_t>(page) < payload_->pageLabels.size())
        return wide(payload_->pageLabels[static_cast<std::size_t>(page)]);
    return L"Page " + std::to_wstring(page + 1);
}

std::vector<AxisPanelState> PlotCanvas::currentAxisPanels() const {
    std::vector<AxisPanelState> result;
    if (!payload_) return result;
    for (const auto index : figureIndices(*payload_, page_, composed_)) {
        const auto limits = automaticAxisLimits(payload_->figures[index]);
        AxisPanelState panel{wide(payload_->figures[index].title),
            limits, std::nullopt};
        const auto found = axisLimits_.find(index);
        if (found != axisLimits_.end()) panel.manual = found->second;
        result.push_back(std::move(panel));
    }
    return result;
}

bool PlotCanvas::setAxisLimits(std::size_t panel, std::optional<AxisLimits> limits) {
    if (!payload_) return false;
    const auto indices = figureIndices(*payload_, page_, composed_);
    if (panel >= indices.size()) return false;
    if (limits && (!std::isfinite(limits->xMin) || !std::isfinite(limits->xMax)
        || !std::isfinite(limits->yMin) || !std::isfinite(limits->yMax)
        || limits->xMin >= limits->xMax || limits->yMin >= limits->yMax
        || (ordinateScale(payload_->figures[indices[panel]]).kind == ValueScaleKind::logarithmic && limits->yMin <= 0))) return false;
    if (limits) axisLimits_[indices[panel]] = *limits;
    else axisLimits_.erase(indices[panel]);
    InvalidateRect(window_, nullptr, FALSE);
    return true;
}

void PlotCanvas::paint(HDC target, bool interactive) {
    RECT client{}; GetClientRect(window_, &client);
    const int width = std::max(1L, client.right), height = std::max(1L, client.bottom);
    const float renderScale = 1.0f;
    Gdiplus::Bitmap buffer(static_cast<int>(width * renderScale),
        static_cast<int>(height * renderScale), PixelFormat32bppPARGB);
    Gdiplus::Graphics graphics(&buffer);
    graphics.SetSmoothingMode(Gdiplus::SmoothingModeAntiAlias);
    graphics.SetTextRenderingHint(Gdiplus::TextRenderingHintAntiAliasGridFit);
    graphics.SetPixelOffsetMode(Gdiplus::PixelOffsetModeHalf);
    graphics.ScaleTransform(renderScale, renderScale);
    graphics.Clear(theme::color(theme::background));
    updateScrollLayout();
    beginCacheFrame();
    hitRegions_.clear();
    if (!payload_ || payload_->figures.empty()) {
        const float center = static_cast<float>(height) * 0.44f;
        drawText(graphics, messageTitle_, {40, center - 30, static_cast<float>(width - 80), 42},
                 fonts_->title, theme::text, Gdiplus::StringAlignmentCenter);
        drawText(graphics, messageDetail_, {60, center + 14, static_cast<float>(width - 120), 40},
                 fonts_->body, theme::muted, Gdiplus::StringAlignmentCenter);
    } else {
        std::vector<std::pair<std::size_t, const Figure*>> displayed;
        for (const auto index : figureIndices(*payload_, page_, composed_))
            displayed.emplace_back(index, &payload_->figures[index]);
        const int count = static_cast<int>(displayed.size());
        const int workspaceWidth = isKrom() ? static_cast<int>(width * .66) : width;
        const auto limits = [this](std::size_t index) -> std::optional<AxisLimits> {
            const auto found = axisLimits_.find(index);
            return found == axisLimits_.end() ? std::nullopt
                                              : std::optional<AxisLimits>(found->second);
        };
        auto renderPanel = [&](std::size_t visualIndex,
                               const std::pair<std::size_t, const Figure*>& entry,
                               const Gdiplus::RectF& card) {
            const Plot3DCameraState* cameraState = nullptr;
            if (entry.second->kind == "three_d")
                cameraState = &cameras3D_[entry.first];
            const auto hoveredControl = hovered3DFigure_ && *hovered3DFigure_ == entry.first
                ? hovered3DControl_ : Plot3DControl::none;
            const PlotTransform transform = drawCachedPanel(
                graphics, entry.first, *entry.second, card, limits(entry.first),
                cameraState, hoveredControl);
            const RECT cardRect = panelRectangle(card);
            PanelHitRegion hit;
            hit.visualIndex = visualIndex; hit.figureIndex = entry.first;
            hit.card = cardRect;
            hit.header = {cardRect.left, cardRect.top, cardRect.right,
                std::min(cardRect.bottom, cardRect.top + 30)};
            hit.close = {std::max(cardRect.left, cardRect.right - 34), cardRect.top,
                cardRect.right, std::min(cardRect.bottom, cardRect.top + 32)};
            hit.resize = {std::max(cardRect.left, cardRect.right - 25),
                std::max(cardRect.top, cardRect.bottom - 25), cardRect.right, cardRect.bottom};
            hit.plot = {static_cast<LONG>(std::floor(transform.left)),
                static_cast<LONG>(std::floor(transform.top)),
                static_cast<LONG>(std::ceil(transform.left + transform.width)),
                static_cast<LONG>(std::ceil(transform.top + transform.height))};
            hit.transform = transform;
            hitRegions_.push_back(hit);
            if (interactive) {
                if (entry.second->kind == "three_d" && active3DFigure_
                    && *active3DFigure_ == entry.first && GetFocus() == window_) {
                    Gdiplus::Pen focus(theme::color(theme::focus), 2.0f);
                    graphics.DrawRectangle(&focus,
                        static_cast<float>(hit.plot.left + 1),
                        static_cast<float>(hit.plot.top + 1),
                        static_cast<float>(std::max(1L, hit.plot.right - hit.plot.left - 2)),
                        static_cast<float>(std::max(1L, hit.plot.bottom - hit.plot.top - 2)));
                }
                if (hoverProbe_ && hoverProbe_->figureIndex == entry.first) {
                    std::optional<CurveProbe> resolved;
                    if (entry.second->kind == "three_d" && cameraState) {
                        const Plot3DCamera camera(*cameraState);
                        const auto probe = resolvePlot3DProbe(*entry.second, camera,
                            {transform.left, transform.top, transform.width, transform.height},
                            hoverProbe_->probe.traceIndex, hoverProbe_->probe.segmentIndex,
                            hoverProbe_->probe.segmentFraction);
                        if (probe) resolved = curveProbe(*probe);
                    } else resolved = resolveCurveProbe(*entry.second, transform,
                        hoverProbe_->probe.traceIndex, hoverProbe_->probe.segmentIndex,
                        hoverProbe_->probe.segmentFraction);
                    if (resolved) {
                        hoverProbe_->probe = *resolved;
                        drawCurveProbe(graphics, *entry.second, transform,
                                       *resolved, false, *fonts_);
                    }
                }
                const auto pinned = pinnedProbes_.find(entry.first);
                if (pinned != pinnedProbes_.end()) {
                    std::optional<CurveProbe> resolved;
                    if (entry.second->kind == "three_d" && cameraState) {
                        const Plot3DCamera camera(*cameraState);
                        const auto probe = resolvePlot3DProbe(*entry.second, camera,
                            {transform.left, transform.top, transform.width, transform.height},
                            pinned->second.traceIndex, pinned->second.segmentIndex,
                            pinned->second.segmentFraction);
                        if (probe) resolved = curveProbe(*probe);
                    } else resolved = resolveCurveProbe(*entry.second, transform,
                        pinned->second.traceIndex, pinned->second.segmentIndex,
                        pinned->second.segmentFraction);
                    if (resolved) {
                        pinned->second = *resolved;
                        drawCurveProbe(graphics, *entry.second, transform,
                                       *resolved, true, *fonts_);
                    }
                }
                if (entry.second->kind != "krom_spectrum") drawPanelChrome(graphics, card, hoveredPanel_ &&
                    *hoveredPanel_ == visualIndex, *fonts_);
            }
        };
        if (count == 0) {
            drawText(graphics, L"Select eigenvalues to show eigenfunctions.",
                {40, 40, static_cast<float>(workspaceWidth - 80), 40}, fonts_->body,
                theme::muted, Gdiplus::StringAlignmentCenter);
        } else if (verticallyScrollable()) {
            const int columns = layoutColumns_;
            const int availableWidth = workspaceWidth - OrnamentalScrollBar::preferredWidth(window_);
            const float gap = 8, margin = 8;
            const float cardWidth = (availableWidth-margin*2-gap*(columns-1))/columns;
            layoutCellWidth_ = cardWidth+gap; layoutCellHeight_ = scrollRowHeight_+gap;
            std::vector<PanelSpan> spans;
            for (const auto& entry : displayed) spans.push_back({entry.second->gridSpan,entry.second->rowSpan});
            const auto saved = graphics.Save();
            graphics.SetClip(Gdiplus::RectF(0,0,static_cast<float>(availableWidth),static_cast<float>(height)));
            for (const auto& placement : packPanelGrid(spans,columns)) {
                Gdiplus::RectF card(margin+placement.column*(cardWidth+gap),
                    margin+placement.row*(scrollRowHeight_+gap)-scrollBar_.position(),
                    cardWidth*placement.columnSpan+gap*(placement.columnSpan-1),
                    scrollRowHeight_*placement.rowSpan+gap*(placement.rowSpan-1));
                if (card.GetBottom() <= 0 || card.Y >= height) {
                    queueCachedPanel(displayed[placement.panel].first,card); continue;
                }
                renderPanel(placement.panel,displayed[placement.panel],card);
            }
            graphics.Restore(saved);
        } else if (count == 5 && displayed.front().second->gridSpan == 2) {
            const float gap = 8.0f, margin = 8.0f;
            const float cardWidth = (workspaceWidth - margin * 2 - gap) / 2.0f;
            const float availableHeight = height - margin * 2 - gap * 2;
            const float topHeight = availableHeight * 0.44f;
            const float componentHeight = (availableHeight - topHeight) / 2.0f;
            layoutColumns_ = 2; layoutRows_ = 3;
            layoutCellWidth_ = cardWidth; layoutCellHeight_ = componentHeight;
            renderPanel(0, displayed.front(), {margin, margin,
                static_cast<float>(workspaceWidth) - margin * 2, topHeight});
            for (int index = 1; index < count; ++index) {
                const int slot = index - 1, column = slot % 2, row = slot / 2 + 1;
                const auto& entry = displayed[static_cast<std::size_t>(index)];
                renderPanel(static_cast<std::size_t>(index), entry,
                    {margin + column * (cardWidth + gap),
                     margin + topHeight + gap + (row - 1) * (componentHeight + gap),
                     cardWidth, componentHeight});
            }
        } else {
            const int columns = adaptivePanelColumns(displayed.size(), workspaceWidth, height);
            std::vector<PanelSpan> spans;
            spans.reserve(displayed.size());
            for (const auto& entry : displayed)
                spans.push_back({entry.second->gridSpan, entry.second->rowSpan});
            const auto placements = packPanelGrid(spans, columns);
            const int rows = packedPanelRows(placements);
            const float gap = 8.0f, margin = 8.0f;
            const float cardWidth = (workspaceWidth - margin * 2 - gap * (columns - 1)) / columns;
            const float cardHeight = (height - margin * 2 - gap * (rows - 1)) / rows;
            layoutColumns_ = columns; layoutRows_ = rows;
            layoutCellWidth_ = cardWidth + gap; layoutCellHeight_ = cardHeight + gap;
            for (const auto& placement : placements) {
                const auto& entry = displayed[placement.panel];
                renderPanel(placement.panel, entry,
                    {margin + placement.column * (cardWidth + gap),
                     margin + placement.row * (cardHeight + gap),
                     cardWidth * placement.columnSpan + gap * (placement.columnSpan - 1),
                     cardHeight * placement.rowSpan + gap * (placement.rowSpan - 1)});
            }
        }
        if (isKrom()) for (std::size_t i = 0; i < payload_->figures.size(); ++i)
            if (payload_->figures[i].kind == "krom_spectrum") {
                renderPanel(displayed.size(), {i,&payload_->figures[i]},
                    {static_cast<float>(workspaceWidth+4),8.f,static_cast<float>(width-workspaceWidth-12),static_cast<float>(height-16)});
                drawText(graphics,L"Click a mode • Ctrl+click adds/removes",
                    {static_cast<float>(workspaceWidth+14),45.f,static_cast<float>(width-workspaceWidth-32),30.f},
                    fonts_->caption,theme::plotMuted,Gdiplus::StringAlignmentCenter);
                drawText(graphics,L"Wheel zooms • Double-click resets",
                    {static_cast<float>(workspaceWidth+14),73.f,static_cast<float>(width-workspaceWidth-32),30.f},
                    fonts_->caption,theme::plotMuted,Gdiplus::StringAlignmentCenter);
            }
        if (pageCount() > 1) {
            const std::wstring page = pageLabel(page_) + L"  •  " +
                std::to_wstring(page_ + 1) + L" / " + std::to_wstring(pageCount());
            const float captionWidth = std::max(230.f,measureMathText(graphics,page,fonts_->caption).Width);
            drawMathText(graphics, page, {workspaceWidth-captionWidth-20,
                static_cast<float>(height - 30), captionWidth, 20}, fonts_->caption, theme::accentStrong,
                Gdiplus::StringAlignmentFar);
        }
    }
    Gdiplus::Graphics screen(target);
    screen.SetInterpolationMode(Gdiplus::InterpolationModeHighQualityBicubic);
    screen.SetPixelOffsetMode(Gdiplus::PixelOffsetModeHalf);
    screen.DrawImage(&buffer, Gdiplus::Rect(0, 0, width, height), 0, 0,
        buffer.GetWidth(), buffer.GetHeight(), Gdiplus::UnitPixel);
    if (verticallyScrollable() && scrollBar_.visible()) scrollBar_.paint(window_,target,scrollViewport_);
    scheduleCacheWarmup();
}

}  // namespace dsl
