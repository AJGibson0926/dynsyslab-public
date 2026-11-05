#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PanelWorkspace.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"

#include <windowsx.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace dsl {
namespace {

constexpr unsigned modifierControl = 1u;
constexpr unsigned modifierAlt = 2u;
constexpr unsigned modifierShift = 4u;

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

template <typename Value>
void moveIndexed(std::map<std::size_t, Value>& values,
                 std::size_t from, std::size_t to) {
    std::map<std::size_t, Value> adjusted;
    for (auto& entry : values) {
        std::size_t index = entry.first;
        if (index == from) index = to;
        else if (from < to && index > from && index <= to) --index;
        else if (to < from && index >= to && index < from) ++index;
        adjusted[index] = std::move(entry.second);
    }
    values = std::move(adjusted);
}

}  // namespace

void PlotCanvas::pointerDown(POINT point, unsigned modifiers, int button) {
    stopScroll();
    SetFocus(window_);
    if (verticallyScrollable() && scrollBar_.visible()
        && scrollBar_.contains(window_,scrollViewport_,point.x,point.y)) {
        scrollBar_.buttonDown(window_,scrollViewport_,point.x,point.y);
        hoverProbe_.reset(); InvalidateRect(window_,nullptr,FALSE); return;
    }
    if (selectKromAt(point,(modifiers & modifierControl) != 0)) return;
    for (const auto& hit : hitRegions_) {
        if (PtInRect(&hit.close, point)) { removePanel(verticallyScrollable() ? hit.figureIndex : hit.visualIndex); return; }
        if (PtInRect(&hit.resize, point)) {
            ensureComposition();
            if (!payload_ || hit.visualIndex >= payload_->figures.size()) return;
            dragMode_ = DragMode::resize; draggedPanel_ = verticallyScrollable() ? hit.figureIndex : hit.visualIndex;
            dragOrigin_ = point;
            dragStartColumnSpan_ = payload_->figures[draggedPanel_].gridSpan;
            dragStartRowSpan_ = payload_->figures[draggedPanel_].rowSpan;
            SetCapture(window_); return;
        }
        if (PtInRect(&hit.card, point) && !PtInRect(&hit.plot, point)) {
            ensureComposition(); dragMode_ = DragMode::reorder;
            draggedPanel_ = verticallyScrollable() ? hit.figureIndex : hit.visualIndex; dragOrigin_ = point;
            SetCapture(window_); return;
        }
    }
    for (const auto& hit : hitRegions_) {
        if (!PtInRect(&hit.plot, point) || !payload_
            || hit.figureIndex >= payload_->figures.size()) continue;
        auto& figure = payload_->figures[hit.figureIndex];
        if (figure.kind == "three_d") {
            active3DFigure_ = hit.figureIndex;
            const Plot3DViewport viewport{hit.transform.left, hit.transform.top,
                hit.transform.width, hit.transform.height};
            const auto control = hitPlot3DControl(viewport, point.x, point.y);
            if (control != Plot3DControl::none && button == 0) {
                Plot3DCamera camera(cameras3D_[hit.figureIndex]);
                camera.apply(control); cameras3D_[hit.figureIndex] = camera.state();
                rasterCache_.reset(); InvalidateRect(window_, nullptr, FALSE); return;
            }
            if (button == 0 && modifiers == 0) {
                const auto pinned = pinnedProbes_.find(hit.figureIndex);
                if (pinned != pinnedProbes_.end()
                    && std::hypot(pinned->second.pixelX - point.x,
                                  pinned->second.pixelY - point.y) <= 10.0) {
                    dragMode_ = DragMode::probe; draggedFigure_ = hit.figureIndex;
                    SetCapture(window_); return;
                }
                const Plot3DCamera camera(cameras3D_[hit.figureIndex]);
                const auto probe = nearestPlot3DProbe(figure, camera, viewport,
                    point.x, point.y, 9.0);
                if (probe) {
                    pinnedProbes_[hit.figureIndex] = curveProbe(*probe);
                    hoverProbe_ = CanvasProbe{hit.figureIndex, curveProbe(*probe)};
                    dragMode_ = DragMode::probe; draggedFigure_ = hit.figureIndex;
                    InvalidateRect(window_, nullptr, FALSE); SetCapture(window_); return;
                }
            }
            dragMode_ = button == 1 || (modifiers & modifierControl)
                ? DragMode::pan3d
                : button == 2 || (modifiers & modifierAlt)
                    ? DragMode::zoom3d : DragMode::orbit3d;
            draggedFigure_ = hit.figureIndex; dragOrigin_ = point;
            dragStartCamera3D_ = cameras3D_[hit.figureIndex];
            dragConstrain3D_ = (modifiers & modifierShift) != 0;
            hoverProbe_.reset(); SetCapture(window_); return;
        }
        const auto pinned = pinnedProbes_.find(hit.figureIndex);
        const std::optional<CurveProbe> pinnedValue = pinned == pinnedProbes_.end()
            ? std::nullopt : std::optional<CurveProbe>(pinned->second);
        const auto decision = inspectPlotPointer(payload_->figures[hit.figureIndex],
            hit.transform, pinnedValue, point.x, point.y, 10.0);
        if (decision.kind == PlotPointerKind::pinned_probe) {
            dragMode_ = DragMode::probe; draggedFigure_ = hit.figureIndex;
            SetCapture(window_); return;
        }
        if (decision.kind == PlotPointerKind::curve && decision.probe) {
            pinnedProbes_[hit.figureIndex] = *decision.probe;
            hoverProbe_ = CanvasProbe{hit.figureIndex, *decision.probe};
            dragMode_ = DragMode::probe; draggedFigure_ = hit.figureIndex;
            InvalidateRect(window_, nullptr, FALSE); SetCapture(window_); return;
        }
        if (backgroundDragMode_ == PlotBackgroundDragMode::rearrange_panels) {
            ensureComposition(); dragMode_ = DragMode::reorder;
            draggedPanel_ = verticallyScrollable() ? hit.figureIndex : hit.visualIndex; dragOrigin_ = point;
            hoverProbe_.reset(); SetCapture(window_); return;
        }
        dragMode_ = DragMode::pan; draggedFigure_ = hit.figureIndex;
        dragOrigin_ = point; dragStartTransform_ = hit.transform;
        hoverProbe_.reset(); SetCapture(window_); return;
    }
}

void PlotCanvas::pointerMove(POINT point) {
    if (verticallyScrollable()) {
        if (scrollBar_.mouseMove(window_,scrollViewport_,point.x,point.y)) {
            hoverProbe_.reset(); InvalidateRect(window_,nullptr,FALSE); SetCursor(LoadCursorW(nullptr,IDC_ARROW)); return;
        }
        if (scrollBar_.visible() && scrollBar_.contains(window_,scrollViewport_,point.x,point.y)) {
            SetCursor(LoadCursorW(nullptr,IDC_ARROW)); return;
        }
    }
    if (dragMode_ == DragMode::reorder && payload_) {
        auto target = hitPanel(point);
        if (verticallyScrollable() && target) {
            const auto hit = std::find_if(hitRegions_.begin(),hitRegions_.end(),
                [&](const auto& region) { return region.visualIndex == *target; });
            if (hit == hitRegions_.end() || payload_->figures[hit->figureIndex].kind == "krom_spectrum") target.reset();
            else target = hit->figureIndex;
        }
        if (target && *target != draggedPanel_ &&
            movePanel(payload_->figures, draggedPanel_, *target)) {
            rasterCache_.reset();
            moveIndexed(pinnedProbes_, draggedPanel_, *target);
            moveIndexed(axisLimits_, draggedPanel_, *target);
            moveIndexed(cameras3D_, draggedPanel_, *target);
            draggedPanel_ = *target;
            hoverProbe_.reset();
            InvalidateRect(window_, nullptr, FALSE);
        }
        SetCursor(LoadCursorW(nullptr, IDC_SIZEALL)); return;
    }
    if (dragMode_ == DragMode::resize && payload_ &&
        draggedPanel_ < payload_->figures.size()) {
        const int columnDelta = static_cast<int>(std::lround(
            static_cast<double>(point.x - dragOrigin_.x) /
            std::max(1.0f, layoutCellWidth_)));
        const int rowDelta = static_cast<int>(std::lround(
            static_cast<double>(point.y - dragOrigin_.y) /
            std::max(1.0f, layoutCellHeight_)));
        auto& figure = payload_->figures[draggedPanel_];
        const int columns = std::clamp(dragStartColumnSpan_ + columnDelta,
            1, std::max(layoutColumns_, 1));
        const int rows = std::clamp(dragStartRowSpan_ + rowDelta,
            1, std::max(layoutRows_, 1));
        if (figure.gridSpan != columns || figure.rowSpan != rows) {
            figure.gridSpan = columns; figure.rowSpan = rows;
            InvalidateRect(window_, nullptr, FALSE);
        }
        SetCursor(LoadCursorW(nullptr, IDC_SIZENWSE)); return;
    }
    if (dragMode_ == DragMode::probe && payload_) {
        const auto pinned = pinnedProbes_.find(draggedFigure_);
        for (const auto& hit : hitRegions_) {
            if (hit.figureIndex != draggedFigure_ || pinned == pinnedProbes_.end()
                || draggedFigure_ >= payload_->figures.size()) continue;
            std::optional<CurveProbe> probe;
            if (payload_->figures[draggedFigure_].kind == "three_d") {
                const Plot3DCamera camera(cameras3D_[draggedFigure_]);
                const auto value = nearestPlot3DProbe(payload_->figures[draggedFigure_],
                    camera, {hit.transform.left, hit.transform.top,
                        hit.transform.width, hit.transform.height}, point.x, point.y,
                    std::numeric_limits<double>::infinity(), pinned->second.traceIndex);
                if (value) probe = curveProbe(*value);
            } else probe = nearestCurveProbe(payload_->figures[draggedFigure_],
                hit.transform, point.x, point.y,
                std::numeric_limits<double>::infinity(), pinned->second.traceIndex);
            if (probe) {
                pinned->second = *probe;
                hoverProbe_ = CanvasProbe{draggedFigure_, *probe};
                InvalidateRect(window_, nullptr, FALSE);
            }
            break;
        }
        SetCursor(LoadCursorW(nullptr, IDC_ARROW)); return;
    }
    if (dragMode_ == DragMode::pan) {
        const auto moved = panPlotTransform(dragStartTransform_,
            point.x - dragOrigin_.x, point.y - dragOrigin_.y);
        axisLimits_[draggedFigure_] = {moved.xMin, moved.xMax,
                                      moved.yMin, moved.yMax};
        hoverProbe_.reset(); InvalidateRect(window_, nullptr, FALSE);
        SetCursor(LoadCursorW(nullptr, IDC_SIZEALL)); return;
    }
    if ((dragMode_ == DragMode::orbit3d || dragMode_ == DragMode::pan3d
         || dragMode_ == DragMode::zoom3d) && payload_
        && draggedFigure_ < payload_->figures.size()) {
        const auto hit = std::find_if(hitRegions_.begin(), hitRegions_.end(),
            [&](const auto& region) { return region.figureIndex == draggedFigure_; });
        if (hit != hitRegions_.end()) {
            Plot3DCamera camera(dragStartCamera3D_);
            const double dx = point.x - dragOrigin_.x;
            const double dy = point.y - dragOrigin_.y;
            if (dragMode_ == DragMode::orbit3d) camera.orbit(dx, dy, dragConstrain3D_);
            else if (dragMode_ == DragMode::pan3d) camera.pan(dx, dy,
                {hit->transform.left, hit->transform.top,
                 hit->transform.width, hit->transform.height});
            else camera.zoomDrag(dy);
            cameras3D_[draggedFigure_] = camera.state();
            rasterCache_.reset(); hoverProbe_.reset();
            InvalidateRect(window_, nullptr, FALSE);
        }
        SetCursor(LoadCursorW(nullptr, dragMode_ == DragMode::zoom3d
            ? IDC_SIZENS : IDC_SIZEALL)); return;
    }
    const auto hovered = hitPanel(point);
    if (hovered != hoveredPanel_) {
        hoveredPanel_ = hovered; InvalidateRect(window_, nullptr, FALSE);
    }
    TRACKMOUSEEVENT tracking{sizeof(tracking), TME_LEAVE, window_, 0};
    TrackMouseEvent(&tracking);
    for (const auto& hit : hitRegions_) {
        if (isKrom() && payload_->figures[hit.figureIndex].kind == "krom_spectrum" && PtInRect(&hit.card,point)) {
            SetCursor(LoadCursorW(nullptr,IDC_ARROW)); return;
        }
        if (PtInRect(&hit.close, point)) {
            if (hoverProbe_) { hoverProbe_.reset(); InvalidateRect(window_, nullptr, FALSE); }
            SetCursor(LoadCursorW(nullptr, IDC_ARROW)); return;
        }
        if (PtInRect(&hit.resize, point)) {
            SetCursor(LoadCursorW(nullptr, IDC_SIZENWSE)); return;
        }
        if (PtInRect(&hit.card, point) && !PtInRect(&hit.plot, point)) {
            SetCursor(LoadCursorW(nullptr, IDC_SIZEALL)); return;
        }
        if (PtInRect(&hit.plot, point) && payload_
            && hit.figureIndex < payload_->figures.size()) {
            if (payload_->figures[hit.figureIndex].kind == "three_d") {
                const Plot3DViewport viewport{hit.transform.left, hit.transform.top,
                    hit.transform.width, hit.transform.height};
                const auto control = hitPlot3DControl(viewport, point.x, point.y);
                const bool changed = !hovered3DFigure_
                    || *hovered3DFigure_ != hit.figureIndex
                    || hovered3DControl_ != control;
                hovered3DFigure_ = hit.figureIndex;
                hovered3DControl_ = control;
                active3DFigure_ = hit.figureIndex;
                if (changed) { rasterCache_.reset(); InvalidateRect(window_, nullptr, FALSE); }
                std::optional<CanvasProbe> candidate;
                if (control == Plot3DControl::none) {
                    const Plot3DCamera camera(cameras3D_[hit.figureIndex]);
                    const auto probe = nearestPlot3DProbe(payload_->figures[hit.figureIndex],
                        camera, viewport, point.x, point.y, 9.0);
                    if (probe) candidate = CanvasProbe{hit.figureIndex, curveProbe(*probe)};
                }
                const bool probeChanged = hoverProbe_.has_value() != candidate.has_value()
                    || (hoverProbe_ && candidate
                        && (hoverProbe_->figureIndex != candidate->figureIndex
                            || hoverProbe_->probe.segmentIndex != candidate->probe.segmentIndex
                            || std::abs(hoverProbe_->probe.pixelX - candidate->probe.pixelX) >= 0.25
                            || std::abs(hoverProbe_->probe.pixelY - candidate->probe.pixelY) >= 0.25));
                if (probeChanged) { hoverProbe_ = candidate; InvalidateRect(window_, nullptr, FALSE); }
                SetCursor(LoadCursorW(nullptr, control != Plot3DControl::none
                    ? IDC_HAND : candidate ? IDC_ARROW : IDC_SIZEALL));
                return;
            }
            const auto pinned = pinnedProbes_.find(hit.figureIndex);
            const std::optional<CurveProbe> pinnedValue = pinned == pinnedProbes_.end()
                ? std::nullopt : std::optional<CurveProbe>(pinned->second);
            const auto decision = inspectPlotPointer(payload_->figures[hit.figureIndex],
                hit.transform, pinnedValue, point.x, point.y);
            const std::optional<CanvasProbe> candidate = decision.kind == PlotPointerKind::curve
                ? std::optional<CanvasProbe>(CanvasProbe{hit.figureIndex, *decision.probe})
                : std::nullopt;
            const bool changed = hoverProbe_.has_value() != candidate.has_value()
                || (hoverProbe_ && candidate
                    && (hoverProbe_->figureIndex != candidate->figureIndex
                        || hoverProbe_->probe.traceIndex != candidate->probe.traceIndex
                        || std::abs(hoverProbe_->probe.pixelX - candidate->probe.pixelX) >= 0.25
                        || std::abs(hoverProbe_->probe.pixelY - candidate->probe.pixelY) >= 0.25));
            if (changed) {
                hoverProbe_ = candidate; InvalidateRect(window_, nullptr, FALSE);
            }
            SetCursor(LoadCursorW(nullptr,
                decision.kind == PlotPointerKind::background ? IDC_SIZEALL : IDC_ARROW));
            return;
        }
    }
    if (hovered3DFigure_ || hovered3DControl_ != Plot3DControl::none) {
        hovered3DFigure_.reset(); hovered3DControl_ = Plot3DControl::none;
        rasterCache_.reset(); InvalidateRect(window_, nullptr, FALSE);
    }
    if (hoverProbe_) { hoverProbe_.reset(); InvalidateRect(window_, nullptr, FALSE); }
    SetCursor(LoadCursorW(nullptr, IDC_ARROW));
}

void PlotCanvas::pointerUp() {
    scrollBar_.cancelDrag(window_);
    if (dragMode_ == DragMode::none) return;
    const bool compositionChanged = dragMode_ == DragMode::reorder
        || dragMode_ == DragMode::resize;
    dragMode_ = DragMode::none;
    if (GetCapture() == window_) ReleaseCapture();
    if (compositionChanged) notifyCompositionChanged();
}

LRESULT CALLBACK PlotCanvas::procedure(HWND window, UINT message,
                                        WPARAM wparam, LPARAM lparam) {
    auto* self = reinterpret_cast<PlotCanvas*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        self = static_cast<PlotCanvas*>(
            reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        self->window_ = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    }
    if (!self) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_SETCURSOR && reinterpret_cast<HWND>(wparam) == window
        && LOWORD(lparam) == HTCLIENT) {
        // pointerMove owns the cursor, including during captured drags. Prevent
        // DefWindowProc from restoring the class arrow before each mouse move.
        return TRUE;
    }
    if (message == WM_PAINT) {
        PAINTSTRUCT paint{}; HDC target = BeginPaint(window, &paint);
        self->paint(target, true); EndPaint(window, &paint); return 0;
    }
    if (message == WM_PRINTCLIENT) {
        self->paint(reinterpret_cast<HDC>(wparam), false); return 0;
    }
    if (message == WM_LBUTTONDBLCLK && self->zoomKromAt({GET_X_LPARAM(lparam),GET_Y_LPARAM(lparam)},0)) return 0;
    if (message == WM_SIZE) { self->stopScroll(); self->updateScrollLayout(); InvalidateRect(window,nullptr,FALSE); return 0; }
    if (message == WM_TIMER && wparam == cacheTimerId_) { self->warmPlotCache(); return 0; }
    if (message == WM_TIMER && wparam == scrollTimerId_) { self->advanceScroll(); return 0; }
    if (message == WM_DESTROY) { KillTimer(window,cacheTimerId_); self->stopScroll(); self->rasterCache_.reset(); return 0; }
    if (message == WM_MOUSEWHEEL) {
        POINT point{GET_X_LPARAM(lparam),GET_Y_LPARAM(lparam)}; ScreenToClient(window,&point);
        if (self->zoom3DAt(point, GET_WHEEL_DELTA_WPARAM(wparam) > 0 ? 1.12 : 1.0 / 1.12)) return 0;
        if (self->zoomKromAt(point,GET_WHEEL_DELTA_WPARAM(wparam) > 0 ? .7 : 1/.7)) return 0;
        if (self->verticallyScrollable()) {
            self->wheelScroll(GET_WHEEL_DELTA_WPARAM(wparam)); return 0;
        }
        if (!self->composed_)
            GET_WHEEL_DELTA_WPARAM(wparam) < 0 ? self->nextPage() : self->previousPage();
        return 0;
    }
    if (message == WM_LBUTTONDOWN) {
        unsigned modifiers = (wparam & MK_CONTROL) ? modifierControl : 0u;
        if (wparam & MK_SHIFT) modifiers |= modifierShift;
        if (GetKeyState(VK_MENU) < 0) modifiers |= modifierAlt;
        self->pointerDown({GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam)}, modifiers, 0); return 0;
    }
    if (message == WM_MBUTTONDOWN || message == WM_RBUTTONDOWN) {
        unsigned modifiers = (wparam & MK_CONTROL) ? modifierControl : 0u;
        if (wparam & MK_SHIFT) modifiers |= modifierShift;
        if (GetKeyState(VK_MENU) < 0) modifiers |= modifierAlt;
        self->pointerDown({GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam)}, modifiers,
            message == WM_MBUTTONDOWN ? 1 : 2); return 0;
    }
    if (message == WM_MOUSEMOVE) {
        self->pointerMove({GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam)}); return 0;
    }
    if (message == WM_LBUTTONUP || message == WM_MBUTTONUP || message == WM_RBUTTONUP) {
        self->pointerUp(); return 0;
    }
    if (message == WM_GETDLGCODE) return DLGC_WANTARROWS | DLGC_WANTCHARS;
    if (message == WM_KEYDOWN) {
        unsigned modifiers = GetKeyState(VK_CONTROL) < 0 ? modifierControl : 0u;
        if (GetKeyState(VK_MENU) < 0) modifiers |= modifierAlt;
        if (GetKeyState(VK_SHIFT) < 0) modifiers |= modifierShift;
        if (self->key3D(wparam, modifiers)) return 0;
    }
    if (message == WM_SETFOCUS) {
        if (!self->active3DFigure_ && self->payload_)
            for (const auto& hit : self->hitRegions_)
                if (hit.figureIndex < self->payload_->figures.size()
                    && self->payload_->figures[hit.figureIndex].kind == "three_d") {
                    self->active3DFigure_ = hit.figureIndex; break;
                }
        InvalidateRect(window, nullptr, FALSE); return 0;
    }
    if (message == WM_KILLFOCUS) {
        InvalidateRect(window, nullptr, FALSE); return 0;
    }
    if (message == WM_CAPTURECHANGED) { self->dragMode_ = DragMode::none; return 0; }
    if (message == WM_MOUSELEAVE) {
        if (self->dragMode_ == DragMode::none) self->hoverProbe_.reset();
        self->hoveredPanel_.reset(); self->hovered3DFigure_.reset();
        self->hovered3DControl_ = Plot3DControl::none; self->rasterCache_.reset();
        InvalidateRect(window, nullptr, FALSE); return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    return DefWindowProcW(window, message, wparam, lparam);
}

}  // namespace dsl
