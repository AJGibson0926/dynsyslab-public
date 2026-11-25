#pragma once

#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Plot3D.hpp"
#include "dynsyslab_viewer/Theme.hpp"
#include "dynsyslab_viewer/OrnamentalScrollBar.hpp"

#include <windows.h>

#include <memory>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace dsl {

constexpr UINT kPlotPageChangedMessage = WM_APP + 45;
constexpr UINT kPlotCompositionChangedMessage = WM_APP + 46;
constexpr UINT kPlotKromChangedMessage = WM_APP + 47;

struct AxisPanelState {
    std::wstring title;
    AxisLimits automatic;
    std::optional<AxisLimits> manual;
};

enum class PlotBackgroundDragMode { rearrange_panels, pan_curves };

struct PlotCanvasState {
    std::shared_ptr<Payload> payload;
    std::map<std::size_t, AxisLimits> axisLimits;
    std::map<std::size_t, CurveProbe> pinnedProbes;
    std::map<std::size_t, Plot3DCameraState> cameras3D;
    int page = 0;
    bool composed = false;
    int scrollY = 0;
};

struct PlotRenderStatistics {
    std::size_t panelRenders = 0, cacheHits = 0, cacheBytes = 0;
};

class PlotCanvas {
public:
    static bool registerClass(HINSTANCE instance);
    bool create(HWND parent, const theme::Fonts& fonts);
    void setPayload(std::shared_ptr<Payload> payload);
    PlotCanvasState captureState() const;
    std::optional<PhaseSpaceDomain> phaseSpaceDomain() const;
    void restoreState(PlotCanvasState state);
    void setMessage(std::wstring title, std::wstring detail = {});
    void nextPage();
    void previousPage();
    void setPage(int page);
    int page() const { return page_; }
    int pageCount() const;
    [[nodiscard]] std::wstring pageLabel(int page) const;
    [[nodiscard]] std::vector<AxisPanelState> currentAxisPanels() const;
    // False rejects invalid ranges, including nonpositive logarithmic Y bounds.
    bool setAxisLimits(std::size_t panel, std::optional<AxisLimits> limits);
    void appendPanel(Figure figure);
    void setBackgroundDragMode(PlotBackgroundDragMode mode);
    PlotBackgroundDragMode backgroundDragMode() const { return backgroundDragMode_; }
    void setValueScaling(const std::string& requested);
    std::string valueScaling() const;
    bool verticallyScrollable() const;
    PlotRenderStatistics renderStatistics() const;
    int scrollPosition() const { return scrollBar_.position(); }
    int scrollMaximum() const { return scrollMaximum_; }
    void setScrollPosition(int pixels);
    std::vector<std::size_t> visibleFigureIndices() const;
    bool isKrom() const { return payload_ && payload_->layout == "krom"; }
    std::optional<Payload::Animation> animation() const { return payload_ ? payload_->animation : std::nullopt; }
    const Payload* payloadIdentity() const { return payload_.get(); }
    std::string complexDisplay() const { return payload_ ? payload_->complexDisplay : ""; }
    void setComplexDisplay(const std::string& display);
    void selectKromModes(std::vector<std::size_t> modes);
    void selectAllKromModes(bool all);
    bool composed() const { return composed_; }
    bool hasPayload() const { return payload_ && !payload_->figures.empty(); }
    bool canDisplayCoordinates(const std::string& coordinates) const;
    std::string displayCoordinates() const;
    void setDisplayCoordinates(const std::string& coordinates);
    HWND handle() const { return window_; }

private:
    HWND window_ = nullptr;
    const theme::Fonts* fonts_ = nullptr;
    std::shared_ptr<Payload> payload_;
    std::wstring messageTitle_ = L"Select a scientific object";
    std::wstring messageDetail_ = L"Viewer loads only the plot data you request.";
    std::map<std::size_t, AxisLimits> axisLimits_;
    std::map<std::size_t, Plot3DCameraState> cameras3D_;
    int page_ = 0;
    OrnamentalScrollBar scrollBar_;
    RECT scrollViewport_{};
    int scrollMaximum_ = 0;
    float scrollRowHeight_ = 340;
    void updateScrollLayout();
    void wheelScroll(int delta);
    void advanceScroll();
    void stopScroll();
    static constexpr UINT_PTR scrollTimerId_ = 0x4453;
    static constexpr UINT_PTR cacheTimerId_ = 0x4454;
    void beginCacheFrame();
    void queueCachedPanel(std::size_t, const Gdiplus::RectF&);
    void warmPlotCache();
    void scheduleCacheWarmup();
    int scrollTarget_ = 0;
    ULONGLONG scrollTick_ = 0;
    struct RasterCache;
    std::shared_ptr<RasterCache> rasterCache_;
    std::size_t panelRenders_ = 0, cacheHits_ = 0;
    PlotTransform drawCachedPanel(Gdiplus::Graphics&, std::size_t, const Figure&,
        const Gdiplus::RectF&, const std::optional<AxisLimits>&,
        const Plot3DCameraState*, Plot3DControl);
    bool composed_ = false;
    std::string preferredCoordinates_;
    std::string preferredValueScale_;
    PlotBackgroundDragMode backgroundDragMode_ = PlotBackgroundDragMode::rearrange_panels;

    struct PanelHitRegion {
        std::size_t visualIndex = 0;
        std::size_t figureIndex = 0;
        RECT card{};
        RECT header{};
        RECT close{};
        RECT resize{};
        RECT plot{};
        PlotTransform transform{};
    };
    struct CanvasProbe {
        std::size_t figureIndex{};
        CurveProbe probe{};
    };
    enum class DragMode { none, reorder, resize, probe, pan, orbit3d, pan3d, zoom3d };
    std::vector<PanelHitRegion> hitRegions_;
    std::optional<std::size_t> hoveredPanel_;
    std::optional<CanvasProbe> hoverProbe_;
    std::map<std::size_t, CurveProbe> pinnedProbes_;
    std::optional<std::size_t> active3DFigure_;
    std::optional<std::size_t> hovered3DFigure_;
    Plot3DControl hovered3DControl_ = Plot3DControl::none;
    DragMode dragMode_ = DragMode::none;
    std::size_t draggedPanel_ = 0;
    std::size_t draggedFigure_ = 0;
    POINT dragOrigin_{};
    PlotTransform dragStartTransform_{};
    Plot3DCameraState dragStartCamera3D_{};
    bool dragConstrain3D_ = false;
    int dragStartColumnSpan_ = 1;
    int dragStartRowSpan_ = 1;
    int layoutColumns_ = 1;
    int layoutRows_ = 1;
    float layoutCellWidth_ = 1.0f;
    float layoutCellHeight_ = 1.0f;

    void paint(HDC target, bool interactive);
    void ensureComposition();
    void removePanel(std::size_t visualIndex);
    void pointerDown(POINT point, unsigned modifiers = 0, int button = 0);
    bool zoomPlotAt(POINT point, int wheelDelta, bool reset = false);
    bool key3D(WPARAM key, unsigned modifiers);
    bool selectKromAt(POINT point, bool control);
    void pointerMove(POINT point);
    void pointerUp();
    std::optional<std::size_t> hitPanel(POINT point) const;
    void notifyCompositionChanged();
    static LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam);
};

}  // namespace dsl
