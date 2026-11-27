#pragma once

#include "dynsyslab_viewer/Json.hpp"

#include <complex>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace dsl {
// Default size of plot text relative to the application UI text.
inline constexpr double defaultPlotFontScale = 2.3;


struct Relation {
    std::string target;
    std::string kind;
    std::string label;
};

struct CatalogItem {
    std::string id;
    std::string kind;
    std::string category;
    std::string label;
    std::string source;
    std::string system;
    std::string subtitle;
    double modified = 0.0;
    double bytes = 0.0;
    std::vector<std::string> views;
    std::vector<std::string> tags;
    std::map<std::string, std::string> properties;
    std::vector<Relation> relations;
    std::string role;
    std::string parentId;
    std::string lineageId;
    bool primary = true;
    int variantCount = 0;
};

struct Catalog {
    std::string projectName;
    std::string projectRoot;
    std::string focusLabel;
    std::size_t hiddenItemCount = 0;
    std::vector<CatalogItem> items;
    std::map<std::string, std::string> viewLabels;

    static Catalog load(const std::wstring& path);
    static Catalog fromJson(const std::string& text);
    [[nodiscard]] std::string toJson() const;
    static Catalog scan(const std::filesystem::path& projectRoot,
                        const std::filesystem::path& archiveRoot);
    static Catalog indexed(const std::filesystem::path& projectRoot,
                           const std::filesystem::path& archiveRoot);
    int indexOf(const std::string& id) const;
};

struct Trace {
    std::string label;
    std::string role;
    std::string color;
    std::string lineStyle;
    std::string encoding;
    bool marker = false;
    bool heatmap = false;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> colorValues;
    std::vector<double> magnitudes;
    std::vector<int> zShape;
    int selectionId = -1;
    bool selected = false;
    // Keep newly added presentation fields at the end so existing aggregate
    // initializers retain their long-standing field order.
    double opacity = 1.0;
    double lineWidthPoints = 0.0; // Zero selects the role's standard pen.
    // Marker shape when `marker` is set: "" = small filled dot, "open-circle" = hollow ring.
    std::string markerShape;
    // A complex heat map keeps its values so the complex display can be switched in place
    // (ComplexHeatmap.hpp): which split panel this is (0 or 1, -1 for none) and the function name.
    std::vector<std::complex<double>> complexValues;
    int complexPart = -1;
    std::string complexName;
    // Fixed colour range of a scalar heat map (e.g. one range for every frame of an animation);
    // empty: the range of this trace's values.
    std::optional<std::pair<double, double>> colorLimits;
    // Brightness reference of a domain-coloured line (brightness = magnitude / scale). Zero: the
    // trace's own largest magnitude. Overlaid runs of one function share one scale.
    double magnitudeScale = 0.0;
};

struct DetailInset {
    double xMin = 0.0;
    double xMax = 0.0;
    double yMin = 0.0;
    double yMax = 0.0;
    std::vector<std::size_t> traceIndices;
    bool automatic = false;
    // Transient policy: when a short burst (e.g. a start-up spike) would flatten the
    // whole panel, the inset shows that burst at full scale and the main panel keeps
    // the meaningful (robust) ordinate range [mainYMin, mainYMax] instead.
    bool transientWindow = false;
    double mainYMin = 0.0;
    double mainYMax = 0.0;
    // Ordinate of the inset itself: symmetric log with this linear core when the
    // episode spans decades beyond the settled level (e.g. a start-up spike).
    bool symmetricLog = false;
    double linearThreshold = 1.0;
    double decadeWeight = 1.0;
};

struct Figure {
    std::string kind;
    std::string title;
    std::string xLabel;
    std::string yLabel;
    bool square = false;
    std::string titleColor;
    int dashboardPage = -1;
    int gridSpan = 1;
    int rowSpan = 1;
    std::vector<Trace> traces;
    std::optional<DetailInset> inset;
    struct CoordinateAxes {
        std::string name, xLabel, yLabel;
        double xScale = 1, yScale = 1, zScale = 1, timeScale = 1;
        std::vector<std::string> traceLabels;
        std::string colorLabel;
        std::string zLabel;
    };
    struct ReferencePoint {
        double x = 0, y = 0;
        std::string label;
        bool equilibrium = false;
    };
    std::string coordinates;
    std::string colorLabel;
    std::optional<CoordinateAxes> alternateCoordinates;
    std::vector<ReferencePoint> references;
    std::string valueScale = "auto";
    std::string observableIdentity;
    std::string zLabel;
    // Optional initial 3-D camera for "three_d" figures; reset returns here.
    std::optional<double> cameraAzimuthDegrees;
    std::optional<double> cameraElevationDegrees;
    // "three_d" figures that are not state-space trajectories (e.g. stem plots):
    // scale each axis to the box independently and omit the trajectory IC marker.
    bool boxAxes3D = false;
    // Index abscissas (singular-value index, mode number): ticks and grid at integers.
    bool integerXTicks = false;
    // Optional right ordinate for a curve plotted as offset + scale * value.
    struct RightAxis { std::string label; double scale = 1, offset = 0; };
    std::optional<RightAxis> rightAxis;
    // Initial abscissa range (e.g. a spectrum's dominant band); the ordinate is fitted inside it.
    // Resetting the axes returns to the full data range.
    std::optional<std::pair<double, double>> focusX;
};

struct Payload {
    std::string source;
    std::string view;
    std::string tileName;
    int columns = 1;
    std::vector<std::string> pageLabels;
    std::vector<std::string> warnings;
    std::vector<Figure> figures;
    std::string layout;
    std::string complexDisplay;
    std::vector<std::size_t> selectedModes;
    // An animation: its pages are frames with the same panels, played in order by the plot
    // toolbar (first, back, play/pause, forward, stop).
    struct Animation {
        double frameSeconds = 0.4;
        bool loop = true;
        bool autoplay = true;
    };
    std::optional<Animation> animation;

    static Payload load(const std::wstring& path);
};

std::optional<DetailInset> automaticDetailInset(
    const Figure& figure, double disparityThreshold = 6.0);

std::wstring wide(const std::string& value);
std::string narrow(const std::wstring& value);
std::string canonicalSystemLabel(const std::string& value);

}  // namespace dsl
