// dynsyslab.research-figure/v1 documents (DocumentPayload.hpp): presentation of supplied research
// results. {"schema":"dynsyslab.research-figure/v1","figure":KIND, ...}:
//   cycle_spectrum_sequence  records:[{label, atoms:[angle], eigenvalues:{method:[[re,im]]}}]
//   return_sections          records:[{label, states:[[...]]}], coordinates:[i,j], labels:[x,y]
//   prediction_comparison    record:{label, truth:[...], predictions:{method:[...]}, observable}, samples
//   regularized_measures     records:[{label, measures:{method:{grid,density}}, atoms}]
//   multiplier_branch        branch, local: {parameter:[], multiplier:[], bracket:[a,b]}, parameter_label
//   perturbation             records:[{label, time, reference_deviation:[[..]], models:[{label,
//                            deviation, error_norm, color, linestyle}]}], coordinates
//   truncation_amplitude     records:[{label, quotients, spectral_radii, perturbations, rms_errors}]
//   cyclic_partition         record:{label, states, geometry:{levels:[{depth, cells, maximum_diameter,
//                            minimum_sample_separation, singleton_terminal_level}]}}, depth
//   coordinate_reconstruction records:[{label, time, reference:[[x1,x2]], decoded, ode?}], labels
//   character_mass           records:[{label, conductor_masses:[{denominator, mass}]}]
//   error_comparison         x:[...], curves:{label:[...]}, x_label, y_label, title
//   clock_character          evolution:{samples:[{time, bare_characters, corrected_characters}]},
//                            linear_reference:[[re,im]], character
//   solenoid_suspension      base, depth, phase_turns, address: the cut-open suspension of the
//                            base-adic adding machine (flow segments glued at the seam) and the Cantor
//                            refinement of one fiber; a finite-prefix schematic, not a phase portrait
//   cycle_clock              periods:[2,4,8] (a divisibility chain), initial_step, steps: finite quotient
//                            clocks of the adding machine, one animation frame per step (exact integers)
// Eigenangles are in turns (angle / 2 pi) per step. Nothing is fitted or certified here.
#include "DocumentPayload.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <set>
#include <stdexcept>

namespace dsl {
namespace {
namespace json = dynsyslab::core::json;
constexpr double kTwoPi = 6.28318530717958647692;
constexpr const char* kRed = "#C52B24";
constexpr const char* kBlue = "#0072B2";
constexpr const char* kGold = "#B57B00";

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.isNull() ? std::nan("") : item.number());
    return out;
}
std::vector<double> column(const json::Value& rows, std::size_t c) {
    std::vector<double> out;
    for (const auto& row : rows.array()) out.push_back(row.array().size() > c ? row.array()[c].number() : std::nan(""));
    return out;
}
std::vector<double> range(std::size_t n) { std::vector<double> out(n); for (std::size_t i = 0; i < n; ++i) out[i] = double(i); return out; }
std::string methodColor(const std::string& method) { return method == "resdmd" ? kRed : method == "mpedmd" ? kBlue : kGold; }

Trace line(std::string label, std::vector<double> x, std::vector<double> y, std::string color, std::string style = "solid") {
    Trace t; t.label = std::move(label); t.x = std::move(x); t.y = std::move(y); t.color = std::move(color); t.lineStyle = std::move(style);
    return t;
}
Trace points(std::string label, std::vector<double> x, std::vector<double> y, std::string color, bool open = false) {
    auto t = line(std::move(label), std::move(x), std::move(y), std::move(color)); t.marker = true;
    if (open) t.markerShape = "open-circle";
    return t;
}
Trace horizontal(const std::string& label, double y, double x0, double x1, const std::string& color) {
    return line(label, {x0, x1}, {y, y}, color, "dash");
}
std::pair<double, double> span(const std::vector<double>& x) {
    double lo = INFINITY, hi = -INFINITY;
    for (const double v : x) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
    return {lo, hi};
}
Figure figure(std::string kind, std::string title, std::string x, std::string y, bool square = false) {
    Figure f; f.kind = std::move(kind); f.title = std::move(title); f.xLabel = std::move(x); f.yLabel = std::move(y); f.square = square;
    return f;
}

void cycleFigures(const json::Value& d, const std::string& kind, Payload& p) {
    if (kind == "cycle_spectrum_sequence") {
        auto f = figure("phase", d.at("title").string("Finite-cycle Koopman eigenvalues"),
                        "Eigenangle / (2 pi), per return-map step", "record");
        const auto& records = d.at("records").array();
        for (std::size_t j = 0; j < records.size(); ++j) {
            std::vector<double> a;
            for (const double atom : numbers(records[j].at("atoms"))) a.push_back(atom / kTwoPi);
            f.traces.push_back(points(j ? "" : "finite-cycle reference", a, std::vector<double>(a.size(), double(j)), "#737373", true));
            double offset = -0.14;
            for (const auto& [method, pairs] : records[j].at("eigenvalues").object()) {
                std::vector<double> x;
                for (const auto& z : pairs.array()) {
                    double t = std::fmod(std::atan2(z.array()[1].number(), z.array()[0].number()) / kTwoPi + 1.0, 1.0);
                    x.push_back(std::abs(t - 1.0) < 1e-8 ? 0.0 : t);
                }
                f.traces.push_back(points(j ? "" : method, x, std::vector<double>(x.size(), j + offset), methodColor(method)));
                offset = 0.14;
            }
        }
        f.integerXTicks = false; f.focusX = std::make_pair(-0.025, 1.025);
        p.figures.push_back(std::move(f));
    } else if (kind == "return_sections") {
        const auto c = numbers(d.at("coordinates")); const auto i = c.empty() ? 0 : std::size_t(c[0]), k = c.size() < 2 ? 1 : std::size_t(c[1]);
        for (const auto& record : d.at("records").array()) {
            auto f = figure("phase", record.at("label").string(), d.at("labels").array().size() > 0 ? d.at("labels").array()[0].string() : "x_1",
                            d.at("labels").array().size() > 1 ? d.at("labels").array()[1].string() : "x_2", true);
            auto t = points("stroboscopic states", column(record.at("states"), i), column(record.at("states"), k), "");
            t.colorValues = range(t.x.size()); t.encoding = "time";
            f.colorLabel = "sample index"; f.traces.push_back(std::move(t)); p.figures.push_back(std::move(f));
        }
        p.columns = std::max<int>(1, static_cast<int>(p.figures.size()));
    } else if (kind == "prediction_comparison") {
        const auto& r = d.at("record");
        auto truth = numbers(r.at("truth"));
        const auto n = std::min<std::size_t>(static_cast<std::size_t>(d.at("samples").number(128)), truth.size());
        truth.resize(n);
        auto a = figure("time", r.at("label").string(), "Return-map steps from initialization", r.at("observable").string("Observable"));
        auto b = figure("time", "Prediction - truth", "Return-map steps from initialization", "Prediction - truth");
        a.traces.push_back(line("Held-out trajectory", range(n), truth, kGold));
        for (const auto& [method, values] : r.at("predictions").object()) {
            auto v = numbers(values); v.resize(n);
            std::vector<double> e(n); for (std::size_t s = 0; s < n; ++s) e[s] = v[s] - truth[s];
            a.traces.push_back(line(method, range(n), v, methodColor(method), "dash"));
            b.traces.push_back(line(method, range(n), e, methodColor(method)));
        }
        p.figures.push_back(std::move(a)); p.figures.push_back(std::move(b));
    } else {  // regularized_measures
        for (const auto& record : d.at("records").array()) {
            auto f = figure("fourier", record.at("label").string(), "Spectral angle / (2 pi), per return-map step", "Signed density");
            f.valueScale = "symlog";
            for (const auto& [method, measure] : record.at("measures").object()) {
                std::vector<double> x; for (const double g : numbers(measure.at("grid"))) x.push_back(g / kTwoPi);
                f.traces.push_back(line(method, x, numbers(measure.at("density")), methodColor(method)));
            }
            p.figures.push_back(std::move(f));
        }
    }
}

void cocycleFigures(const json::Value& d, const std::string& kind, Payload& p) {
    if (kind == "multiplier_branch") {
        for (const char* key : {"branch", "local"}) {
            const auto& r = d.at(key);
            const auto x = numbers(r.at("parameter")); const auto [x0, x1] = span(x);
            auto f = figure("time", key == std::string("branch") ? "Physical periodic-orbit branch" : "Refined crossing bracket",
                            d.at("parameter_label").string("Parameter"), "Real transverse multiplier");
            f.traces.push_back(line("multiplier", x, numbers(r.at("multiplier")), kBlue));
            f.traces.push_back(horizontal("flip candidate: multiplier = -1", -1.0, x0, x1, kRed));
            const auto bracket = numbers(r.at("bracket"));
            if (bracket.size() == 2) f.references.push_back({bracket[0], -1.0, "bracket start"}), f.references.push_back({bracket[1], -1.0, "bracket end"});
            p.figures.push_back(std::move(f));
        }
    } else if (kind == "perturbation") {
        const auto c = numbers(d.at("coordinates"));
        const std::vector<std::size_t> comps = c.size() >= 2 ? std::vector<std::size_t>{std::size_t(c[0]), std::size_t(c[1])} : std::vector<std::size_t>{0, 1};
        for (const auto& record : d.at("records").array()) {
            const auto t = numbers(record.at("time"));
            for (const auto comp : comps) {
                auto f = figure("time", record.at("label").string() + " | dx" + std::to_string(comp + 1), "Forcing periods since perturbation",
                                "dx" + std::to_string(comp + 1));
                f.traces.push_back(line("independent nonlinear ODE", t, column(record.at("reference_deviation"), comp), kRed));
                for (const auto& model : record.at("models").array())
                    f.traces.push_back(line(model.at("label").string(), t, column(model.at("deviation"), comp),
                                            model.at("color").string(kBlue), model.at("linestyle").string("dash")));
                p.figures.push_back(std::move(f));
            }
            auto e = figure("time", record.at("label").string() + " | prediction error", "Forcing periods since perturbation",
                            "Deviation prediction error");
            e.valueScale = "log";
            for (const auto& model : record.at("models").array()) {
                auto err = numbers(model.at("error_norm"));
                if (!err.empty()) err[0] = std::nan("");  // initialization is not an out-of-sample forecast
                for (auto& v : err) if (!(v > 0)) v = std::nan("");
                e.traces.push_back(line(model.at("label").string(), t, err, model.at("color").string(kBlue), model.at("linestyle").string("dash")));
            }
            p.figures.push_back(std::move(e));
        }
    } else if (kind == "truncation_amplitude") {
        auto a = figure("phase", "Averaging the cocycle can create false instability", "Retained quotient period M (log10)",
                        "Spectral radius over the original period N");
        auto b = figure("phase", "Full cocycle: limits of the linear approximation", "Initial angle displacement (log10)",
                        "RMS deviation prediction error");
        a.valueScale = b.valueScale = "log";
        constexpr const char* colors[] = {kBlue, kGold, kRed, "#009E73"};
        std::size_t i = 0;
        for (const auto& record : d.at("records").array()) {
            const auto logs = [](std::vector<double> v) { for (auto& x : v) x = x > 0 ? std::log10(x) : std::nan(""); return v; };
            auto ta = line(record.at("label").string(), logs(numbers(record.at("quotients"))), numbers(record.at("spectral_radii")), colors[i % 4]);
            auto tb = line(record.at("label").string(), logs(numbers(record.at("perturbations"))), numbers(record.at("rms_errors")), colors[i % 4]);
            a.traces.push_back(std::move(ta)); b.traces.push_back(std::move(tb)); ++i;
        }
        if (!a.traces.empty()) { const auto [x0, x1] = span(a.traces.front().x); a.traces.push_back(horizontal("stability boundary", 1.0, x0, x1, "#555555")); }
        p.figures.push_back(std::move(a)); p.figures.push_back(std::move(b));
    } else {  // cyclic_partition
        const auto& r = d.at("record"); const auto& levels = r.at("geometry").at("levels").array();
        const auto depth = static_cast<std::size_t>(d.at("depth").number(3));
        if (depth >= levels.size()) throw std::invalid_argument("depth must be an available partition level");
        const auto cells = static_cast<std::size_t>(levels[depth].at("cells").number());
        auto a = figure("phase", r.at("label").string() + " | residues modulo " + std::to_string(cells), "x_1", "x_2", true);
        auto t = points("states", column(r.at("states"), 0), column(r.at("states"), 1), "");
        for (std::size_t s = 0; s < t.x.size(); ++s) t.colorValues.push_back(double(s % std::max<std::size_t>(1, cells)));
        t.encoding = "time"; a.colorLabel = "time-index residue"; a.traces.push_back(std::move(t));
        auto b = figure("time", "Shrinking finite cells (not an infinite-limit proof)", "Binary depth n", "Euclidean distance");
        b.valueScale = "log"; b.integerXTicks = true;
        std::vector<double> dx, diameter, sx, separation;
        for (const auto& level : levels) {
            if (level.at("singleton_terminal_level").boolean()) continue;
            dx.push_back(level.at("depth").number()); diameter.push_back(level.at("maximum_diameter").number());
            if (!level.at("minimum_sample_separation").isNull()) { sx.push_back(level.at("depth").number()); separation.push_back(level.at("minimum_sample_separation").number()); }
        }
        b.traces.push_back(line("largest sampled cell diameter", dx, diameter, kBlue));
        b.traces.push_back(line("smallest sampled intercell distance", sx, separation, kGold, "dash"));
        p.figures.push_back(std::move(a)); p.figures.push_back(std::move(b));
    }
}

void reconstructionFigures(const json::Value& d, const std::string& kind, Payload& p) {
    const auto label = [&](std::size_t i) { return d.at("labels").array().size() > i ? d.at("labels").array()[i].string() : "x_" + std::to_string(i + 1); };
    if (kind == "coordinate_reconstruction") {
        p.columns = 3;
        for (const auto& r : d.at("records").array()) {
            const auto t = numbers(r.at("time"));
            for (std::size_t c = 0; c < 2; ++c) {
                auto f = figure("time", r.at("label").string(), "Time", label(c));
                if (!r.at("ode").isNull()) f.traces.push_back(line("independent ODE", t, column(r.at("ode"), c), kGold));
                f.traces.push_back(line("closed form", t, column(r.at("reference"), c), kRed));
                f.traces.push_back(line("character reconstruction", t, column(r.at("decoded"), c), kBlue, "dash"));
                p.figures.push_back(std::move(f));
            }
            auto f = figure("phase", r.at("label").string(), label(0), label(1), true);
            f.traces.push_back(line("closed form", column(r.at("reference"), 0), column(r.at("reference"), 1), kRed));
            f.traces.push_back(line("character reconstruction", column(r.at("decoded"), 0), column(r.at("decoded"), 1), kBlue, "dash"));
            p.figures.push_back(std::move(f));
        }
    } else if (kind == "character_mass") {
        std::set<double> orders;
        const auto& records = d.at("records").array();
        for (const auto& r : records) for (const auto& g : r.at("conductor_masses").array()) if (g.at("denominator").number() > 1) orders.insert(g.at("denominator").number());
        if (orders.empty()) throw std::invalid_argument("no nonconstant character mass supplied");
        const std::vector<double> rows(orders.begin(), orders.end());
        Trace heat; heat.heatmap = true; heat.label = "log10 squared coefficient mass";
        heat.zShape = {static_cast<int>(rows.size()), static_cast<int>(records.size())};
        heat.z.assign(rows.size() * records.size(), std::nan(""));
        for (std::size_t j = 0; j < records.size(); ++j)
            for (const auto& g : records[j].at("conductor_masses").array()) {
                const auto at = std::find(rows.begin(), rows.end(), g.at("denominator").number());
                if (at != rows.end() && g.at("mass").number() > 0)
                    heat.z[std::size_t(at - rows.begin()) * records.size() + j] = std::log10(g.at("mass").number());
            }
        heat.x = range(records.size()); heat.y = rows;
        auto f = figure("heatmap", d.at("title").string("Character contributions across finite cycles"), "cycle (record index)",
                        "character order / reduced denominator");
        f.colorLabel = "log10 squared coefficient mass (constant part omitted)"; f.integerXTicks = true;
        f.traces.push_back(std::move(heat)); p.figures.push_back(std::move(f));
    } else if (kind == "error_comparison") {
        auto f = figure("phase", d.at("title").string("Reconstruction diagnostics"), d.at("x_label").string("Resolution") + " (log10)",
                        d.at("y_label").string("Error (stated norm)"));
        f.valueScale = "log";
        constexpr const char* colors[] = {kRed, kBlue, kGold, "#00865A"};
        std::size_t i = 0;
        for (const auto& [name, values] : d.at("curves").object()) {
            std::vector<double> x, y; const auto xs = numbers(d.at("x")); const auto ys = numbers(values);
            for (std::size_t k = 0; k < xs.size() && k < ys.size(); ++k) if (ys[k] > 0 && xs[k] > 0) { x.push_back(std::log10(xs[k])); y.push_back(ys[k]); }
            auto t = line(name, x, y, colors[i++ % 4]); f.traces.push_back(std::move(t));
        }
        p.figures.push_back(std::move(f));
    } else {  // clock_character
        const auto& evolution = d.at("evolution").at("result").isNull() ? d.at("evolution") : d.at("evolution").at("result");
        const auto character = static_cast<std::size_t>(d.at("character").number(0));
        std::vector<double> time, bare, corrected, reference, error;
        const auto& rows = evolution.at("samples").array(); const auto& linear = d.at("linear_reference").array();
        if (rows.empty() || linear.size() != rows.size()) throw std::invalid_argument("the linear reference must use the clock sample mesh");
        for (std::size_t k = 0; k < rows.size(); ++k) {
            const auto pair = [&](const json::Value& v) { return std::complex<double>(v.array()[0].number(), v.array()[1].number()); };
            const auto b = pair(rows[k].at("bare_characters").array().at(character)), c = pair(rows[k].at("corrected_characters").array().at(character));
            const auto r = pair(linear[k]);
            time.push_back(rows[k].at("time").number()); bare.push_back(b.real()); corrected.push_back(c.real());
            reference.push_back(r.real()); error.push_back(std::abs(c - r));
        }
        auto a = figure("time", d.at("title").string("Characters under a change of speed"), "Physical time", "Real part");
        a.traces.push_back(line("original character", time, bare, kGold));
        a.traces.push_back(line("corrected eigenfunction", time, corrected, kRed));
        a.traces.push_back(line("constant-frequency evolution", time, reference, kBlue, "dash"));
        auto b = figure("time", "Complex prediction error", "Physical time", "|corrected - reference|");
        b.traces.push_back(line("error", time, error, kBlue));
        p.figures.push_back(std::move(a)); p.figures.push_back(std::move(b));
    }
}

// Solenoid and odometer schematics (finite prefixes; the transverse placement is illustrative).
void schematicFigures(const json::Value& d, const std::string& kind, Payload& p) {
    if (kind == "solenoid_suspension") {
        const int base = static_cast<int>(d.at("base").number(2)), depth = static_cast<int>(d.at("depth").number(3));
        const double phase = d.at("phase_turns").number(0.4);
        const long long address = static_cast<long long>(d.at("address").number(3));
        if (base < 2 || depth < 1 || depth > 5 || std::pow(base, depth) > 32 || !(phase >= 0 && phase < 1))
            throw std::invalid_argument("solenoid_suspension needs base >= 2, 1 <= depth, base^depth <= 32 and phase in [0, 1)");
        const int n = static_cast<int>(std::lround(std::pow(base, depth))), radix = 2 * base - 1;
        const auto interval = [&](long long a, int levels) {
            double y = 0; long long power = 1;
            for (int j = 0; j < levels; ++j, power *= base) y += 2.0 * ((a / power) % base) / std::pow(radix, j + 1);
            return std::make_pair(y, y + std::pow(radix, -levels));
        };
        std::vector<double> ys;
        for (int a = 0; a < n; ++a) { const auto [lo, hi] = interval(a, depth); ys.push_back(0.5 * (lo + hi)); }
        const auto chosen = static_cast<int>(((address % n) + n) % n), successor = (chosen + 1) % n;
        auto flow = figure("phase", "Flow segments; the two edges are glued", "base phase / 2 pi (time runs to the right)",
                           "address modulo " + std::to_string(n) + " (schematic placement)");
        Trace rows = line("cylinders", {}, {}, "#AEB9BF");
        for (const double y : ys) { rows.x.insert(rows.x.end(), {0.0, 1.0, std::nan("")}); rows.y.insert(rows.y.end(), {y, y, std::nan("")}); }
        flow.traces.push_back(std::move(rows));
        flow.traces.push_back(points("one fiber", std::vector<double>(ys.size(), phase), ys, kRed));
        flow.traces.push_back(line("seam: right edge at address " + std::to_string(chosen) + " continues at the left edge at "
                                   + std::to_string(successor), {phase, 1.0, std::nan(""), 0.0, 0.22},
                                   {ys[chosen], ys[chosen], std::nan(""), ys[successor], ys[successor]}, kBlue));
        flow.traces.push_back(points("glued points", {1.0, 0.0}, {ys[chosen], ys[successor]}, kBlue));
        auto fiber = figure("phase", "Fiber: Cantor refinement", "refinement depth", "schematic transverse placement");
        Trace cells = line("cylinder sets", {}, {}, kRed); cells.lineWidthPoints = 6.0;
        for (int level = 1; level <= depth + 1; ++level)
            for (long long a = 0; a < static_cast<long long>(std::pow(base, level)); ++a) {
                const auto [lo, hi] = interval(a, level);
                cells.x.insert(cells.x.end(), {double(level), double(level), std::nan("")});
                cells.y.insert(cells.y.end(), {lo, hi, std::nan("")});
            }
        fiber.traces.push_back(std::move(cells)); fiber.integerXTicks = true;
        p.figures.push_back(std::move(flow)); p.figures.push_back(std::move(fiber));
        p.tileName = std::to_string(base) + "-adic solenoid, cut open (finite-prefix schematic)";
    } else {  // cycle_clock
        std::vector<long long> periods;
        for (const auto& period : d.at("periods").array()) periods.push_back(static_cast<long long>(period.number()));
        if (periods.empty() || periods.size() > 6) throw std::invalid_argument("cycle_clock needs one to six periods");
        for (std::size_t k = 0; k < periods.size(); ++k)
            if (periods[k] < 2 || periods[k] > 32 || (k && (periods[k] <= periods[k - 1] || periods[k] % periods[k - 1])))
                throw std::invalid_argument("cycle_clock periods must be an increasing divisibility chain in [2, 32]");
        const auto start = static_cast<long long>(d.at("initial_step").number(0));
        const auto steps = static_cast<int>(d.at("steps").number(static_cast<double>(periods.back())));
        for (int frame = 0; frame < steps; ++frame) {
            const long long step = start + frame;
            std::string label = "step " + std::to_string(step) + ": residues";
            for (std::size_t k = 0; k < periods.size(); ++k) {
                const long long n = periods[k], residue = ((step % n) + n) % n;
                label += " " + std::to_string(residue) + "/" + std::to_string(n);
                auto clock = figure("phase", "period " + std::to_string(n) + ": residue " + std::to_string(residue),
                                    "cos(2 pi r / n)", "sin(2 pi r / n)", true);
                clock.dashboardPage = frame;
                Trace circle = line("", {}, {}, "#AEB9BF");
                for (int a = 0; a <= 120; ++a) { circle.x.push_back(std::cos(kTwoPi * a / 120)); circle.y.push_back(std::sin(kTwoPi * a / 120)); }
                clock.traces.push_back(std::move(circle));
                std::vector<double> cx, cy;
                for (long long r = 0; r < n; ++r) { cx.push_back(std::cos(kTwoPi * double(r) / double(n))); cy.push_back(std::sin(kTwoPi * double(r) / double(n))); }
                clock.traces.push_back(points("residues", cx, cy, "#17232E", true));
                clock.traces.push_back(points("current", {cx[residue]}, {cy[residue]}, kRed));
                p.figures.push_back(std::move(clock));
            }
            p.pageLabels.push_back(label);
        }
        p.columns = static_cast<int>(periods.size());
        p.animation = Payload::Animation{};
        p.tileName = "Finite quotient clocks of the adding machine";
    }
}

}  // namespace

Payload buildResearchFigurePayload(const json::Value& document) {
    const auto kind = document.at("figure").string();
    Payload payload{{}, kind, document.at("title").string(kind), 2};
    if (kind == "cycle_spectrum_sequence" || kind == "return_sections" || kind == "prediction_comparison" || kind == "regularized_measures")
        cycleFigures(document, kind, payload);
    else if (kind == "multiplier_branch" || kind == "perturbation" || kind == "truncation_amplitude" || kind == "cyclic_partition")
        cocycleFigures(document, kind, payload);
    else if (kind == "coordinate_reconstruction" || kind == "character_mass" || kind == "error_comparison" || kind == "clock_character")
        reconstructionFigures(document, kind, payload);
    else if (kind == "solenoid_suspension" || kind == "cycle_clock") schematicFigures(document, kind, payload);
    else throw std::invalid_argument("unknown research figure: " + kind);
    if (payload.figures.empty()) throw std::invalid_argument("the research figure document has no records");
    return payload;
}

}  // namespace dsl
