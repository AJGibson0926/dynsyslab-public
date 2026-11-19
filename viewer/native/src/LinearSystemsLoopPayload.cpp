#include "LinearSystemsLoopPayload.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace dsl {
namespace {

using dynsyslab::core::NpzRecord;
namespace json = dynsyslab::core::json;
using Cx = std::complex<double>;

constexpr const char* kBlue = "#004DB8";
constexpr const char* kRed = "#D8242A";
constexpr const char* kGold = "#DE9B00";
constexpr const char* kGrey = "#5B6270";
constexpr double kPi = 3.14159265358979323846;
constexpr double kNyquistRadius = 10.0;  // |L| beyond this is not drawn (zoom with axis ranges)
const double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string fmt(double v) { std::ostringstream t; t.precision(3); t << v; return t.str(); }
std::string plain(std::string s) { std::replace(s.begin(), s.end(), '_', ' '); return s; }  // '_' would subscript

struct Loop {
    std::string name, label, color, style;
    double opacity = 1.0;
    std::vector<double> w;
    std::vector<Cx> l;
    const json::Value* summary = nullptr;
};

std::string describe(const std::string& name, const json::Value& s) {
    if (name == "lqr_input") return "LQR, plant input";
    if (name == "kalman_output") return "Kalman filter, plant output";
    if (name == "lqg_input") return "LQG (KLQE->KLQR), plant input";
    if (name == "lqg_output") return "LQG (KLQE->KLQR), plant output";
    if (name.rfind("ltr_", 0) == 0) return "LTR q^2 = " + fmt(s.at("ltr_weight").number()) + ", plant input";
    return name;
}

std::vector<Loop> loops(const NpzRecord& r, const json::Value& meta) {
    std::vector<Loop> out;
    const auto& names = meta.at("loop_names").array();
    std::size_t ltr = 0, ltrCount = 0;
    for (const auto& n : names) if (n.string().rfind("ltr_", 0) == 0) ++ltrCount;
    for (const auto& n : names) {
        Loop loop; loop.name = n.string();
        const std::string prefix = "loop__" + loop.name + "__";
        if (!r.contains(prefix + "response")) continue;
        loop.summary = &meta.at("loop_robustness").at(loop.name);
        loop.w = r.at(prefix + "frequencies").doubles();
        const auto v = r.at(prefix + "response").doubles();
        for (std::size_t k = 0; k + 1 < v.size(); k += 2) loop.l.emplace_back(v[k], v[k + 1]);
        const double alpha = loop.summary->at("minimum_return_difference").number();
        loop.label = describe(loop.name, *loop.summary) + " (min|1+L| " + fmt(alpha) + ")";
        if (loop.name == "lqr_input") { loop.color = kBlue; loop.style = "solid"; }
        else if (loop.name == "kalman_output") { loop.color = kGold; loop.style = "solid"; }
        else if (loop.name == "lqg_input") { loop.color = kRed; loop.style = "solid"; }
        else if (loop.name == "lqg_output") { loop.color = kRed; loop.style = "dash"; }
        else { loop.color = kGold; loop.style = "dash"; loop.opacity = 0.35 + 0.6 * double(++ltr) / double(std::max<std::size_t>(ltrCount, 1)); }
        out.push_back(std::move(loop));
    }
    return out;
}

bool inputSide(const Loop& l) { return l.name.find("_input") != std::string::npos; }

Trace line(const std::string& label, const char* color, const std::string& style, double opacity = 1.0) {
    Trace t; t.label = label; t.color = color; t.lineStyle = style; t.opacity = opacity; return t;
}

Trace criticalPoint() {
    Trace t = line("-1", kGrey, "none"); t.marker = true; t.x = {-1.0}; t.y = {0.0}; return t;
}
Trace guaranteeDisk() {  // |1 + L| = 1: the LQR / Kalman loops never enter it
    Trace t = line("|1 + L| = 1 (LQR, Kalman guarantee)", kGrey, "dot");
    for (int k = 0; k <= 360; ++k) { const double a = 2 * kPi * k / 360; t.x.push_back(-1 + std::cos(a)); t.y.push_back(std::sin(a)); }
    return t;
}

Figure nyquistFigure(const std::vector<Loop>& all, bool input) {
    Figure f{input ? "nyquist_input" : "nyquist_output",
             std::string("Nyquist, loops broken at the plant ") + (input ? "input" : "output") +
             " (|L| > " + fmt(kNyquistRadius) + " not drawn; mirror image dotted)", "Re L(j omega)", "Im L(j omega)", true};
    f.traces.push_back(guaranteeDisk());
    for (const auto& loop : all) {
        if (inputSide(loop) != input) continue;
        Trace t = line(loop.label, loop.color.c_str(), loop.style, loop.opacity);
        Trace m = line("", loop.color.c_str(), "dot", 0.45 * loop.opacity);
        for (const auto& l : loop.l) {
            const bool keep = std::abs(l) <= kNyquistRadius;
            t.x.push_back(keep ? l.real() : kNaN); t.y.push_back(keep ? l.imag() : kNaN);
        }
        for (std::size_t k = loop.l.size(); k-- > 0;) {
            const auto l = std::conj(loop.l[k]);
            const bool keep = std::abs(l) <= kNyquistRadius;
            m.x.push_back(keep ? l.real() : kNaN); m.y.push_back(keep ? l.imag() : kNaN);
        }
        f.traces.push_back(std::move(t));
        if (!m.x.empty()) f.traces.push_back(std::move(m));
    }
    f.traces.push_back(criticalPoint());
    return f;
}

std::vector<double> log10w(const std::vector<double>& w) {
    std::vector<double> x;
    for (double v : w) x.push_back(v > 0 ? std::log10(v) : kNaN);
    return x;
}

std::vector<double> unwrap(const std::vector<Cx>& v) {
    std::vector<double> out; double prev = 0, offset = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const double raw = std::arg(v[i]) * 180 / kPi;
        if (i) { double step = raw + offset - prev; while (step > 180) { offset -= 360; step -= 360; } while (step < -180) { offset += 360; step += 360; } }
        prev = raw + offset; out.push_back(prev);
    }
    return out;
}

double db(double v) { return 20.0 * std::log10(std::max(v, 1e-300)); }

Payload bodePlant(Payload p, const NpzRecord& r, const json::Value& meta) {
    if (!r.contains("bode_magnitude_db")) throw std::runtime_error("archive has no Bode data (rerun linear_systems with inputs and outputs)");
    const auto x = log10w(r.at("bode_frequencies").doubles());
    const auto& mag = r.at("bode_magnitude_db"); const auto& ph = r.at("bode_phase_deg");
    const std::size_t channels = mag.shape.size() > 1 ? mag.shape[1] : 1;
    const auto m = mag.doubles(), a = ph.doubles();
    const char* colors[] = {kBlue, kRed, kGold};
    Figure fm{"bode_magnitude", "Bode magnitude of the plant channels", "log10 omega (rad per model time unit)", "|G| (dB)", false};
    Figure fp{"bode_phase", "Bode phase of the plant channels (unwrapped)", "log10 omega (rad per model time unit)", "phase (deg)", false};
    const auto& names = meta.at("bode_channels").array();
    const auto& tfs = meta.at("transfer_functions").array();
    for (std::size_t c = 0; c < channels; ++c) {
        std::string label = c < names.size() ? plain(names[c].string()) : "channel " + std::to_string(c);
        if (c < tfs.size())
            label += ", gain " + fmt(tfs[c].at("gain").number()) + ", relative degree " + fmt(tfs[c].at("relative_degree").number());
        Trace tm = line(label, colors[c % 3], c < 3 ? "solid" : "dash"), tp = tm;
        tm.x = x; tp.x = x;
        for (std::size_t k = 0; k < x.size(); ++k) { tm.y.push_back(m[k * channels + c]); tp.y.push_back(a[k * channels + c]); }
        fm.traces.push_back(std::move(tm)); fp.traces.push_back(std::move(tp));
    }
    p.figures.push_back(std::move(fm)); p.figures.push_back(std::move(fp));
    p.columns = 1;
    return p;
}

Payload loopBode(Payload p, const std::vector<Loop>& all) {
    Figure fm{"loop_bode_magnitude", "Open-loop magnitude |L(j omega)|", "log10 omega (rad per model time unit)", "|L| (dB)", false};
    Figure fp{"loop_bode_phase", "Open-loop phase (unwrapped)", "log10 omega (rad per model time unit)", "phase (deg)", false};
    for (const auto& loop : all) {
        Trace tm = line(loop.label, loop.color.c_str(), loop.style, loop.opacity), tp = tm;
        tm.x = tp.x = log10w(loop.w);
        for (const auto& l : loop.l) tm.y.push_back(db(std::abs(l)));
        tp.y = unwrap(loop.l);
        fm.traces.push_back(std::move(tm)); fp.traces.push_back(std::move(tp));
    }
    p.figures.push_back(std::move(fm)); p.figures.push_back(std::move(fp));
    return p;
}

Payload sensitivity(Payload p, const std::vector<Loop>& all) {
    Figure fs{"sensitivity", "Sensitivity |S| = 1/|1 + L| (peak = 1 / min|1 + L|; 0 dB line = LQR guarantee)",
              "log10 omega (rad per model time unit)", "|S| (dB)", false};
    Figure ft{"complementary_sensitivity", "Complementary sensitivity |T| = |L / (1 + L)|",
              "log10 omega (rad per model time unit)", "|T| (dB)", false};
    for (const auto& loop : all) {
        Trace ts = line(loop.label, loop.color.c_str(), loop.style, loop.opacity), tt = ts;
        ts.x = tt.x = log10w(loop.w);
        for (const auto& l : loop.l) { ts.y.push_back(db(1.0 / std::abs(1.0 + l))); tt.y.push_back(db(std::abs(l / (1.0 + l)))); }
        fs.traces.push_back(std::move(ts)); ft.traces.push_back(std::move(tt));
    }
    p.figures.push_back(std::move(fs)); p.figures.push_back(std::move(ft));
    return p;
}

}  // namespace

bool isClassicalLinearSystemsView(const std::string& view) {
    return view == "bode" || view == "nyquist" || view == "loop_bode" || view == "sensitivity";
}

Payload buildClassicalLinearSystemsPayload(const std::filesystem::path& source, const std::string& view,
                                           const NpzRecord& record, const json::Value& meta) {
    Payload p{source.u8string(), view, meta.at("label").string(source.stem().u8string()), 1};
    if (view == "bode") return bodePlant(std::move(p), record, meta);
    if (meta.at("loop_names").array().empty())
        throw std::runtime_error("archive has no loop analysis (give linear_systems LQR and/or LQE weights with a single input)");
    const auto all = loops(record, meta);
    if (view == "nyquist") {
        p.figures.push_back(nyquistFigure(all, true));
        if (std::any_of(all.begin(), all.end(), [](const Loop& l) { return !inputSide(l); }))
            p.figures.push_back(nyquistFigure(all, false));
        p.columns = int(p.figures.size());
        return p;
    }
    if (view == "loop_bode") return loopBode(std::move(p), all);
    if (view == "sensitivity") return sensitivity(std::move(p), all);
    throw std::runtime_error("unknown classical linear-systems view: " + view);
}

}  // namespace dsl
