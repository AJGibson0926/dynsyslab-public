#pragma once
#include "dynsyslab_viewer/Theme.hpp"
#include <algorithm>
#include <cmath>
namespace dsl {
// Cyclic hue encodes argument; brightness is sqrt(|psi| / the mode's shared maximum).
inline Gdiplus::Color domainColor(double phase, double magnitude = 1) {
    if (!std::isfinite(phase) || !std::isfinite(magnitude)) return theme::color(theme::plotBackground);
    constexpr double pi = 3.14159265358979323846;
    const double hue = (phase+pi)/(2*pi)*6, value = std::sqrt(std::clamp(magnitude,0.,1.));
    const double x = 1-std::abs(std::fmod(hue,2.)-1);
    double r=0,g=0,b=0;
    switch (static_cast<int>(std::floor(hue))%6) {
        case 0: r=1; g=x; break; case 1: r=x; g=1; break;
        case 2: g=1; b=x; break; case 3: g=x; b=1; break;
        case 4: r=x; b=1; break; default: r=1; b=x; break;
    }
    return Gdiplus::Color(255,static_cast<BYTE>(255*r*value),static_cast<BYTE>(255*g*value),static_cast<BYTE>(255*b*value));
}
}
