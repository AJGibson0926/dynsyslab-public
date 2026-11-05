#pragma once
#include "dynsyslab_viewer/Theme.hpp"
#include <gdiplus.h>
namespace dsl {
float scientificTextScale(HFONT);
Gdiplus::SizeF measureMathText(Gdiplus::Graphics&,const std::wstring&,HFONT);
void drawMathText(Gdiplus::Graphics&,const std::wstring&,const Gdiplus::RectF&,
    HFONT,COLORREF,Gdiplus::StringAlignment=Gdiplus::StringAlignmentNear);
}
