#include "MathTextDrawing.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include <algorithm>
#include <cmath>
namespace dsl {
namespace {
struct Glyph { std::wstring text; float x{},y{},size{}; };
struct Rule { float x{},y{},width{}; };
struct Box {
    float width{},above{},below{};
    std::vector<Glyph> glyphs; std::vector<Rule> rules;
    void append(const Box& child,float x,float y) {
        for(auto glyph:child.glyphs) { glyph.x+=x;glyph.y+=y;glyphs.push_back(glyph); }
        for(auto rule:child.rules) { rule.x+=x;rule.y+=y;rules.push_back(rule); }
        above=std::max(above,child.above-y);below=std::max(below,child.below+y);
    }
};
Box layout(Gdiplus::Graphics& g,const std::vector<MathTextSpan>& spans,float size) {
    Box row;
    for(const auto& span:spans) {
        Box box;
        if(span.role==MathTextRole::fraction) {
            const auto numerator=layout(g,span.numerator,size*.83f);
            const auto denominator=layout(g,span.denominator,size*.83f);
            box.width=std::max(numerator.width,denominator.width)+size*.3f;
            const float axis=-size*.26f,gap=size*.16f;
            box.append(numerator,(box.width-numerator.width)/2,axis-gap-numerator.below);
            box.append(denominator,(box.width-denominator.width)/2,axis+gap+denominator.above);
            box.rules.push_back({0,axis,box.width});
        } else if(span.role!=MathTextRole::normal) {
            box=layout(g,span.numerator.empty()?parseMathText(span.text):span.numerator,size*.7f);
        } else {
            Gdiplus::Font font(L"Cambria Math",size,Gdiplus::FontStyleRegular,Gdiplus::UnitPixel);
            Gdiplus::RectF bounds;
            g.MeasureString(span.text.c_str(),-1,&font,Gdiplus::PointF{},
                Gdiplus::StringFormat::GenericTypographic(),&bounds);
            box.width=bounds.Width;box.above=size*.82f;box.below=size*.27f;
            box.glyphs.push_back({span.text,0,-box.above,size});
        }
        float offset=0;
        if(span.role==MathTextRole::subscript) offset=size*.26f;
        if(span.role==MathTextRole::superscript) offset=-size*.45f;
        row.append(box,row.width,offset);row.width+=box.width;
    }
    return row;
}
float fontSize(HFONT font) {
    LOGFONTW info{};GetObjectW(font,sizeof(info),&info);
    return std::max(14.0f,static_cast<float>(std::abs(info.lfHeight)));
}
}
float scientificTextScale(HFONT font) { return fontSize(font) / 14.f; }
Gdiplus::SizeF measureMathText(Gdiplus::Graphics& g,const std::wstring& text,HFONT font) {
    const auto box=layout(g,parseMathText(text),fontSize(font));
    return {box.width,box.above+box.below};
}
void drawMathText(Gdiplus::Graphics& g,const std::wstring& text,
    const Gdiplus::RectF& area,HFONT font,COLORREF color,Gdiplus::StringAlignment alignment) {
    const auto box=layout(g,parseMathText(text),fontSize(font));
    float x=area.X;
    if(alignment==Gdiplus::StringAlignmentCenter) x+=(area.Width-box.width)/2;
    if(alignment==Gdiplus::StringAlignmentFar) x+=area.Width-box.width;
    const float baseline=area.Y+(area.Height-box.above-box.below)/2+box.above;
    Gdiplus::SolidBrush brush(theme::color(color));
    Gdiplus::StringFormat format(Gdiplus::StringFormat::GenericTypographic());
    format.SetFormatFlags(Gdiplus::StringFormatFlagsNoWrap);
    for(const auto& glyph:box.glyphs) {
        Gdiplus::Font display(L"Cambria Math",glyph.size,Gdiplus::FontStyleRegular,Gdiplus::UnitPixel);
        g.DrawString(glyph.text.c_str(),-1,&display,Gdiplus::PointF{x+glyph.x,baseline+glyph.y},
            &format,&brush);
    }
    Gdiplus::Pen pen(theme::color(color),1.0f);
    for(const auto& rule:box.rules)
        g.DrawLine(&pen,x+rule.x,baseline+rule.y,x+rule.x+rule.width,baseline+rule.y);
}
}
