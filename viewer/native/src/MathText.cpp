#include "dynsyslab_viewer/MathText.hpp"
#include <cwctype>
#include <cctype>
#include <map>
namespace dsl {
namespace {
std::wstring symbol(const std::wstring& name) {
    static const std::map<std::wstring,std::wstring> values{
        {L"alpha",L"α"},{L"beta",L"β"},{L"gamma",L"γ"},{L"delta",L"δ"},
        {L"Delta",L"Δ"},{L"kappa",L"κ"},{L"lambda",L"λ"},{L"mu",L"μ"},
        {L"nu",L"ν"},{L"xi",L"ξ"},{L"Xi",L"Ξ"},{L"rho",L"ρ"},{L"sigma",L"σ"},{L"tau",L"τ"},
        {L"theta",L"\u03b8"},{L"Theta",L"\u0398"},
        {L"omega",L"ω"},{L"chi",L"χ"},{L"Pi",L"Π"},{L"psi",L"ψ"},
        {L"Psi",L"Ψ"},{L"phi",L"φ"},{L"nabla",L"∇"},{L"times",L"×"},
        {L"cdot",L"·"},{L"infty",L"∞"},{L"partial",L"∂"}};
    const auto found=values.find(name);
    return found==values.end()?name:found->second;
}
struct Parser {
    const std::wstring& input;
    std::size_t at=0;
    std::wstring group() {
        if(at>=input.size()) return {};
        if(input[at]!=L'{') return input.substr(at++,1);
        const auto start=++at; int depth=1;
        while(at<input.size()&&depth) {
            if(input[at]==L'{') ++depth;
            if(input[at]==L'}') --depth;
            if(depth) ++at;
        }
        const auto result=input.substr(start,at-start);
        if(at<input.size()) ++at;
        return result;
    }
    std::vector<MathTextSpan> parse() {
        std::vector<MathTextSpan> result;
        auto text=[&](std::wstring value) {
            if(!result.empty()&&result.back().role==MathTextRole::normal)
                result.back().text+=value;
            else result.push_back({std::move(value)});
        };
        while(at<input.size()) {
            const auto ch=input[at++];
            if(ch==L'$'||ch==L'{'||ch==L'}') continue;
            if(ch==L'_'||ch==L'^') {
                auto children=parseMathText(group()); std::wstring plain;
                for(const auto& child:children) plain+=child.text;
                result.push_back({plain,ch==L'_'?MathTextRole::subscript:
                    MathTextRole::superscript,std::move(children)});
            } else if(ch==L'\\') {
                const auto start=at;
                while(at<input.size()&&std::iswalpha(input[at])) ++at;
                const auto name=input.substr(start,at-start);
                if(name==L"frac") {
                    const auto top=group(), bottom=group();
                    result.push_back({{},MathTextRole::fraction,
                        parseMathText(top),parseMathText(bottom)});
                } else if(name==L"mathrm"||name==L"text"||name==L"operatorname") text(group());
                else if(name==L"dot") {
                    for(const auto& part:parseMathText(group())) text(part.text);
                    text(L"\u0307");
                } else if(name.empty()) {
                    if(at<input.size()) { const auto escaped=input[at++];
                        text(escaped==L','||escaped==L';'||escaped==L' '?L" ":
                            std::wstring(1,escaped)); }
                } else text(symbol(name));
            } else text(std::wstring(1,ch));
        }
        return result;
    }
};
}
std::vector<MathTextSpan> parseMathText(const std::wstring& source) {
    return Parser{source}.parse();
}
std::string scientificMathLabel(std::string source) {
    const std::pair<const char*,const char*> replacements[]{
        {" (dimensionless)",""},{" (dimensional)",""},
        {"dr/d\\tau","\\frac{dr}{d\\tau}"},{"dr/dtau","\\frac{dr}{d\\tau}"},
        {"dx/dt","\\frac{dx}{dt}"},{"dy/dt","\\frac{dy}{dt}"},{"dz/dt","\\frac{dz}{dt}"},
        {"dR/dt","\\frac{dR}{dt}"},{"R/R_{0}","\\frac{R}{R_{0}}"},
        {"R/R0","\\frac{R}{R_{0}}"}};
    for(const auto& entry:replacements) {
        std::size_t at=0;
        while((at=source.find(entry.first,at))!=std::string::npos) {
            source.replace(at,std::char_traits<char>::length(entry.first),entry.second);
            at+=std::char_traits<char>::length(entry.second);
        }
    }
    std::size_t at = 0;
    while ((at = source.find("tau", at)) != std::string::npos) {
        const bool left = at == 0 || (!std::isalnum(static_cast<unsigned char>(source[at-1])) && source[at-1] != '\\');
        const bool right = at + 3 == source.size() || !std::isalnum(static_cast<unsigned char>(source[at+3]));
        if (left && right) { source.replace(at, 3, "\\tau"); at += 4; }
        else at += 3;
    }
    return source;
}
}
