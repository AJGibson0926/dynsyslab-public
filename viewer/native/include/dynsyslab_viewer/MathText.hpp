#pragma once

#include <string>
#include <vector>

namespace dsl {

enum class MathTextRole { normal, subscript, superscript, fraction };

struct MathTextSpan {
    std::wstring text;
    MathTextRole role{MathTextRole::normal};
    std::vector<MathTextSpan> numerator;
    std::vector<MathTextSpan> denominator;
};

[[nodiscard]] std::vector<MathTextSpan> parseMathText(
    const std::wstring& source);
[[nodiscard]] std::string scientificMathLabel(std::string source);

}  // namespace dsl
