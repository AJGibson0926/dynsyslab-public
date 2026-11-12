#include "dynsyslab_compute/ScientificParameterFormatting.hpp"

#include <iostream>
#include <string>

namespace {

bool expect(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= expect(dsl::compute::scientificParameterLatex("equilibrium_radius_m")
        == L"R_{0}", "canonical RPE parameter should use mathematical notation");
    ok &= expect(dsl::compute::scientificParameterLatex("custom_gain_value")
        == L"custom_{gv}", "unknown names should retain a compact deterministic subscript");
    const auto scientific = dsl::compute::scientificParameterNumber(100000.0, false);
    if (scientific != L"1 \u00d7 10\u2075") {
        std::cerr << "scientific display code points:";
        for (const wchar_t value : scientific) std::cerr << ' ' << static_cast<unsigned>(value);
        std::cerr << '\n';
    }
    ok &= expect(scientific == L"1 \u00d7 10\u2075",
        "short display should use four-significant-figure scientific notation");
    ok &= expect(dsl::compute::scientificParameterNumber(1.0 / 3.0, false)
        == L"0.3333", "short display should retain four significant figures");
    ok &= expect(dsl::compute::scientificParameterUnit("dimensionless").empty(),
        "dimensionless groups should not reserve a unit column");
    ok &= expect(dsl::compute::scientificParameterUnit("m^3") == L"m\u00b3",
        "physical units should use compact mathematical exponents");
    double parsed = 0.0;
    ok &= expect(dsl::compute::parseScientificParameterInput(L"1e-3", parsed)
        && parsed == 1e-3, "C-style scientific notation should be editable");
    ok &= expect(dsl::compute::parseScientificParameterInput(L"1 x 10^{-3}", parsed)
        && parsed == 1e-3, "TeX-style scientific notation should be editable");
    ok &= expect(dsl::compute::parseScientificParameterInput(L"1 \u00d7 10\u207b\u00b3", parsed)
        && parsed == 1e-3, "rendered scientific notation should round-trip through editing");
    ok &= expect(dsl::compute::parseScientificParameterInput(L"2.5*10^(-4)", parsed)
        && parsed == 2.5e-4, "calculator-style scientific notation should be editable");
    ok &= expect(!dsl::compute::parseScientificParameterInput(L"1 x ten", parsed),
        "invalid scientific notation should be rejected");
    return ok ? 0 : 1;
}
