#include "dynsyslab_compute/SolverExecutionPane.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <utility>

namespace dsl::compute {
namespace {

void selectText(HWND combo, const std::string& value) {
    const auto wanted = dynsyslab::core::wide(value);
    LRESULT index = SendMessageW(combo, CB_FINDSTRINGEXACT,
        static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(wanted.c_str()));
    if (index == CB_ERR) {
        SendMessageW(combo, CB_ADDSTRING, 0,
            reinterpret_cast<LPARAM>(wanted.c_str()));
        index = SendMessageW(combo, CB_FINDSTRINGEXACT,
            static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(wanted.c_str()));
    }
    if (index != CB_ERR) SendMessageW(combo, CB_SETCURSEL, index, 0);
}

void setValue(HWND window, const dsl::json::Value& value,
              const char* fallback) {
    const std::string rendered = value.isNull()
        ? fallback : dsl::json::dump(value);
    SetWindowTextW(window, dynsyslab::core::wide(rendered).c_str());
}

}  // namespace

bool SolverExecutionPane::create(HWND parent, const dsl::theme::Fonts& fonts) {
    parent_ = parent;
    auto label = [&](const wchar_t* text) {
        HWND result = makeLabel(parent, text, fonts.caption);
        labels_.push_back(result); controls_.push_back(result); return result;
    };
    section_ = makeLabel(parent, L"Solver and execution", fonts.heading);
    controls_.push_back(section_);
    label(L"Method"); method_ = makeCombo(parent, 0, fonts.body);
    label(L"Internal step"); internalStep_ = makeInlineValue(parent, 0, L"0", fonts.body);
    label(L"rtol"); relativeTolerance_ = makeInlineValue(parent, 0, L"1e-8", fonts.body);
    label(L"atol"); absoluteTolerance_ = makeInlineValue(parent, 0, L"1e-10", fonts.body);
    label(L"Maximum step"); maximumStep_ = makeInlineValue(parent, 0, L"0", fonts.body);
    label(L"Maximum steps"); maximumInternalSteps_ = makeInlineValue(parent, 0, L"10000000", fonts.body);
    label(L"Backend"); backend_ = makeCombo(parent, 0, fonts.body);
    for (HWND control : {method_, internalStep_, relativeTolerance_,
         absoluteTolerance_, maximumStep_, maximumInternalSteps_, backend_})
        controls_.push_back(control);
    setComboItems(method_, {"RK45", "BDF", "RK4", "Euler",
                            "method_of_steps", "metropolis",
                            "RK23", "DOP853", "Radau", "Adams", "LSODA"});
    setComboItems(backend_, {"auto", "native"});
    resetForFamily("ode");
    return true;
}

void SolverExecutionPane::show(bool visible, bool enabled) {
    for (HWND control : controls_) {
        ShowWindow(control, visible ? SW_SHOW : SW_HIDE);
        EnableWindow(control, enabled ? TRUE : FALSE);
    }
}

void SolverExecutionPane::layout(int x, int y, int width, int height) {
    const int gap = scalePixels(parent_, 7);
    const int labelH = scalePixels(parent_, 15);
    const int headerH = scalePixels(parent_, 19);
    const int rowH = scalePixels(parent_, 27);
    auto field = [&](std::size_t label, HWND control, int left, int top,
                     int fieldWidth, int captionWidth) {
        MoveWindow(labels_.at(label), left, top + (rowH - labelH) / 2,
            captionWidth, labelH, TRUE);
        moveFieldControl(control, left + captionWidth, top,
            std::max(1, fieldWidth - captionWidth), rowH);
    };
    MoveWindow(section_, x, y, width, headerH, TRUE);
    int top = y + headerH + scalePixels(parent_, 3);
    const int first = (width - 2 * gap) / 3;
    field(0, method_, x, top, first, scalePixels(parent_, 50));
    field(6, backend_, x + first + gap, top, first,
        scalePixels(parent_, 48));
    field(1, internalStep_, x + 2 * (first + gap), top,
        width - 2 * (first + gap), scalePixels(parent_, 78));
    top += rowH + gap;
    const int quarter = (width - 3 * gap) / 4;
    field(2, relativeTolerance_, x, top, quarter, scalePixels(parent_, 32));
    field(3, absoluteTolerance_, x + quarter + gap, top, quarter,
        scalePixels(parent_, 32));
    field(4, maximumStep_, x + 2 * (quarter + gap), top, quarter,
        scalePixels(parent_, 76));
    field(5, maximumInternalSteps_, x + 3 * (quarter + gap), top,
        width - 3 * (quarter + gap), scalePixels(parent_, 82));
    (void)height;
}

void SolverExecutionPane::resetForFamily(const std::string& family) {
    const std::string method = family == "pde" ? "RK4"
        : family == "dde" ? "method_of_steps"
        : family == "lattice" ? "metropolis" : "RK45";
    baseSettings_ = "{}";
    selectText(method_, method); selectText(backend_, "auto");
    SetWindowTextW(internalStep_, family == "pde" ? L"0.01"
        : family == "dde" ? L"0.05" : L"0");
    SetWindowTextW(relativeTolerance_, L"1e-8");
    SetWindowTextW(absoluteTolerance_, L"1e-10");
    SetWindowTextW(maximumStep_, L"0");
    SetWindowTextW(maximumInternalSteps_, L"10000000");
}

void SolverExecutionPane::load(const std::string& solverJson) {
    auto root = dsl::json::parse(jsonObject(solverJson));
    root.object().erase("control_evaluation");
    root.object().erase("control_sample_time");
    root.object().erase("parallel");
    baseSettings_ = dsl::json::dump(root);
    selectText(method_, root.at("method").string("RK45"));
    setValue(internalStep_, root.at("step_size"), "0");
    setValue(relativeTolerance_, root.at("rtol"), "1e-8");
    setValue(absoluteTolerance_, root.at("atol"), "1e-10");
    setValue(maximumStep_, root.at("max_step"), "0");
    setValue(maximumInternalSteps_, root.at("max_internal_steps"), "10000000");
    selectText(backend_, root.at("execution_backend").string("auto"));
}

std::string SolverExecutionPane::settings() const {
    std::vector<std::pair<std::string, std::string>> values{
        {"method", "\"" + jsonEscape(selectedText(method_)) + "\""},
        {"step_size", jsonObject(utf8Text(internalStep_), "0")},
        {"rtol", jsonObject(utf8Text(relativeTolerance_), "1e-8")},
        {"atol", jsonObject(utf8Text(absoluteTolerance_), "1e-10")},
        {"max_step", jsonObject(utf8Text(maximumStep_), "0")},
        {"max_internal_steps", jsonObject(utf8Text(maximumInternalSteps_), "10000000")},
        {"execution_backend", "\"" + jsonEscape(selectedText(backend_)) + "\""}};
    return mergeSettingsObject(baseSettings_, values);
}

void SolverExecutionPane::focus() {
    if (method_) SetFocus(method_);
}

}  // namespace dsl::compute
