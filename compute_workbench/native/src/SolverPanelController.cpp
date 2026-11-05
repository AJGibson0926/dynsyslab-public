#include "dynsyslab_compute/SolverPanel.hpp"

#include "ControllerSettingsDialog.hpp"
#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_compute/WaveformPreview.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace dsl::compute {
namespace {
std::string baseTargetKind(std::string value) {
    if(value.rfind("observable from state ",0)==0)value.erase(0,22);
    else if(value.rfind("observable ",0)==0)value.erase(0,11);
    return value;
}


int choice(HWND combo) {
    return static_cast<int>(SendMessageW(combo, CB_GETCURSEL, 0, 0));
}

void visible(HWND window, bool show) {
    if (!!IsWindowVisible(window) != show)
        ShowWindow(window, show ? SW_SHOW : SW_HIDE);
}

void fieldVisible(HWND label, HWND field, bool show) {
    visible(label, show); visible(field, show);
}

std::string arrayLiteral(std::string value, const char* fallback = "[]") {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return fallback;
    if (value[begin] != '[') value = '[' + value + ']';
    const auto parsed = dsl::json::parse(value);
    if (!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))
        throw std::invalid_argument("Expected a comma-separated list of values.");
    return dsl::json::dump(parsed);
}

std::string scalarOrArray(std::string value, const char* fallback) {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return fallback;
    try { return dsl::json::dump(dsl::json::parse(value)); }
    catch (...) { return arrayLiteral(value, fallback); }
}

std::vector<double> numbers(std::string value, double fallback) {
    const auto parsed = dsl::json::parse(arrayLiteral(std::move(value)));
    std::vector<double> result;
    for (const auto& item : parsed.array()) result.push_back(item.number());
    if (result.empty()) result.push_back(fallback);
    return result;
}

std::string compact(const std::vector<double>& values) {
    std::ostringstream result;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) result << ", ";
        result << values[index];
    }
    return result.str();
}

std::string zeros(std::size_t dimension) {
    std::vector<double> values(std::max<std::size_t>(1, dimension), 0.0);
    return '[' + compact(values) + ']';
}

std::string display(const dsl::json::Value& value, const char* fallback) {
    if (value.isNull()) return fallback;
    std::string encoded = dsl::json::dump(value);
    if (encoded.size() >= 2 && encoded.front() == '[' && encoded.back() == ']')
        encoded = encoded.substr(1, encoded.size() - 2);
    return encoded;
}

double numeric(HWND window, double fallback) {
    try { return std::stod(utf8Text(window)); }
    catch (...) { return fallback; }
}

}  // namespace

void SolverPanel::changeControllerKind() {
    const std::string kind = selectedText(controller_);
    if (kind == "none") controlPageActive_ = false;
    SetWindowTextW(controllerSettings_, L"{}");
    if (kind == "waveform") {
        SetWindowTextW(waveformFrequencies_, L"1");
        SetWindowTextW(waveformAmplitudes_, L"1");
        SetWindowTextW(waveformPhases_, L"0");
        SetWindowTextW(waveformOffset_, L"0");
        SendMessageW(waveformCoordinates_, CB_SETCURSEL, 0, 0);
        updateWaveformSettings();
    }
    updatePanelPageVisibility();
    if (layoutWidth_ > 0)
        layoutEditor(layoutX_, layoutY_, layoutWidth_, layoutHeight_);
    RECT upper{layoutX_, layoutY_, layoutX_ + layoutWidth_,
        layoutY_ + layoutHeight_};
    RedrawWindow(parent_, &upper, nullptr,
        RDW_INVALIDATE | RDW_ALLCHILDREN);
}

void SolverPanel::updateWaveformControls() {
    try {
        const auto settings = dsl::json::parse(jsonObject(utf8Text(controllerSettings_)));
        const bool model = !settings.at("model_amplitudes").isNull()
            || !settings.at("model_amplitude").isNull()
            || !settings.at("frequencies").isNull();
        SendMessageW(waveformCoordinates_, CB_SETCURSEL, model ? 1 : 0, 0);
        const auto& frequencies = settings.at(model ? "frequencies" : "frequencies_hz");
        const auto& amplitudes = settings.at(model ? "model_amplitudes" : "amplitudes");
        SetWindowTextW(waveformFrequencies_, dsl::wide(display(frequencies, "1")).c_str());
        SetWindowTextW(waveformAmplitudes_, dsl::wide(display(amplitudes, "1")).c_str());
        SetWindowTextW(waveformPhases_, dsl::wide(display(settings.at("phases_rad"), "0")).c_str());
        SetWindowTextW(waveformOffset_, dsl::wide(display(
            settings.at(model ? "model_offset" : "offset"), "0")).c_str());
    } catch (...) {
        SetWindowTextW(waveformFrequencies_, L"1");
        SetWindowTextW(waveformAmplitudes_, L"1");
        SetWindowTextW(waveformPhases_, L"0");
        SetWindowTextW(waveformOffset_, L"0");
    }
    setWaveformPreview(waveformPreview_, utf8Text(controllerSettings_),
        numeric(start_, 0.0), numeric(stop_, 10.0));
}

void SolverPanel::updateWaveformSettings() {
    const bool model = choice(waveformCoordinates_) == 1;
    std::vector<std::pair<std::string, std::string>> values{
        {model ? "frequencies" : "frequencies_hz",
            arrayLiteral(utf8Text(waveformFrequencies_))},
        {model ? "model_amplitudes" : "amplitudes",
            arrayLiteral(utf8Text(waveformAmplitudes_))},
        {"phases_rad", arrayLiteral(utf8Text(waveformPhases_))},
        {model ? "model_offset" : "offset",
            jsonObject(utf8Text(waveformOffset_), "0")},
        {model ? "frequencies_hz" : "frequencies", "null"},
        {model ? "amplitudes" : "model_amplitudes", "null"},
        {model ? "offset" : "model_offset", "null"},
        {"frequency", "null"}, {"frequency_hz", "null"},
        {"amplitude", "null"}, {"model_amplitude", "null"}};
    try {
        const std::string settings = mergeSettingsObject(
            jsonObject(utf8Text(controllerSettings_)), values);
        SetWindowTextW(controllerSettings_, dsl::wide(settings).c_str());
        setWaveformPreview(waveformPreview_, settings,
            numeric(start_, 0.0), numeric(stop_, 10.0));
    } catch (const std::exception& error) {
        setWaveformPreview(waveformPreview_,
            "{\"frequencies\":[],\"model_amplitudes\":[]}", 0.0, 1.0);
        MessageBoxW(parent_, dsl::wide(error.what()).c_str(),
            L"Invalid waveform", MB_OK | MB_ICONWARNING);
    }
}

void SolverPanel::adjustWaveformComponents(bool add) {
    try {
        auto frequencies = numbers(utf8Text(waveformFrequencies_), 1.0);
        auto amplitudes = numbers(utf8Text(waveformAmplitudes_), 1.0);
        auto phases = numbers(utf8Text(waveformPhases_), 0.0);
        const std::size_t count = std::max({frequencies.size(), amplitudes.size(), phases.size()});
        frequencies.resize(count, frequencies.back()); amplitudes.resize(count, 0.0);
        phases.resize(count, 0.0);
        if (add) {
            frequencies.push_back(frequencies.back() + 1.0);
            amplitudes.push_back(0.0); phases.push_back(0.0);
        } else if (count > 1) {
            frequencies.pop_back(); amplitudes.pop_back(); phases.pop_back();
        }
        SetWindowTextW(waveformFrequencies_, dsl::wide(compact(frequencies)).c_str());
        SetWindowTextW(waveformAmplitudes_, dsl::wide(compact(amplitudes)).c_str());
        SetWindowTextW(waveformPhases_, dsl::wide(compact(phases)).c_str());
        updateWaveformSettings();
    } catch (const std::exception& error) {
        MessageBoxW(parent_, dsl::wide(error.what()).c_str(),
            L"Invalid waveform components", MB_OK | MB_ICONWARNING);
    }
}

void SolverPanel::updateTargetControls() {
    updateControllerControls();
}

void SolverPanel::updateControllerControls() {
    const std::string type = selectedText(controller_);
    const bool waveform = type == "waveform";
    const bool lqr = type == "lqr";
    const bool koopman = type == "klqr" || type == "kmpc";
    const bool riccati = lqr || type == "sdre" || type == "mpc" || koopman;
    const bool predictive = type == "mpc" || type == "kmpc";
    const bool targeted = type == "pid" || type == "sdre"
        || type == "mpc" || koopman;
    const bool feedback = type != "none" && type != "waveform"
        && type != "random_pulse";
    const bool configured = type != "none" && !lqr;

    const std::vector<HWND> conditional{
        controllerEmpty_, targetSection_, waveformPreview_, waveformCoordinatesLabel_,
        waveformCoordinates_, waveformFrequenciesLabel_, waveformFrequencies_,
        waveformAmplitudesLabel_, waveformAmplitudes_, waveformPhasesLabel_,
        waveformPhases_, waveformOffsetLabel_, waveformOffset_, waveformAdd_,
        waveformRemove_, lqrALabel_, lqrA_, lqrBLabel_, lqrB_, labels_[14], krom_,
        labels_[16], target_, labels_[17], targetState_, labels_[18], q_, labels_[19], r_,
        labels_[20], representation_, labels_[21], horizon_, labels_[22], controlLimit_,
        controlEvaluationLabel_, controlEvaluation_, controlSampleTimeLabel_,
        controlSampleTime_,
        controllerSettingsButton_, targetKindLabel_, targetKind_, targetSourceLabel_,
        targetSource_, targetFrequenciesLabel_, targetFrequencies_,
        targetAmplitudesLabel_, targetAmplitudes_, targetPhasesLabel_, targetPhases_};
    for (HWND window : conditional) visible(window, false);
    if (!editorVisible_ || !controlPageActive_ || !hasControl()) return;
    visible(controllerEmpty_, type == "none");
    visible(controllerSettingsButton_, configured);
    fieldVisible(labels_[14], krom_, koopman);
    fieldVisible(labels_[18], q_, riccati);
    fieldVisible(labels_[19], r_, riccati);
    fieldVisible(labels_[20], representation_, koopman);
    fieldVisible(labels_[21], horizon_, predictive);
    fieldVisible(labels_[22], controlLimit_, predictive);
    fieldVisible(controlEvaluationLabel_, controlEvaluation_, feedback);
    fieldVisible(controlSampleTimeLabel_, controlSampleTime_, feedback
        && selectedText(controlEvaluation_) == "sampled");
    fieldVisible(lqrALabel_, lqrA_, lqr);
    fieldVisible(lqrBLabel_, lqrB_, lqr);
    for (HWND window : {waveformPreview_, waveformCoordinatesLabel_, waveformCoordinates_,
         waveformFrequenciesLabel_, waveformFrequencies_, waveformAmplitudesLabel_,
         waveformAmplitudes_, waveformPhasesLabel_, waveformPhases_,
         waveformOffsetLabel_, waveformOffset_, waveformAdd_, waveformRemove_})
        visible(window, waveform);
    visible(targetSection_, targeted);
    fieldVisible(targetKindLabel_, targetKind_, targeted);
    if (targeted) {
        const std::string targetKind = baseTargetKind(selectedText(targetKind_));
        fieldVisible(labels_[17], targetState_, targetKind == "fixed point"
            || targetKind == "waveform");
        fieldVisible(targetSourceLabel_, targetSource_, targetKind == "trajectory");
        fieldVisible(labels_[16], target_, targetKind == "canonical preset");
        const bool harmonic = targetKind == "waveform";
        fieldVisible(targetFrequenciesLabel_, targetFrequencies_, harmonic);
        fieldVisible(targetAmplitudesLabel_, targetAmplitudes_, harmonic);
        fieldVisible(targetPhasesLabel_, targetPhases_, harmonic);
    }
    if (waveform) updateWaveformControls();
    if (layoutWidth_ > 0)
        layoutEditor(layoutX_, layoutY_, layoutWidth_, layoutHeight_);
}

void SolverPanel::layoutControllerArea(int x, int y, int width, int height) {
    const int inner = scalePixels(parent_, 6), labelH = scalePixels(parent_, 15);
    const int headerH = scalePixels(parent_, 19), rowH = scalePixels(parent_, 27);
    auto place = [&](HWND label, HWND field, int px, int py, int pw, int captionW) {
        MoveWindow(label, px, py + (rowH - labelH) / 2, captionW, labelH, TRUE);
        moveFieldControl(field, px + captionW, py, pw - captionW, rowH);
    };
    auto placeEvaluation = [&](int px, int py, int pw) {
        const bool sampled = selectedText(controlEvaluation_) == "sampled";
        const int evaluationW = sampled ? (pw - inner) * 63 / 100 : pw;
        place(controlEvaluationLabel_, controlEvaluation_, px, py, evaluationW,
            scalePixels(parent_, 104));
        if (sampled)
            place(controlSampleTimeLabel_, controlSampleTime_,
                px + evaluationW + inner, py, pw - evaluationW - inner,
                scalePixels(parent_, 62));
    };
    MoveWindow(controllerSection_, x, y, width, headerH, TRUE);
    MoveWindow(controllerSeparator_, x, y + height, width, 1, TRUE);
    const int top = y + headerH + scalePixels(parent_, 3);
    const std::string type = selectedText(controller_);
    const bool waveform = type == "waveform", lqr = type == "lqr";
    const bool koopman = type == "klqr" || type == "kmpc";
    const bool predictive = type == "mpc" || type == "kmpc";
    const bool targeted = type == "pid" || type == "sdre" || type == "mpc" || koopman;
    const int selectorW = width * 28 / 100;
    place(labels_[13], controller_, x, top, selectorW, scalePixels(parent_, 58));

    if (type == "none") {
        MoveWindow(controllerEmpty_, x, top + rowH + inner, width,
            height - headerH - rowH - inner, TRUE); return;
    }
    if (waveform) {
        place(labels_[13], controller_, x, top, width,
            scalePixels(parent_, 66));
        const int coordinatesY = top + rowH + inner;
        place(waveformCoordinatesLabel_, waveformCoordinates_, x,
            coordinatesY, width, scalePixels(parent_, 66));

        const int bottom = y + height - scalePixels(parent_, 2);
        const int buttonsY = bottom - rowH;
        const int offsetY = buttonsY - inner - rowH;
        const int phasesY = offsetY - inner - rowH;
        const int amplitudesY = phasesY - inner - rowH;
        const int frequenciesY = amplitudesY - inner - rowH;
        const int previewY = coordinatesY + rowH + inner;
        MoveWindow(waveformPreview_, x, previewY, width,
            std::max(scalePixels(parent_, 96), frequenciesY - inner - previewY), TRUE);
        place(waveformFrequenciesLabel_, waveformFrequencies_, x,
            frequenciesY, width, scalePixels(parent_, 24));
        place(waveformAmplitudesLabel_, waveformAmplitudes_, x,
            amplitudesY, width, scalePixels(parent_, 24));
        place(waveformPhasesLabel_, waveformPhases_, x,
            phasesY, width, scalePixels(parent_, 24));
        place(waveformOffsetLabel_, waveformOffset_, x,
            offsetY, width, scalePixels(parent_, 24));
        const int third = (width - 2 * inner) / 3;
        MoveWindow(waveformAdd_, x, buttonsY, third, rowH, TRUE);
        MoveWindow(waveformRemove_, x + third + inner, buttonsY,
            third, rowH, TRUE);
        MoveWindow(controllerSettingsButton_, x + 2 * (third + inner),
            buttonsY, width - 2 * (third + inner), rowH, TRUE); return;
    }
    if (lqr) {
        const int half = (width - inner) / 2;
        place(lqrALabel_, lqrA_, x, top + rowH + inner, half,
            scalePixels(parent_, 22));
        place(lqrBLabel_, lqrB_, x + half + inner, top + rowH + inner,
            width - half - inner, scalePixels(parent_, 22));
        place(labels_[18], q_, x, top + 2 * (rowH + inner), half,
            scalePixels(parent_, 22));
        place(labels_[19], r_, x + half + inner, top + 2 * (rowH + inner),
            width - half - inner, scalePixels(parent_, 22));
        placeEvaluation(x, top + 3 * (rowH + inner), width); return;
    }

    const int settingsX = x + selectorW + inner;
    if (koopman) {
        const int kromW = width * 43 / 100;
        place(labels_[14], krom_, settingsX, top, kromW, scalePixels(parent_, 62));
        MoveWindow(controllerSettingsButton_, settingsX + kromW + inner, top,
            width - selectorW - kromW - 2 * inner, rowH, TRUE);
    } else {
        MoveWindow(controllerSettingsButton_, settingsX, top,
            width - selectorW - inner, rowH, TRUE);
    }
    const int bodyY = top + rowH + inner;
    const int targetX = targeted ? x + width * 52 / 100 : x + width;
    const int controlsW = targeted ? targetX - x - inner : width;
    int next = bodyY;
    if (type == "pid") {
        // PID gains and active window live in the structured settings editor.
    } else {
        const int half = (controlsW - inner) / 2;
        place(labels_[18], q_, x, bodyY, half, scalePixels(parent_, 20));
        place(labels_[19], r_, x + half + inner, bodyY,
            controlsW - half - inner, scalePixels(parent_, 20));
        next = bodyY + rowH + inner;
        if (koopman) {
            place(labels_[20], representation_, x, next, controlsW,
                scalePixels(parent_, 54)); next += rowH + inner;
        }
        if (predictive) {
            place(labels_[21], horizon_, x, next, half, scalePixels(parent_, 58));
            place(labels_[22], controlLimit_, x + half + inner, next,
                controlsW - half - inner, scalePixels(parent_, 40));
            next += rowH + inner;
        }
    }
    placeEvaluation(x, next, controlsW);
    if (!targeted) return;
    const int targetW = x + width - targetX;
    MoveWindow(targetSection_, targetX, bodyY, targetW, headerH, TRUE);
    const int targetTop = bodyY + headerH + scalePixels(parent_, 2);
    place(targetKindLabel_, targetKind_, targetX, targetTop, targetW,
        scalePixels(parent_, 36));
    const int detailY = targetTop + rowH + inner;
    const std::string kind = baseTargetKind(selectedText(targetKind_));
    if (kind == "fixed point")
        place(labels_[17], targetState_, targetX, detailY, targetW,
            scalePixels(parent_, 25));
    else if (kind == "trajectory")
        place(targetSourceLabel_, targetSource_, targetX, detailY, targetW,
            scalePixels(parent_, 62));
    else if (kind == "canonical preset")
        place(labels_[16], target_, targetX, detailY, targetW,
            scalePixels(parent_, 48));
    else {
        place(labels_[17], targetState_, targetX, detailY, targetW,
            scalePixels(parent_, 25));
        const int third = (targetW - 2 * inner) / 3;
        place(targetFrequenciesLabel_, targetFrequencies_, targetX,
            detailY + rowH + inner, third, scalePixels(parent_, 18));
        place(targetAmplitudesLabel_, targetAmplitudes_, targetX + third + inner,
            detailY + rowH + inner, third, scalePixels(parent_, 18));
        place(targetPhasesLabel_, targetPhases_, targetX + 2 * (third + inner),
            detailY + rowH + inner, targetW - 2 * (third + inner),
            scalePixels(parent_, 18));
    }
}

void SolverPanel::editControllerSettings() {
    const auto result = showControllerSettingsDialog(parent_, selectedText(controller_),
        utf8Text(controllerSettings_), *fonts_);
    if (result) {
        SetWindowTextW(controllerSettings_, dsl::wide(*result).c_str());
        if (selectedText(controller_) == "waveform") updateWaveformControls();
    }
}

std::string SolverPanel::controlSettings() const {
    const std::string type = selectedText(controller_);
    std::string settings = jsonObject(utf8Text(controllerSettings_));
    if (type == "waveform" || type == "random_pulse" || type == "pid"
        || type == "none") return settings;
    const bool riccati = type == "lqr" || type == "sdre" || type == "mpc"
        || type == "klqr" || type == "kmpc";
    std::vector<std::pair<std::string, std::string>> values;
    if (riccati) {
        values.push_back({"Q", scalarOrArray(utf8Text(q_), "1.0")});
        values.push_back({"R", scalarOrArray(utf8Text(r_), "1.0")});
    }
    if (type == "lqr") {
        if (!utf8Text(lqrA_).empty()) values.push_back({"A", arrayLiteral(utf8Text(lqrA_))});
        if (!utf8Text(lqrB_).empty()) values.push_back({"B", arrayLiteral(utf8Text(lqrB_))});
    }
    const bool koopman = type == "klqr" || type == "kmpc";
    if (koopman) {
        const int model = choice(krom_);
        const std::string source = model > 0
            && static_cast<std::size_t>(model - 1) < visibleKroms_.size()
            ? registry_->kroms[visibleKroms_[static_cast<std::size_t>(model - 1)]].source : "";
        values.insert(values.begin(), {"krom", "\"" + jsonEscape(source) + "\""});
        values.push_back({"control_representation",
            "\"" + jsonEscape(selectedText(representation_)) + "\""});
    }
    if (type == "mpc" || type == "kmpc") {
        values.push_back({"horizon", jsonObject(utf8Text(horizon_), "20")});
        if (!utf8Text(controlLimit_).empty())
            values.push_back({"control_limit", scalarOrArray(utf8Text(controlLimit_), "null")});
    }
    return mergeSettingsObject(settings, values);
}

std::string SolverPanel::targetSettingsBase() const {
    const std::string controller = selectedText(controller_);
    if (controller == "none" || controller == "waveform"
        || controller == "random_pulse") return "null";
    if (controller == "lqr") {
        const std::size_t index = systemIndex();
        const std::size_t dimension = index < registry_->systems.size()
            ? static_cast<std::size_t>(std::max(1, registry_->systems[index].stateDimension))
            : 1;
        return "{\"kind\":\"constant\",\"state\":" + zeros(dimension) + "}";
    }
    const std::string kind = baseTargetKind(selectedText(targetKind_));
    if (kind == "trajectory") {
        const int selected = choice(targetSource_);
        const std::string source = selected > 0
            && static_cast<std::size_t>(selected - 1) < visibleTargetTrajectories_.size()
            ? registry_->trajectories[visibleTargetTrajectories_[
                static_cast<std::size_t>(selected - 1)]].source : "";
        return "{\"kind\":\"trajectory\",\"source\":\""
            + jsonEscape(source) + "\"}";
    }
    if (kind == "waveform")
        return "{\"kind\":\"waveform\",\"offset\":"
            + arrayLiteral(utf8Text(targetState_)) + ",\"frequencies\":"
            + arrayLiteral(utf8Text(targetFrequencies_)) + ",\"amplitudes\":"
            + arrayLiteral(utf8Text(targetAmplitudes_)) + ",\"phases_rad\":"
            + arrayLiteral(utf8Text(targetPhases_)) + "}";
    if (kind == "canonical preset")
        return "{\"kind\":\"rpe_canonical\",\"name\":\""
            + jsonEscape(selectedText(target_)) + "\"}";
    return "{\"kind\":\"constant\",\"state\":"
        + arrayLiteral(utf8Text(targetState_)) + "}";
}

std::string SolverPanel::targetSettings() const {
    auto target=json::parse(targetSettingsBase());
    const auto kind=selectedText(targetKind_);
    if(!target.isNull()&&kind.rfind("observable ",0)==0) {
        if(kind.rfind("observable from state ",0)==0 || baseTargetKind(kind)=="canonical preset")
            target=json::Value::Object{{"space","observable"},{"kind","from_state"},{"state_target",target}};
        else target.at("space")="observable";
    }
    return json::dump(target);
}

}  // namespace dsl::compute
