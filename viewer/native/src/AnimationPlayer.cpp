#include "AnimationPlayer.hpp"

#include <sstream>

namespace dsl {
namespace {
std::string escapeHtml(const std::string& text) {
    std::string out;
    for (const char c : text) {
        if (c == '&') out += "&amp;";
        else if (c == '<') out += "&lt;";
        else if (c == '>') out += "&gt;";
        else if (c == '"') out += "&quot;";
        else out += c;
    }
    return out;
}
std::string jsString(const std::string& text) {
    std::string out = "\"";
    for (const char c : text) {
        if (c == '"' || c == '\\') out += '\\';
        if (c == '\n') { out += "\\n"; continue; }
        if (c == '<') { out += "\\u003c"; continue; }
        out += c;
    }
    return out + "\"";
}
}  // namespace

std::string base64(const std::string& bytes) {
    static const char* table = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve((bytes.size() + 2) / 3 * 4);
    std::size_t i = 0;
    for (; i + 2 < bytes.size(); i += 3) {
        const unsigned v = (unsigned(static_cast<unsigned char>(bytes[i])) << 16)
                         | (unsigned(static_cast<unsigned char>(bytes[i + 1])) << 8) | unsigned(static_cast<unsigned char>(bytes[i + 2]));
        out += table[(v >> 18) & 63]; out += table[(v >> 12) & 63]; out += table[(v >> 6) & 63]; out += table[v & 63];
    }
    if (i < bytes.size()) {
        unsigned v = unsigned(static_cast<unsigned char>(bytes[i])) << 16;
        if (i + 1 < bytes.size()) v |= unsigned(static_cast<unsigned char>(bytes[i + 1])) << 8;
        out += table[(v >> 18) & 63]; out += table[(v >> 12) & 63];
        out += i + 1 < bytes.size() ? table[(v >> 6) & 63] : '=';
        out += '=';
    }
    return out;
}

std::string animationPlayerHtml(const std::string& title, const std::vector<AnimationPlayerFrame>& frames,
                                double frameSeconds, bool loop, bool autoplay) {
    std::ostringstream h;
    h << "<!doctype html>\n<html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\">"
      << "<title>" << escapeHtml(title) << "</title>\n<style>"
      << ":root{--bg:#f8fafc;--fg:#0f172a;--muted:#475569;--line:#cbd5e1;--btn:#ffffff}"
      << "@media (prefers-color-scheme:dark){:root{--bg:#0f172a;--fg:#e2e8f0;--muted:#94a3b8;--line:#334155;--btn:#1e293b}}"
      << "body{margin:0;background:var(--bg);color:var(--fg);font:14px system-ui,sans-serif}"
      << ".player{max-width:1400px;margin:0 auto;padding:12px}"
      << "h1{font-size:16px;font-weight:600;margin:0 0 8px}"
      << "img{width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:6px;background:#fff}"
      << ".bar{display:flex;flex-wrap:wrap;gap:6px;align-items:center;margin-top:8px}"
      << "button,select{background:var(--btn);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:4px 10px;font:inherit;cursor:pointer}"
      << "input[type=range]{flex:1;min-width:120px}.label{color:var(--muted);min-width:120px;text-align:right}"
      << "</style></head><body><div class=\"player\">\n<h1>" << escapeHtml(title) << "</h1>\n"
      << "<img id=\"frame\" alt=\"animation frame\">\n<div class=\"bar\">"
      << "<button id=\"first\" title=\"First frame\">|&lt;</button><button id=\"back\" title=\"Previous frame\">&lt;</button>"
      << "<button id=\"play\" title=\"Play or pause\">Play</button><button id=\"forward\" title=\"Next frame\">&gt;</button>"
      << "<button id=\"stop\" title=\"Stop and rewind\">Stop</button>"
      << "<input id=\"scrub\" type=\"range\" min=\"0\" max=\"" << (frames.empty() ? 0 : frames.size() - 1) << "\" value=\"0\">"
      << "<select id=\"speed\"><option value=\"0.5\">0.5x</option><option value=\"1\" selected>1x</option>"
      << "<option value=\"2\">2x</option><option value=\"4\">4x</option></select>"
      << "<span class=\"label\" id=\"label\"></span></div></div>\n<script>\nconst frames=[";
    for (std::size_t k = 0; k < frames.size(); ++k)
        h << (k ? "," : "") << "\"data:image/png;base64," << base64(frames[k].png) << "\"";
    h << "];\nconst labels=[";
    for (std::size_t k = 0; k < frames.size(); ++k) h << (k ? "," : "") << jsString(frames[k].label);
    h << "];\nconst seconds=" << frameSeconds << ",loop=" << (loop ? "true" : "false") << ";\n"
      << "let k=0,timer=null;const img=document.getElementById('frame'),scrub=document.getElementById('scrub'),"
      << "label=document.getElementById('label'),play=document.getElementById('play'),speed=document.getElementById('speed');\n"
      << "function show(i){k=(i+frames.length)%frames.length;img.src=frames[k];scrub.value=k;"
      << "label.textContent=(labels[k]||'')+'  '+(k+1)+' / '+frames.length;}\n"
      << "function pause(){if(timer){clearInterval(timer);timer=null;}play.textContent='Play';}\n"
      << "function start(){pause();timer=setInterval(()=>{if(k+1>=frames.length&&!loop){pause();return;}show(k+1);},"
      << "seconds*1000/Number(speed.value));play.textContent='Pause';}\n"
      << "document.getElementById('first').onclick=()=>{pause();show(0);};"
      << "document.getElementById('back').onclick=()=>{pause();show(k-1);};"
      << "document.getElementById('forward').onclick=()=>{pause();show(k+1);};"
      << "document.getElementById('stop').onclick=()=>{pause();show(0);};"
      << "play.onclick=()=>{timer?pause():start();};scrub.oninput=()=>{pause();show(Number(scrub.value));};"
      << "speed.onchange=()=>{if(timer)start();};\n"
      << "frames.forEach(f=>{const i=new Image();i.src=f;});show(0);" << (autoplay ? "start();" : "") << "\n</script></body></html>\n";
    return h.str();
}

}  // namespace dsl
