// Post-processing method "animation": strings views of existing archives into a Viewer animation
// (dynsyslab.animation/v1). The archive holds only the manifest; Viewer builds frame k from page k
// of every panel (viewer/native/src/AnimationPayload.cpp).
//
//   {"name": "animation", "settings": {
//      "panels": [{"view": "domain_field", "label": "angle"}, {"view": "time"}],  one per source, in order
//      "frame_seconds": 0.4, "loop": true, "autoplay": true, "columns": 2, "title": "..."}}
// The stage's sources are the panel archives. A missing panel entry uses view "domain_field" for
// grid-field archives and "time" otherwise.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::filesystem::path runAnimation(const std::vector<std::filesystem::path>& sources, const J& settings,
                                   const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("animation needs panel archives as its sources");
    const auto& panels = settings.at("panels");
    const auto count = panels.isNull() ? std::size_t{0} : panels.array().size();
    if (count > sources.size()) throw std::invalid_argument("animation lists more panels than sources");
    J::Array out;
    for (std::size_t i = 0; i < sources.size(); ++i) {
        if (!std::filesystem::exists(sources[i])) throw std::invalid_argument("animation panel archive is missing: " + sources[i].u8string());
        const auto meta = json::parse(c::readNpzIndex(sources[i]).metadataJson);
        const auto view = i < count ? panels.array()[i].at("view").string() : std::string();
        const auto fallback = meta.at("format").string() == "dynsyslab.grid-field/v1" ? "domain_field" : "time";
        out.push_back(J::Object{{"source", std::filesystem::absolute(sources[i]).generic_u8string()},
                                {"view", view.empty() ? fallback : view},
                                {"label", i < count ? panels.array()[i].at("label").string() : std::string()}});
    }
    J::Object meta{{"format", "dynsyslab.animation/v1"}, {"method", "animation"}, {"label", label},
                   {"title", settings.at("title").string(label)}, {"panels", out},
                   {"frame_seconds", settings.at("frame_seconds").number(0.4)},
                   {"loop", settings.at("loop").boolean(true)}, {"autoplay", settings.at("autoplay").boolean(true)},
                   {"source_archives", [&] { J::Array a; for (const auto& s : sources) a.emplace_back(std::filesystem::absolute(s).generic_u8string()); return a; }()}};
    if (!settings.at("columns").isNull()) meta["columns"] = settings.at("columns").number();
    c::NpzRecord archive;
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    const auto directory = objectFolder(root, "animation");
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    c::writeNpzAtomicNew(path, archive);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
