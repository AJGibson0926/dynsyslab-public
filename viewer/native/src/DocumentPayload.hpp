#pragma once
// Viewer payloads for JSON result documents (no archive): the presentation of results that native
// queries already computed. The Viewer only arranges the supplied numbers; it fits nothing.
//
//   periodic-analysis spectral_compare result ("method": "periodic-spectral-comparison"):
//       one page per observable, family and packet centre: the spectral measure (symlog) and the
//       RMS packet amplitude across orbits; view "frequency" (default) or "parameter" as abscissa.
//   dynsyslab.phase-cover-geometry/v1 (analytic query phase_cover): one panel per tower level, the
//       sampled orbit coloured by mapped base phase with its compatible phase points; options
//       "base_phase_turns" and "digits" in the document's "display" object.
//   dynsyslab.research-figure/v1 {"figure": kind, ...}: the research figures of
//       ResearchFigurePayload.cpp (finite-cycle spectra, cocycles, coordinate reconstruction,
//       clock changes).
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

[[nodiscard]] bool isDocumentPayload(const dynsyslab::core::json::Value& document);
[[nodiscard]] bool isDocumentPayloadFile(const std::filesystem::path& path);
[[nodiscard]] Payload buildDocumentPayload(const std::filesystem::path& source, const std::string& view);
[[nodiscard]] Payload buildDocumentPayload(const dynsyslab::core::json::Value& document, const std::string& view,
                                           const std::string& source = {});
// dynsyslab.research-figure/v1 documents (ResearchFigurePayload.cpp).
[[nodiscard]] Payload buildResearchFigurePayload(const dynsyslab::core::json::Value& document);

}  // namespace dsl
