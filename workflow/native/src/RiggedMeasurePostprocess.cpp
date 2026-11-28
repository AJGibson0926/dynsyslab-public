// Post-processing method "rigged_measure": Rigged DMD of observables on an existing krom.
//
// Rigged DMD (Colbrook, Drysdale & Horning, "Rigged Dynamic Mode Decomposition", arXiv:2405.00782,
// 2024) smooths the resolvent of a unitary compression of the Koopman operator to give, for one
// observable, its spectral measure on the unit circle and, at each angle, a generalized eigenfunction
// (wave packet). It applies to any krom whose spectral fit is a unitary compression, which
// measure-preserving EDMD gives (Colbrook, SIAM J. Numer. Anal. 61, 2023); the operator is checked,
// not the method's name.
//
// Settings:
//   krom      the krom to analyse (default: the stage's source); a child krom finds its spectral fit
//             through its parents.
//   measures  [{name, angles, bandwidth, order (1..8, default 1),
//               coefficients (one per dictionary feature) | reconstruction {name, output}}]
// Output: a child krom whose spectrum's continuous part holds the new measures, each with its
// density and its generalized eigenfunctions, after any the krom already had.
#include "ObjectPublishing.hpp"
#include "PostprocessInternal.hpp"
#include "KoopmanSpectralSettings.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/Text.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runRiggedMeasure(const std::filesystem::path& source, const json::Value& settings,
                                       const std::filesystem::path& root, const std::string& label) {
    namespace c = dynsyslab::core;
    const auto requested = settings.at("krom").string();
    const auto kromPath = requested.empty() ? source : std::filesystem::u8path(requested);
    if (settings.at("measures").array().empty()) throw std::invalid_argument("rigged_measure needs at least one measure");
    // Loaded from a file, the krom is published as its child; it carries its operator (compression),
    // inherited from its parents when it is a child itself.
    auto model = c::loadKoopmanModel(kromPath);
    const auto before = model.spectrum().measures().size();
    spectral_settings::riggedMeasures(model, settings.at("measures"));
    json::Value::Array names;  // the new measures, as band/observable/output
    const auto& measures = model.spectrum().measures();
    for (auto it = measures.begin() + std::ptrdiff_t(before); it != measures.end(); ++it)
        names.emplace_back(c::measureKey(it->band, it->observable, it->output));
    const json::Value meta = json::Value::Object{{"method", "rigged_measure"}, {"label", label},
        {"measures", json::Value(std::move(names))},
        {"assumptions", "Rigged DMD of a unitary compression (Colbrook, Drysdale & Horning 2024)."}};
    const auto directory = kindFolder(root, c::object_kind::krom);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__rigged-measure__" + timestamp()) + ".npz");
    c::saveKoopmanModel(path, model, withDataObject(json::dump(meta), c::object_kind::krom, {kromPath}));
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
