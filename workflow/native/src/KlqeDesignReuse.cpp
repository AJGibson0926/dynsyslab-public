// KLQE design reuse (estimator key "design_archive"): the gain, covariance and observer poles of
// a published klqe stage are reloaded instead of redesigned, after verifying that the archive was
// designed for the same K-ROM (its krom_identity equals the bound model's identity) and that its
// measurement matrix equals the one bound now. A mismatch refuses the job.
#include "KlqeDesignReuse.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow {

json::Value reuseKlqeDesign(core::KoopmanLQEDesign& design, const std::string& path, const core::KoopmanModel& model) {
    const auto archive = core::readNpz(core::wide(path));
    const auto meta = json::parse(archive.metadataJson());
    const auto identity = meta.at("krom_identity").string();
    const auto bound = core::koopmanModelIdentity(model);
    if (identity.empty()) throw std::invalid_argument("KLQE design archive records no krom_identity: " + path);
    if (identity != bound) throw std::invalid_argument("KLQE design archive was designed for another K-ROM: " + path);
    for (const char* name : {"estimator_gain", "estimator_covariance", "observer_poles", "measurement_matrix"})
        if (!archive.contains(name)) throw std::invalid_argument(std::string("KLQE design archive lacks ") + name);
    const auto C = archive.at("measurement_matrix").doubles();
    if (C.size() != design.measurementMatrix.size()) throw std::invalid_argument("KLQE design archive measurement differs in shape");
    double difference = 0.0, scale = 0.0;
    for (std::size_t k = 0; k < C.size(); ++k) {
        difference = std::max(difference, std::abs(C[k] - design.measurementMatrix[k]));
        scale = std::max(scale, std::abs(design.measurementMatrix[k]));
    }
    if (difference > 1e-9 * std::max(1.0, scale)) throw std::invalid_argument("KLQE design archive measurement matrix differs");
    auto gain = archive.at("estimator_gain").doubles();
    if (gain.size() != design.solution.gain.size()) throw std::invalid_argument("KLQE design archive gain differs in shape");
    design.solution.gain = std::move(gain);
    design.solution.covariance = archive.at("estimator_covariance").doubles();
    design.solution.observerPoles = archive.at("observer_poles").complexDoubles();
    return json::Value::Object{{"reused_design", path}, {"krom_identity", identity}};
}

}  // namespace dynsyslab::workflow
