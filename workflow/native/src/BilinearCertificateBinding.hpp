#pragma once
// Binds a derived bilinear_certificate archive to a KLQR at deployment. The gain is used
// only when the archive is certified and was computed from this exact K-ROM and control
// map. docs/architecture/CERTIFIED_BILINEAR_CONTROL_TOOLS.md
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/LQR.hpp>

#include <string>

namespace dynsyslab::workflow::detail {

// Stated wherever a certificate is written or used.
inline constexpr const char* kBilinearCertificateScope =
    "Holds for the frozen-A bilinear model with any residual inside the stated sector, on the "
    "certified ellipsoid, with exact lifted state and the certified input hold. The sector is a "
    "data envelope; it is not a guarantee for unseen states or for the physical plant.";

struct BoundBilinearCertificate {
    core::CertifiedKoopmanFeedback feedback;
    core::json::Value evidence;
};

// Reads settings.bilinear_certificate (archive path). Throws unless it is certified and
// matches the model identity, representation, generator and map bit for bit.
[[nodiscard]] BoundBilinearCertificate bindBilinearCertificate(const core::json::Value& settings,
    const core::KoopmanModel& model, const std::string& representation);

// The certificate assumes a hold of its sample time, and an actuator cap must not cut
// into the certified input box.
void requireCertifiedActuation(const core::json::Value& settings, bool sampled, double sampleTime);

}  // namespace dynsyslab::workflow::detail
