#include "BilinearCertificateBinding.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanBilinearBridge.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/ControlFactory.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {
namespace c = dynsyslab::core;
namespace json = c::json;

std::vector<double> required(const c::NpzRecord& archive, const char* name, std::size_t size) {
    if (!archive.contains(name)) throw std::invalid_argument(std::string("bilinear certificate lacks ") + name);
    auto values = archive.at(name).doubles();
    if (values.size() != size) throw std::invalid_argument(std::string("bilinear certificate array has the wrong size: ") + name);
    return values;
}
}  // namespace

BoundBilinearCertificate bindBilinearCertificate(const json::Value& settings, const c::KoopmanModel& model,
                                                 const std::string& representation) {
    const auto path = settings.at("bilinear_certificate").string();
    if (path.empty()) throw std::invalid_argument("bilinear_certificate must be an archive path");
    const auto archive = c::readNpz(c::wide(path));
    const auto meta = json::parse(archive.metadataJson());
    if (meta.at("method").string() != "koopman-bilinear-certificate")
        throw std::invalid_argument("bilinear_certificate is not a bilinear certificate archive");
    const auto& record = meta.at("bilinear_certificate");
    const auto status = record.at("status").string();
    if (status != "certified" || !record.at("certified").boolean(false))
        throw std::invalid_argument("bilinear_certificate has status '" + status + "' and carries no usable gain");
    if (record.at("model_identity").string() != c::koopmanModelIdentity(model))
        throw std::invalid_argument("bilinear_certificate belongs to a different K-ROM");
    if (record.at("control_representation").string() != representation)
        throw std::invalid_argument("bilinear_certificate was computed for another control_representation");

    const auto split = c::continuousBilinearSplit(model, representation);
    const auto n = split.stateDimension, m = split.inputDimension;
    if (required(archive, "continuous_generator", n * n) != split.generator
        || required(archive, "continuous_offsets", n * m) != split.offsets)
        throw std::invalid_argument("bilinear_certificate generator or input offset differs from this K-ROM");
    std::vector<double> operators;
    for (const auto& matrix : split.operators) operators.insert(operators.end(), matrix.begin(), matrix.end());
    if (required(archive, "continuous_operators", m * n * n) != operators)
        throw std::invalid_argument("bilinear_certificate bilinear operators differ from this K-ROM");

    BoundBilinearCertificate bound;
    bound.feedback.gain = required(archive, "certificate_K", m * n);
    bound.feedback.targetLifted = required(archive, "target_lifted", n);
    bound.feedback.targetInput = required(archive, "target_input", m);
    bound.feedback.sampleTime = record.at("sample_time").number(0.0);
    json::Value::Object evidence{
        {"archive", json::Value(path)}, {"status", json::Value(status)},
        {"control_law", json::Value("u = u* + K (z - z*), z = split(psi)")},
        {"operators", record.at("operators")}, {"sector", record.at("sector")},
        {"input_bounds", record.at("input_bounds")}, {"state_radius", record.at("state_radius")},
        {"contraction", record.at("contraction")}, {"sample_time", record.at("sample_time")},
        {"scope", json::Value(kBilinearCertificateScope)}};
    bound.evidence = json::Value(std::move(evidence));
    return bound;
}

void requireCertifiedActuation(const json::Value& settings, bool sampled, double sampleTime) {
    const auto archive = c::readNpz(c::wide(settings.at("bilinear_certificate").string()));
    const auto meta = json::parse(archive.metadataJson());
    const auto& record = meta.at("bilinear_certificate");
    const double dt = record.at("sample_time").number(0.0);
    if (!sampled || !(std::abs(sampleTime - dt) <= 1e-9 * dt))
        throw std::invalid_argument("bilinear_certificate needs sampled control at its certified sample time");
    if (const auto cap = controlCap(settings.at("control_limit"))) {
        const auto targetInput = archive.at("target_input").doubles();
        const auto& bounds = record.at("input_bounds").array();
        if (targetInput.size() != 1 || bounds.size() != 1)
            throw std::invalid_argument("control_limit with a certificate supports one input");
        const double bound = bounds.front().number();
        if (cap->first > targetInput[0] - bound || cap->second < targetInput[0] + bound)
            throw std::invalid_argument("control_limit cuts into the certified input box, which would void the certificate");
    }
}

}  // namespace dynsyslab::workflow::detail
