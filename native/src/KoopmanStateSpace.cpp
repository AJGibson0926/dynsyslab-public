#include <dynsyslab/core/KoopmanStateSpace.hpp>
#include <dynsyslab/core/Observability.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cmath>
#include <stdexcept>
#include <variant>

namespace dynsyslab::core {
namespace {

std::vector<double> splitTransition(const KoopmanModel& model, double dt) {
    const auto mu = model.discreteEigenvalues(dt);
    const std::size_t n = model.splitDimension();
    std::vector<double> A(n * n, 0.0);
    std::size_t j = 0;
    for (std::size_t k = 0; k < mu.size(); ++k) {
        A[j * n + j] = mu[k].real();
        if (std::abs(model.continuousEigenvalues()[k].imag()) > 1e-12) {
            A[j * n + j + 1] = -mu[k].imag(); A[(j + 1) * n + j] = mu[k].imag();
            A[(j + 1) * n + j + 1] = mu[k].real(); j += 2;
        } else ++j;
    }
    return A;
}

std::vector<double> inputMatrix(const KoopmanModel& model, const KoopmanStateSpaceOptions& o) {
    const auto stored = model.controlMaps().find(o.controlMap);  // "exact" is built in, not stored
    if (stored != model.controlMaps().end())
        if (const auto* constant = std::get_if<ConstantControlMap>(&stored->second))
            return model.realControlMap(constant->matrix);
    if (!o.operatingState.empty() && !o.operatingLifted.empty())
        throw std::invalid_argument("give either an operating state or operating lifted values, not both");
    if (!o.operatingState.empty())
        return model.controlMapSplit(o.operatingState.data(), o.controlMap);
    if (!o.operatingLifted.empty()) {
        if (o.operatingLifted.size() != model.splitDimension())
            throw std::invalid_argument("operating lifted values must match the KROM split dimension");
        return model.realControlMap(model.controlMapFromLifted(model.valuesFromSplit(o.operatingLifted), o.controlMap));
    }
    throw std::invalid_argument("control map '" + o.controlMap +
        "' is state dependent; an LTI realization needs an explicit operating point");
}

}  // namespace

StateSpaceModel koopmanStateSpace(const KoopmanModel& model, const KoopmanStateSpaceOptions& o) {
    const std::size_t n = model.splitDimension();
    std::vector<double> C;
    std::size_t outputs = 0;
    for (const auto& name : o.measurements) {
        const auto rows = reconstructionSplitMatrix(model, model.reconstruction(name));
        C.insert(C.end(), rows.begin(), rows.end());
        outputs += rows.size() / n;
    }
    const bool discrete = o.domain == TimeDomain::discrete;
    if (discrete && !o.controlMap.empty()) {
        const auto found = model.discreteInputMaps().find(o.controlMap);
        if (found != model.discreteInputMaps().end()) {
            const auto& map = found->second;
            if (o.sampleTime != 0.0 && std::abs(o.sampleTime - map.sampleInterval) > 1e-12 * map.sampleInterval)
                throw std::invalid_argument("sampleTime differs from the stored discrete input map");
            return makeStateSpace(splitTransition(model, map.sampleInterval), model.realControlMap(map.matrix), C, {},
                                  n, map.inputDimension, outputs, TimeDomain::discrete, map.sampleInterval);
        }
    }
    std::vector<double> B;
    std::size_t inputs = 0;
    if (!o.controlMap.empty()) { B = inputMatrix(model, o); inputs = B.size() / n; }
    auto continuous = makeStateSpace(model.generatorSplit(), B, C, {}, n, inputs, outputs);
    if (!discrete) return continuous;
    if (!(o.sampleTime > 0.0)) throw std::invalid_argument("discrete KROM realization needs a positive sampleTime");
    return discretize(continuous, o.sampleTime);
}

}  // namespace dynsyslab::core
