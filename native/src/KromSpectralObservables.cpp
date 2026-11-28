#include <dynsyslab/core/KromSpectralObservables.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <memory>
#include <stdexcept>

namespace dynsyslab::core {

const FunctionExpansion& frozenDictionaryExpansion(const KoopmanModel& model) {
    if(!model.compression())throw std::invalid_argument("This krom has no compressed operator (a spectral fit's)");
    const auto features=model.compression()->dictionary.size();
    for(const auto& rep:model.representations())
        if(const auto* expansion=dynamic_cast<const FunctionExpansion*>(rep.function.get()))
            if(!expansion->coordinates().window&&expansion->dictionary().size()==features)return *expansion;
    throw std::invalid_argument("This krom has no eigenfunction expansion on its frozen dictionary");
}

std::string kromObservableIdentity(const KoopmanModel& model) {
    return model.compression() ? koopmanCompressionIdentity(*model.compression()) : koopmanModelIdentity(model);
}

Observable kromSpectralObservable(const KoopmanModel& model, ObservableDecompositionReference reference, std::string name) {
    if (reference.modelIdentity.empty() || reference.modelIdentity != kromObservableIdentity(model))
        throw std::invalid_argument("Observable reference does not match the krom's identity");
    ObservableDescriptor descriptor;
    descriptor.name = std::move(name); descriptor.label = descriptor.name;
    descriptor.complexValued = true; descriptor.decompositions = {reference};
    if (!model.compression()) {  // a krom without an operator: its named readouts
        if (reference.kind != "reconstruction")
            throw std::invalid_argument("A krom without a compressed operator offers only its named readouts");
        if (model.requiresExternalCoordinates())
            throw std::invalid_argument("This krom's readouts need its external coordinate provider");
        return Observable(descriptor, model.stateDimension(), model.reconstruction(reference.component).outputs,
                          [model, reference](const double* x) { return model.reconstruct(reference.component, x); });
    }
    const auto compression = std::make_shared<const KoopmanCompression>(*model.compression());
    descriptor.stateCoordinateSystem = compression->coordinateContract;
    const auto dimension = compression->dictionary.dimension, features = compression->dictionary.size();
    // Coefficients on the frozen dictionary (outputs x features) of the referenced functions.
    std::vector<std::complex<double>> coefficients;
    std::size_t outputs = 1;
    if (reference.kind == "reconstruction") {
        outputs = compression->readouts.at(reference.component).outputs;
    } else if (reference.kind == "eigenfunction") {
        if (reference.index >= model.liftedDimension()) throw std::invalid_argument("Eigenfunction reference index out of range");
        const auto& all = frozenDictionaryExpansion(model).coefficients();
        coefficients.assign(all.begin() + std::ptrdiff_t(reference.index * features),
                            all.begin() + std::ptrdiff_t((reference.index + 1) * features));
    } else if (reference.kind == "wave-packet") {
        // The measure by key (component/observable/output) when the reference names one, else by position.
        const auto& measures = model.spectrum().measures();
        const SpectralMeasure* measure = nullptr;
        for (const auto& m : measures)
            if (!reference.component.empty() && measureKey(m.band, m.observable, m.output) == reference.component) measure = &m;
        if (reference.component.empty() && reference.index < measures.size()) measure = &measures[reference.index];
        const auto* packets = measure ? measure->packets.get() : nullptr;
        if (!packets) throw std::invalid_argument("Wave packet reference: no measure with generalized eigenfunctions there");
        coefficients = packets->coefficients();
        outputs = packets->outputCount();
    } else {
        throw std::invalid_argument("Spectral measures are associated data, not pointwise observables");
    }
    return Observable(descriptor, dimension, outputs,
        [compression, reference, coefficients = std::move(coefficients), dimension, features, outputs](const double* point) {
            const std::vector<double> x(point, point + dimension);
            Observable::Values out;
            if (reference.kind == "reconstruction") {
                for (const double v : compression->reconstruct(reference.component, x)) out.emplace_back(v);
                return out;
            }
            const auto phi = compression->dictionary.evaluate(x);
            out.assign(outputs, {});
            for (std::size_t o = 0; o < outputs; ++o)
                for (std::size_t k = 0; k < features; ++k) out[o] += coefficients[o * features + k] * phi[k];
            return out;
        });
}

Observable existingKromObservable(const std::filesystem::path& krom, ObservableDecompositionReference reference, std::string name) {
    reference.artifact = krom.generic_u8string();
    return kromSpectralObservable(loadKoopmanModel(krom), std::move(reference), std::move(name));
}

}  // namespace dynsyslab::core
