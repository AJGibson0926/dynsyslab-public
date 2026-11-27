// K-ROM archives (dynsyslab.krom/v1, KoopmanArchive.hpp): one data object holding a KoopmanModel, or
// the parts a step added to a parent krom it references. Loading follows parents and lays each
// archive's parts over its parent's chosen modes.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/core/SpectrumArchive.hpp>
#include <dynsyslab/core/System.hpp>
#include "KoopmanArchiveMetadata.hpp"
#include "KoopmanBinding.hpp"
#include "KoopmanInputMaps.hpp"
#include "KromArchiveParts.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

json::Value metadataOf(const NpzRecord& archive) {
    return archive.metadataJson().empty() ? json::Value{} : json::parse(archive.metadataJson());
}

std::vector<std::size_t> allModes(std::size_t count) {
    std::vector<std::size_t> modes(count);
    for (std::size_t j = 0; j < count; ++j) modes[j] = j;
    return modes;
}

// One krom archive on `system` (resolved once, for the whole chain), its parents resolved.
KoopmanModel readNode(const NpzRecord& archive, const std::string& label, const SystemConfiguration& system, int depth) {
    const auto metadata = metadataOf(archive);
    const auto& block = metadata.at("krom");
    if (metadata.at("schema").string() != kKromSchema || block.isNull())
        throw std::invalid_argument(label + " is not a krom (" + std::string(kKromSchema) + "); regenerate it");
    if (depth > 64) throw std::invalid_argument("krom parents form a cycle at " + label);
    const auto origin = dataOriginFromJson(block.at("origin"));
    const auto& parentBlock = block.at("parent");
    std::optional<KoopmanModel> model;
    if (parentBlock.isNull()) {
        auto spectrum = readSpectrumArrays(archive, origin);
        auto samples = detail::readKromSamples(archive, spectrum.size());
        auto functions = detail::readKromRepresentations(block.at("representations"), archive, samples, spectrum.size());
        model.emplace(std::move(spectrum), std::move(samples), origin, std::move(functions),
                      std::size_t(metadata.at("input_dimension").number(0)));
        if (!system.identifier.empty()) model->bindSystem(system);
    } else {
        const auto path = std::filesystem::u8path(parentBlock.at("path").string());
        if (!std::filesystem::is_regular_file(extendedLengthPath(path)))
            throw std::invalid_argument("the parent krom of " + label + " is missing: " + path.u8string());
        if (sha256File(path) != parentBlock.at("sha256").string())
            throw std::invalid_argument("the parent krom of " + label + " has changed: " + path.u8string());
        std::vector<std::size_t> modes;
        for (const auto& mode : parentBlock.at("modes").array()) modes.push_back(std::size_t(mode.number()));
        auto base = readNode(readNpz(path), path.u8string(), system, depth + 1);
        model.emplace(modes.empty() ? std::move(base) : base.extract(modes));
        for (auto& function : detail::readKromRepresentations(block.at("representations"), archive, model->samples(), model->liftedDimension())) {
            const auto name = representationName(*function);
            const auto names = model->representationNames();
            if (std::find(names.begin(), names.end(), name) != names.end()) model->replaceRepresentation(std::move(function));
            else model->addRepresentation(std::move(function), false);
        }
        model->setOrigin(origin);
    }
    const auto first = block.at("default_representation").string();
    if (!first.empty()) model->useRepresentation(first);
    for (const auto& ref : koopmanDatasetReferencesFromMetadata(metadata)) model->addDatasetReference(ref);
    detail::readKromReadouts(archive, *model);
    detail::readNamedInputMaps(archive, *model);
    // The spectrum's evidence (pair residuals, pseudospectrum, continuous part) and the fit's
    // diagnostics: a root's are read with its spectrum; a child's own replace what it inherited.
    if (!parentBlock.isNull() && hasSpectralEvidence(archive)) {
        auto evidence = model->spectrum();
        readSpectralEvidence(archive, evidence);
        model->setSpectralEvidence(evidence);
    }
    if (archive.contains("spectral_diagnostics_json"))
        model->setSpectralDiagnostics(spectralDiagnosticsFromJson(json::parse(archive.at("spectral_diagnostics_json").string())));
    return std::move(*model);
}

}  // namespace

bool isKromMetadata(const std::string& metadataJson) {
    if (metadataJson.empty()) return false;
    const auto metadata = json::parse(metadataJson);
    return metadata.at("schema").string() == kKromSchema && !metadata.at("krom").isNull();
}

bool isKrom(const NpzRecord& archive) { return archive.contains("metadata_json") && isKromMetadata(archive.metadataJson()); }

KoopmanModel loadKoopmanModel(const std::filesystem::path& path, const SystemConfiguration& system,
                              const std::vector<std::size_t>& selectedModes, const ModelBindingOptions& binding) {
    const auto archive = readNpz(path);
    const auto resolved = detail::resolveArchiveBinding(archive, system, binding);
    auto model = readNode(archive, path.u8string(), resolved.system, 0);
    model.setBindingEvidenceJson(resolved.evidence);
    model.setParent({std::filesystem::absolute(path).lexically_normal(), sha256File(path), allModes(model.liftedDimension())});
    // The binding's choices are this model's own: a krom written from it records them.
    if (!binding.interpolationDomain.empty()) {
        model.useRepresentation("interpolation");
        const auto own = domainName(model.interpolation().coordinates());
        if (own != binding.interpolationDomain)
            throw std::invalid_argument("these eigenfunctions are interpolated on " + own + ", not " + binding.interpolationDomain
                                        + "; the domain belongs to the eigenfunctions and is fixed when they are made");
    }
    if (!binding.interpolationMethod.isNull()) {
        std::size_t neighbors = model.interpolation().data().neighbors;
        auto method = interpolationMethodFromJson(binding.interpolationMethod, &neighbors);
        model.setInterpolationMethod(std::move(method), neighbors);
    }
    if (!binding.representation.empty()) model.useRepresentation(binding.representation);
    return selectedModes.empty() ? std::move(model) : model.extract(selectedModes);
}

KoopmanModel loadKoopmanModel(const std::filesystem::path& path, const std::vector<std::size_t>& selectedModes) {
    return loadKoopmanModel(path, detail::systemFromMetadata(readNpz(path)), selectedModes);
}

KoopmanModel loadKoopmanModel(const std::filesystem::path& path, const ModelBindingOptions& binding) {
    return loadKoopmanModel(path, detail::systemFromMetadata(readNpz(path)), {}, binding);
}

KoopmanModel loadKoopmanModelOnDomain(const std::filesystem::path& path, const std::string& domain, const json::Value& method) {
    ModelBindingOptions binding;
    binding.interpolationDomain = domain;
    binding.interpolationMethod = method;
    return loadKoopmanModel(path, detail::systemFromMetadata(readNpz(path)), {}, binding);
}

KoopmanModel readKoopmanModel(const NpzRecord& archive) {
    return readNode(archive, "this archive", detail::systemFromMetadata(archive), 0);
}

KoopmanModel spectralKromModel(const NpzRecord& archive, const SystemConfiguration& system) {
    // A spectral K-ROM in its own frozen dictionary: psi_j(q) = sum_k E_kj phi_k(q). The domain follows
    // its Hankel window: one delay of a solve's full state is the physical state; one delay of named
    // coordinate columns is those catalog observables; longer windows are data coordinates.
    const auto metadata = metadataOf(archive);
    const auto krom = readSpectralKromArchive(archive);
    const std::size_t features = krom.dictionary.size(), modes = krom.multipliers.size();
    if (!features || !modes || krom.eigenfunctions.size() != features * modes)
        throw std::invalid_argument("spectral K-ROM eigenfunction coefficients have the wrong shape");
    std::vector<std::complex<double>> coefficients(modes * features);
    for (std::size_t j = 0; j < modes; ++j)
        for (std::size_t k = 0; k < features; ++k) coefficients[j * features + k] = krom.eigenfunctions[k * modes + j];
    auto spectrum = Spectrum::fromDiscreteMultipliers(krom.multipliers, krom.sampleInterval);
    FunctionDomain domain; domain.dimension = krom.dictionary.dimension;
    json::Value contract;
    try { contract = json::parse(krom.coordinateContract); } catch (const std::exception&) {}
    const auto& channels = contract.at("channels");
    bool oneDelay = contract.at("delays").number(0) == 1 && !channels.isNull() && channels.array().size() == domain.dimension;
    std::vector<std::string> observables;
    for (std::size_t k = 0; oneDelay && k < channels.array().size(); ++k) {
        const auto& channel = channels.array()[k];
        oneDelay = channel.at("lag_samples").number(0) == 0 && channel.at("measurement_column").number(-1) == double(k);
        observables.push_back(channel.at("observable").string());
    }
    const bool fromSolve = metadata.at("coordinate_transform_kind").string().empty() && metadata.at("value_space").string().empty();
    const auto physical = createSystem(system)->descriptor().stateDimension;
    if (oneDelay && fromSolve && physical == domain.dimension) domain = FunctionDomain::physicalState(domain.dimension);
    else if (oneDelay && !fromSolve && std::none_of(observables.begin(), observables.end(), [](const std::string& n) { return n.empty(); }))
        domain = FunctionDomain::fromCatalog(observables, domain.dimension);
    FunctionExpansion functions(std::make_shared<const ExpansionDictionary>(ExpansionDictionary::fromSpectral(krom.dictionary)),
                                std::move(coefficients), std::move(domain));
    return KoopmanModel(std::move(functions), std::move(spectrum), system);
}

NpzRecord makeKoopmanArchive(const KoopmanModel& model, const std::string& metadataJson,
                             std::map<std::string, NpyArray> extraArrays,
                             const std::vector<ChannelDescriptor>& inputChannels) {
    const auto& parent = model.parent();
    const bool child = parent.has_value();
    NpzRecord archive;
    archive.arrays = std::move(extraArrays);
    if (!child) {
        writeSpectrumArrays(archive, model.spectrum());
        detail::writeKromSamples(archive, model.samples());
    }
    json::Value::Array order;
    for (const auto& name : model.representationNames()) order.emplace_back(name);
    json::Value::Object block{{"origin", dataOriginJson(model.origin())},
        {"input_channels", channelDescriptorsJson(inputChannels)},
        {"representations", detail::writeKromRepresentations(archive, model, child)},
        {"order", order}, {"default_representation", model.representation().name}, {"parent", json::Value()}};
    if (child) {
        json::Value::Array modes;
        for (const auto mode : parent->modes) modes.emplace_back(double(mode));
        block["parent"] = json::Value::Object{
            {"path", std::filesystem::absolute(parent->path).lexically_normal().generic_u8string()},
            {"sha256", parent->sha256}, {"modes", modes}};
    }
    detail::writeKromReadouts(archive, model, child);
    detail::writeNamedInputMaps(archive, model, child);
    if (child && model.ownsSpectralEvidence()) writeSpectralEvidence(archive, model.spectrum());
    if (!child || model.ownsSpectralEvidence())
        archive.arrays["spectral_diagnostics_json"] = NpyArray::fromString(json::dump(spectralDiagnosticsJson(model.spectralDiagnostics())));

    json::Value metadata = metadataJson.empty() ? json::Value(json::Value::Object{}) : json::parse(metadataJson);
    if (!std::holds_alternative<json::Value::Object>(metadata.storage())) throw std::invalid_argument("KROM metadata must be a JSON object");
    auto references = model.datasetReferences();
    for (const auto& ref : koopmanDatasetReferencesFromMetadata(metadata)) mergeKoopmanDatasetReference(references, ref);
    metadata.object()["dataset_references"] = koopmanDatasetReferencesJson(references);
    if (model.hasPhysicalSystem()) {
        metadata.object()["system_configuration"] = detail::koopmanSystemJson(model.systemConfiguration());
        metadata.object()["system"] = model.systemConfiguration().identifier;
    } else {
        metadata.object()["input_dimension"] = static_cast<double>(model.inputDimension());  // no system to give it
    }
    metadata.object()["schema"] = kKromSchema;
    metadata.object()["krom"] = json::Value(std::move(block));
    metadata.object()["eigenfunction_representation"] = model.representation().name;
    metadata.object()["binding_evidence"] = json::parse(model.bindingEvidenceJson());
    metadata.object()["krom_mode_count"] = static_cast<double>(model.liftedDimension());
    const bool analytic = !model.scatteredEigenfunctions() && !model.representation().function->coordinates().window
                          && model.expansion().dictionary().exactGradient();
    metadata.object()["gradient_backend"] = analytic ? "analytic" : "finite-difference";
    // Every krom declares its kind; a writer's envelope keeps its provenance.
    const auto envelope = dataObjectFromMetadata(json::dump(metadata));
    DataObjectInfo info = envelope ? *envelope : DataObjectInfo{object_kind::krom, model.origin(), {}};
    info.kind = object_kind::krom;
    if (!envelope && child)
        info.producedBy.inputs.push_back({std::filesystem::absolute(parent->path).lexically_normal().generic_u8string(), parent->sha256, "parent"});
    metadata.object()["data_object"] = dataObjectJson(info);
    archive.arrays["metadata_json"] = NpyArray::fromString(json::dump(metadata));
    return archive;
}

void saveKoopmanModel(const std::filesystem::path& path, const KoopmanModel& model, const std::string& metadataJson,
                      std::map<std::string, NpyArray> extraArrays) {
    writeNpzAtomic(path, makeKoopmanArchive(model, metadataJson, std::move(extraArrays)));
}

}  // namespace dynsyslab::core
