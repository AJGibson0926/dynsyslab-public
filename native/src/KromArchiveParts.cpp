#include "KromArchiveParts.hpp"
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/HankelRepresentations.hpp>

#include <cctype>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {

std::string safeKey(std::string value) {
    for (char& item : value)
        if (!std::isalnum(static_cast<unsigned char>(item)) && item != '_' && item != '-') item = '_';
    return value;
}

std::vector<std::string> splitNames(const std::string& source) {
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin < source.size()) {
        const auto end = source.find('\n', begin);
        result.push_back(source.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
}

}  // namespace

void writeKromSamples(NpzRecord& archive, const EigenfunctionSamples& samples) {
    const std::size_t n = samples.samples(), m = samples.modes;
    if (!n) return;  // e.g. a fitted expansion: no training samples
    // A method's own array of the same name and shape is the more precise record (e.g. the true
    // sample times beside an interpolation's sample order) and is kept; any other clash is an error.
    const auto put = [&](const char* key, NpyArray value) {
        if (archive.contains(key)) {
            if (archive.at(key).shape != value.shape)
                throw std::logic_error(std::string("a K-ROM's samples and a method's array share the name ") + key);
            return;
        }
        archive.arrays[key] = std::move(value);
    };
    put("eigenfunction_values", NpyArray::fromComplex({n, m}, samples.values));
    put("training_time", NpyArray::fromDoubles({n}, samples.time));
    std::vector<double> indices(samples.sampleIndices.begin(), samples.sampleIndices.end());
    put("training_sample_indices", NpyArray::fromDoubles({n}, indices));
}

EigenfunctionSamples readKromSamples(const NpzRecord& archive, std::size_t modes) {
    EigenfunctionSamples samples;
    samples.modes = modes;
    // All three or none: a method's own array named like one of them is its property, not samples.
    if (!archive.contains("eigenfunction_values") || !archive.contains("training_time")
        || !archive.contains("training_sample_indices")) return samples;
    samples.time = archive.at("training_time").doubles();
    for (const double index : archive.at("training_sample_indices").doubles()) samples.sampleIndices.push_back(std::size_t(index));
    samples.values = archive.at("eigenfunction_values").complexDoubles();
    return samples;
}

json::Value::Array writeKromRepresentations(NpzRecord& archive, const KoopmanModel& model, bool ownOnly) {
    const std::size_t m = model.liftedDimension();
    json::Value::Array representations;
    for (const auto& rep : model.representations()) {
        if (ownOnly && !model.ownsRepresentation(rep.name)) continue;
        const auto& function = *rep.function;
        const auto* expansion = dynamic_cast<const FunctionExpansion*>(&function);
        if (const auto& window = function.coordinates().window) {
            // A window map: its layout and the matrix T (row j is eigenfunction j's own map).
            if (!expansion) throw std::invalid_argument("a window-domain representation must be a window map");
            archive.arrays["window_map"] = NpyArray::fromComplex({m, window->rows()}, windowMapMatrix(*expansion));
            auto description = windowMapDescription(*window);
            description.object()["name"] = rep.name;
            representations.emplace_back(std::move(description));
            continue;
        }
        json::Value::Object description{{"kind", function.kind()}, {"name", rep.name},
            {"domain_name", domainName(function.coordinates())}, {"domain", functionDomainJson(function.coordinates())}};
        if (expansion) {
            // One dictionary; row j of the coefficients is eigenfunction j's own expansion.
            const std::string array = archive.contains("expansion_coefficients")
                ? "expansion_coefficients__" + safeKey(rep.name) : std::string("expansion_coefficients");
            description["dictionary"] = expansion->dictionary().describe();
            description["coefficients"] = array;
            archive.arrays[array] = NpyArray::fromComplex({m, expansion->dictionary().size()}, expansion->coefficients());
        } else if (const auto* interpolation = dynamic_cast<const FunctionInterpolation*>(&function)) {
            // The base points belong to the interpolation; its values are the samples' (eigenfunction_values).
            if (interpolation->data().samples != model.samples().samples()
                || interpolation->data().trainingValues != model.samples().values)
                throw std::invalid_argument("an archived interpolation must be of the K-ROM's own samples");
            description["method"] = interpolationMethodJson(interpolation->data().method, interpolation->data().neighbors);
            description["points"] = "interpolation_points";
            archive.arrays["interpolation_points"] = NpyArray::fromDoubles(
                {interpolation->data().samples, function.coordinateDimension()}, interpolation->points());
        } else {
            throw std::invalid_argument("the " + rep.name + " representation cannot be archived");
        }
        representations.emplace_back(json::Value(std::move(description)));
    }
    return representations;
}

std::vector<std::shared_ptr<const FunctionRepresentation>> readKromRepresentations(
    const json::Value& descriptions, const NpzRecord& archive, const EigenfunctionSamples& samples, std::size_t modes) {
    std::vector<std::shared_ptr<const FunctionRepresentation>> out;
    for (const auto& description : descriptions.array()) {
        const auto kind = description.at("kind").string();
        if (kind == "interpolation") {
            const auto& points = archive.at("interpolation_points");
            InterpolationData data;
            data.samples = samples.samples();
            data.stateDimension = points.shape.at(1);
            data.trainingStates = points.doubles();
            data.trainingValues = samples.values;
            data.method = interpolationMethodFromJson(description.at("method"), &data.neighbors);
            data.provenance = "eigenfunction interpolation";
            out.push_back(std::make_shared<const FunctionInterpolation>(std::move(data), functionDomainFromJson(description.at("domain"))));
        } else if (kind == "window_map" || description.at("family").string() == "linear-window") {
            out.push_back(windowMapFromMap(hankelWindowLayoutFromJson(description), archive.at("window_map").complexDoubles(), modes));
        } else if (kind == "expansion") {
            auto dictionary = std::make_shared<const ExpansionDictionary>(ExpansionDictionary::fromJson(description.at("dictionary")));
            out.push_back(std::make_shared<const FunctionExpansion>(std::move(dictionary),
                archive.at(description.at("coefficients").string("expansion_coefficients")).complexDoubles(),
                functionDomainFromJson(description.at("domain"))));
        } else {
            throw std::invalid_argument("unknown eigenfunction representation " + kind);
        }
    }
    return out;
}

void writeKromReadouts(NpzRecord& archive, const KoopmanModel& model, bool ownOnly) {
    std::string names;
    for (const auto& [name, readout] : model.reconstructions()) {
        if (ownOnly && !model.ownsReconstruction(name)) continue;
        const auto key = safeKey(name);
        if (archive.contains("reconstruction__" + key)) throw std::invalid_argument("readout archive keys collide: " + name);
        archive.arrays["reconstruction__" + key] = NpyArray::fromComplex({readout.outputs, readout.modes}, readout.matrix);
        if (!readout.offset.empty())
            archive.arrays["reconstruction_offset__" + key] = NpyArray::fromDoubles({readout.offset.size()}, readout.offset);
        if (readout.residualRms != 0 || readout.residualMaximum != 0 || readout.relativeResidual != 0)
            archive.arrays["reconstruction_error__" + key] = NpyArray::fromDoubles(
                {3}, {readout.residualRms, readout.residualMaximum, readout.relativeResidual});
        names += (names.empty() ? "" : "\n") + name;
    }
    if (!names.empty()) archive.arrays["reconstruction_names"] = NpyArray::fromString(names);
}

void readKromReadouts(const NpzRecord& archive, KoopmanModel& model) {
    if (!archive.contains("reconstruction_names")) return;
    for (const auto& name : splitNames(archive.at("reconstruction_names").string())) {
        const auto key = safeKey(name);
        const auto& value = archive.at("reconstruction__" + key);
        if (value.shape.size() != 2) throw std::invalid_argument("archived readout must be outputs x modes");
        if (value.shape[1] != model.liftedDimension())
            throw std::invalid_argument("archived readout " + name + " does not match the K-ROM modes");
        ReconstructionResult loaded{value.complexDoubles(), value.shape[0], model.liftedDimension(), 0.0, 0.0, 0.0};
        if (archive.contains("reconstruction_offset__" + key)) loaded.offset = archive.at("reconstruction_offset__" + key).doubles();
        if (archive.contains("reconstruction_error__" + key)) {
            const auto error = archive.at("reconstruction_error__" + key).doubles();
            if (error.size() == 3) { loaded.residualRms = error[0]; loaded.residualMaximum = error[1]; loaded.relativeResidual = error[2]; }
        }
        model.setReconstruction(name, std::move(loaded));
    }
}

}  // namespace dynsyslab::core::detail
