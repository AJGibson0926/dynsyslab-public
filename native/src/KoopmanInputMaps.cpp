// A K-ROM's input maps in its archive: control_<name>_* arrays and input_map_manifest_json, each map
// bound to the identity of the eigenfunctions it acts on (koopmanModelIdentity).
#include "KoopmanInputMaps.hpp"
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <stdexcept>

namespace dynsyslab::core::detail {

void writeNamedInputMaps(NpzRecord& a, const KoopmanModel& model, bool ownOnly) {
    json::Value::Array maps;
    const auto identity = koopmanModelIdentity(model);
    for (const auto& entry : model.controlMaps()) {
        if (ownOnly && !model.ownsControlMap(entry.first)) continue;
        if (const auto* constant = std::get_if<ConstantControlMap>(&entry.second)) {
            a.arrays["control_" + entry.first + "_matrix"] =
                NpyArray::fromComplex({constant->modes, constant->inputs}, constant->matrix);
            if (!constant->discreteMatrix.empty())
                a.arrays["control_" + entry.first + "_discrete_matrix"] =
                    NpyArray::fromComplex({constant->modes, constant->inputs}, constant->discreteMatrix);
            if (constant->sampleTime > 0.0)
                a.arrays["control_" + entry.first + "_sample_time"] = NpyArray::fromDoubles({1}, {constant->sampleTime});
        } else {
            const auto& bilinear = std::get<BilinearControlMap>(entry.second);
            a.arrays["control_" + entry.first + "_operators"] =
                NpyArray::fromComplex({bilinear.inputs, bilinear.modes, bilinear.modes}, bilinear.operators);
            a.arrays["control_" + entry.first + "_offsets"] =
                NpyArray::fromComplex({bilinear.modes, bilinear.inputs}, bilinear.offsets);
            if (!bilinear.conjugateOperators.empty())
                a.arrays["control_" + entry.first + "_conjugate_operators"] =
                    NpyArray::fromComplex({bilinear.inputs, bilinear.modes, bilinear.modes}, bilinear.conjugateOperators);
        }
        maps.push_back(json::Value::Object{
            {"name", entry.first}, {"model_identity", identity},
            {"kind", std::holds_alternative<ConstantControlMap>(entry.second) ? "continuous-constant" : "continuous-bilinear"},
            {"provenance", std::visit([](const auto& m) { return m.provenance; }, entry.second)}});
    }
    std::size_t index = 0;
    for (const auto& entry : model.discreteInputMaps()) {
        if (ownOnly && !model.ownsDiscreteInputMap(entry.first)) continue;
        const auto array = "discrete_input_map_" + std::to_string(index++);
        const auto& m = entry.second;
        a.arrays[array] = NpyArray::fromComplex({m.eigenfunctionCount, m.inputDimension}, m.matrix);
        maps.push_back(json::Value::Object{{"name", entry.first}, {"model_identity", identity}, {"kind", "discrete"},
                                           {"array", array}, {"sample_interval", m.sampleInterval},
                                           {"provenance", m.provenance}});
    }
    if (!maps.empty() || !ownOnly) a.arrays["input_map_manifest_json"] = NpyArray::fromString(json::dump(maps));
}

void readNamedInputMaps(const NpzRecord& a, KoopmanModel& model) {
    if (!a.contains("input_map_manifest_json")) return;
    const auto manifest = json::parse(a.at("input_map_manifest_json").string());
    for (const auto& entry : manifest.array()) {
        const auto name = entry.at("name").string(), kind = entry.at("kind").string(), prefix = "control_" + name;
        const auto provenance = entry.at("provenance").string();
        if (kind == "discrete")
            model.setDiscreteInputMap(name, {a.at(entry.at("array").string()).complexDoubles(), model.eigenfunctionCount(),
                                             model.inputDimension(), entry.at("sample_interval").number(), provenance,
                                             entry.at("model_identity").string()});
        else if (kind == "continuous-constant") {
            ConstantControlMap m{a.at(prefix + "_matrix").complexDoubles(), model.eigenfunctionCount(),
                                 model.inputDimension(), provenance};
            m.modelIdentity = entry.at("model_identity").string();
            if (a.contains(prefix + "_discrete_matrix")) m.discreteMatrix = a.at(prefix + "_discrete_matrix").complexDoubles();
            if (a.contains(prefix + "_sample_time")) m.sampleTime = a.at(prefix + "_sample_time").doubles().at(0);
            model.setControlMap(name, std::move(m));
        } else if (kind == "continuous-bilinear") {
            BilinearControlMap m;
            m.modes = model.eigenfunctionCount();
            m.inputs = model.inputDimension();
            m.provenance = provenance;
            m.modelIdentity = entry.at("model_identity").string();
            m.operators = a.at(prefix + "_operators").complexDoubles();
            m.offsets = a.at(prefix + "_offsets").complexDoubles();
            if (a.contains(prefix + "_conjugate_operators"))
                m.conjugateOperators = a.at(prefix + "_conjugate_operators").complexDoubles();
            model.setControlMap(name, std::move(m));
        } else {
            throw std::invalid_argument("Unsupported input map contract");
        }
    }
}

}  // namespace dynsyslab::core::detail
