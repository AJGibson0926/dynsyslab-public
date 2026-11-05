#include <dynsyslab/core/KoopmanInspection.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
namespace dynsyslab::core {
json::Value inspectKoopmanModel(const KoopmanModel &model) {
    const auto &s = model.spectrum();
    const auto &f = model.fitDiagnostics();
    json::Value::Array spectrum, maps;
    for (std::size_t j = 0; j < s.size(); ++j) {
        const auto primary = s.values()[j];
        const auto rate = s.continuousRate(j);
        spectrum.emplace_back(json::Value::Object{
            {"index", double(j)},
            {"value", json::Value::Array{primary.real(), primary.imag()}},
            {"continuous_rate",
             rate ? json::Value(json::Value::Array{rate->real(), rate->imag()}) : json::Value{}}});
    }
    for (const auto &entry : model.controlMaps())
        maps.emplace_back(
            json::Value::Object{{"name", entry.first}, {"time_domain", "continuous-time"}});
    for (const auto &entry : model.discreteInputMaps())
        maps.emplace_back(json::Value::Object{{"name", entry.first},
                                              {"time_domain", "discrete-time"},
                                              {"sample_interval", entry.second.sampleInterval}});
    json::Value::Object result{
        {"schema", "dynsyslab.krom-info/v1"},
        {"model_identity",koopmanModelIdentity(model)},
        {"binding_evidence",json::parse(model.bindingEvidenceJson())},
        {"representation", model.eigenfunctionRepresentation()},
        {"eigenfunction_count", double(model.eigenfunctionCount())},
        {"real_lifted_dimension", double(model.realLiftedDimension())},
        {"coordinates", eigenfunctionCoordinatesJson(model.coordinates())},
        {"physical_system_attached", model.hasPhysicalSystem()},
        {"input_dimension", double(model.inputDimension())},
        {"spectrum_domain",
         s.domain() == EigenvalueDomain::discreteTime ? "discrete-time" : "continuous-time"},
        {"sample_interval", s.sampleInterval()},
        {"spectrum", spectrum},
        {"input_maps", maps},
        {"dataset_references", koopmanDatasetReferencesJson(model.datasetReferences())},
        {"fit_evidence", json::Value::Object{{"algorithm", f.algorithm},
                                             {"sample_count", double(f.sampleCount)},
                                             {"residual_definition", f.residualDefinition},
                                             {"residuals_available", !f.residualRms.empty()}}},
        {"validation", "Dataset pointers do not imply successful validation. Read the separate "
                       "validation result."}};
    if (model.physicalStateDimension())
        result["physical_state_dimension"] = double(*model.physicalStateDimension());
    else
        result["physical_state_dimension"] = json::Value{};
    if (!model.scatteredEigenfunctions())
        result["polynomial_feature_count"] =
            double(model.polynomialEigenfunctions().dictionary().featureCount());
    return result;
}
} // namespace dynsyslab::core
