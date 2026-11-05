#include "AnalyticArchive.hpp"
#include <dynsyslab/workflow/AnalyticQuery.hpp>
namespace dynsyslab::workflow::detail {
void appendAnalyticComparison(dynsyslab::core::NpzRecord &archive,
                              const dynsyslab::core::SystemDescriptor &descriptor,
                              const dynsyslab::core::SystemConfiguration &configuration,
                              const std::vector<dynsyslab::core::AnalyticSolutionFamily> &families,
                              const std::vector<double> &time, const std::vector<double> &values,
                              const std::vector<double> &initial, double initialTime,
                              const json::Value &job) {
    using namespace dynsyslab::core;
    AnalyticComparison comparison;
    const auto selection = job.at("solver").at("analytic_solution").string("auto");
    if (selection == "off") {
        comparison.status = "disabled";
        comparison.reason = "Disabled in solver settings";
    } else if (job.at("controller").at("kind").string("none") != "none") {
        comparison.reason = "No analytic family is registered for this controller law";
    } else {
        try {
            comparison =
                compareAnalyticSolution(families, initial, initialTime, time, values, selection);
        } catch (const std::exception &error) {
            comparison.status = "failed";
            comparison.reason = error.what();
        }
    }
    json::Value::Array errors;
    for (double error : comparison.maximumAbsoluteError)
        errors.emplace_back(error);
    auto solutionIdentity = analyticDescriptorJson(comparison.descriptor);
    solutionIdentity.object().erase(
        "references"); // Papers belong to the system catalogue, never solution archives.
    json::Value::Object info{
        {"schema", "dynsyslab.analytic-comparison/v1"},
        {"status", comparison.status},
        {"reason", comparison.reason},
        {"selection", selection},
        {"solution", solutionIdentity},
        {"coordinate_system", descriptor.preferredCoordinates},
        {"initial_time", initialTime},
        {"equations_fingerprint", equationsFingerprint(configuration)},
        {"maximum_absolute_error", errors},
        {"error_interpretation",
         "sampled numerical minus analytic comparison; not a certified error bound"}};
    if (!comparison.values.empty()) {
        archive.arrays["analytic_values"] =
            NpyArray::fromDoubles({time.size(), initial.size()}, comparison.values);
        archive.arrays["analytic_valid"] = NpyArray::fromDoubles({time.size()}, comparison.valid);
        info["values_array"] = "analytic_values";
        info["valid_array"] = "analytic_valid";
    }
    auto metadata = json::parse(archive.metadataJson());
    metadata.object()["analytic_solution"] = info;
    archive.arrays["metadata_json"] = NpyArray::fromString(json::dump(metadata));
}
void appendAnalyticComparison(dynsyslab::core::NpzRecord &archive,
                              const dynsyslab::core::ODESystem &system,
                              const dynsyslab::core::TimeSeries &series,
                              const std::vector<double> &initial, double initialTime,
                              const json::Value &job) {
    appendAnalyticComparison(
        archive, system.descriptor(), system.configuration(), system.analyticSolutions(),
        std::vector<double>(series.timeData(), series.timeData() + series.sampleCount()),
        std::vector<double>(series.valuesData(),
                            series.valuesData() + series.sampleCount() * series.componentCount()),
        initial, initialTime, job);
}
} // namespace dynsyslab::workflow::detail
