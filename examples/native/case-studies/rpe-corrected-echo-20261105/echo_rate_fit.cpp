// Case-study client for the standard SINDy fitter. Predict echo-jet rates from physical state.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <fstream>
#include <iostream>

namespace core = dynsyslab::core;
namespace fs = std::filesystem;
using Json = core::json::Value;

core::TimeSeries readSeries(const core::NpzRecord& record, const std::string& name) {
    const auto& values = record.at(name);
    if (values.shape.size() != 2 || values.shape[1] != 2)
        throw std::runtime_error("Expected a two-channel archive: " + name);
    return {record.at("time").doubles(), values.doubles(), 2};
}

core::TimeSeries subsample(const core::TimeSeries& full, std::size_t stride) {
    std::vector<double> times, values;
    for (std::size_t k = 0; k < full.sampleCount(); k += stride) {
        times.push_back(full.time(k));
        values.push_back(full.value(k, 0));
        values.push_back(full.value(k, 1));
    }
    return {std::move(times), std::move(values), 2};
}

// Higher-order finite-difference control for the observation labels on the same clock.
// These four antisymmetric weights cancel the third, fifth and seventh moments.
core::TimeSeries eighthOrderLabels(const core::TimeSeries& jet, std::size_t stride) {
    const double weights[] = {4.0/5, -1.0/5, 4.0/105, -1.0/280};
    std::vector<double> times, values;
    for (std::size_t k = 4; k + 4 < jet.sampleCount(); k += stride) {
        times.push_back(jet.time(k));
        for (std::size_t component = 0; component < 2; ++component) {
            long double sum = 0;
            for (std::size_t offset = 1; offset <= 4; ++offset)
                sum += weights[offset-1] * (jet.value(k+offset, component) - jet.value(k-offset, component));
            values.push_back(static_cast<double>(sum / jet.deltaTime()));
        }
    }
    return {std::move(times), std::move(values), 2};
}

core::Observable laurentLibrary(const Json& terms) {
    const auto r = core::AnalyticExpression::state(2, 0);
    const auto v = core::AnalyticExpression::state(2, 1);
    std::vector<core::AnalyticExpression> features;
    std::vector<std::string> names;
    for (const auto& term : terms.array()) {
        const auto& powers = term.array();
        const auto p = static_cast<int>(powers.at(0).number());
        const auto q = static_cast<int>(powers.at(1).number());
        names.push_back("r^" + std::to_string(p) + " v^" + std::to_string(q));
        auto expression = core::AnalyticExpression::constant(1, 2);
        if (p) expression = expression * r.power(p);
        if (q) expression = expression * v.power(q);
        features.push_back(std::move(expression));
    }
    return core::Observable::fromExpressions(
        {"physical_state_library", "Physical-state Laurent library", "", "model", names, false},
        std::move(features));
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("Usage: echo_rate_fit job.json result.json");
        const auto job = core::json::load(fs::path(argv[1]).wstring());
        const auto stateRecord = core::readNpz(job.at("state_archive").string());
        const auto jetRecord = core::readNpz(job.at("echo_jet_archive").string());
        const auto stateMeta = core::json::parse(stateRecord.metadataJson());
        const auto jetMeta = core::json::parse(jetRecord.metadataJson());
        if (stateMeta.at("source_archive").string() != jetMeta.at("source_archive").string())
            throw std::runtime_error("State and echo must come from the same physical trajectory");
        const auto fullState = readSeries(stateRecord, "values");
        const auto fullJet = readSeries(jetRecord, "values");
        core::requireAlignedTimeSeries(fullState, fullJet);
        const auto interior = fullState.sampleCount() - 8;
        const auto stride = static_cast<std::size_t>(job.at("fit_stride").number());
        if (!stride) throw std::runtime_error("fit_stride must be positive");
        const auto inputs = subsample(core::trimTimeSeries(fullState, 4, interior), stride);
        const auto fdRates = subsample(core::trimTimeSeries(core::differentiateTimeSeries(fullJet), 2, interior), stride);
        const auto highOrderRates = eighthOrderLabels(fullJet, stride);
        const auto exactRates = subsample(core::trimTimeSeries(
            readSeries(jetRecord, "coordinate_derivatives"), 4, interior), stride);
        Json::Array fits;
        for (const auto& spec : job.at("fits").array()) {
            const auto name = spec.at("name").string();
            std::cout << name << std::endl;
            const auto library = laurentLibrary(spec.at("terms"));
            core::SindyDynamicsOptions options;
            options.threshold = spec.at("threshold").number();
            options.singularValueTolerance = spec.at("svd_cutoff").number();
            const auto labels = spec.at("labels").string();
            if (labels != "analytic" && labels != "fd8" && labels != "finite_difference")
                throw std::runtime_error("Unknown echo derivative source");
            const auto& rates = labels == "analytic" ? exactRates : (labels == "fd8" ? highOrderRates : fdRates);
            core::SindyDynamicsFitDiagnostics diagnostics;
            const auto model = core::fitSindyDynamics(inputs, rates, library, options, &diagnostics);
            fits.emplace_back(Json::Object{{"spec",spec},{"model",model.definition()},
                {"sparsity_iterations",diagnostics.toJson()}});
        }
        if (fs::exists(argv[2])) throw std::runtime_error("Output already exists");
        std::ofstream out(argv[2]);
        out << core::json::dump(Json::Object{{"inputs",Json::Array{"r","v"}},
            {"outputs",Json::Array{"de/dtau","d2e/dtau2"}},{"source_job",job},
            {"fit_rows",double(inputs.sampleCount())},{"fits",fits}});
        if (!out) throw std::runtime_error("Failed to write fit output");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
