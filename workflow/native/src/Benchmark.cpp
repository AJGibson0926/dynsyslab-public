#include <dynsyslab/workflow/Benchmark.hpp>

#include <dynsyslab/core/BuildInfo.hpp>
#include <dynsyslab/core/Decomposition.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/JobStore.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/System.hpp>

#include "WorkflowSupport.hpp"

#include <algorithm>
#include <chrono>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <thread>

namespace dynsyslab::workflow {
namespace {

using Kernel = std::function<void()>;

struct Case { std::string name; std::string category; Kernel run; };

dynsyslab::core::NativeSolution duffing(dynsyslab::core::IntegrationMethod method,
                                        double stop, double step) {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing";
    configuration.preset = "double_well_reference";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = stop; options.outputStep = step; options.method = method;
    options.relativeTolerance = 1e-9; options.absoluteTolerance = 1e-12;
    return dynsyslab::core::integrate(*system, {-1.0, 0.1}, 1, options);
}

Case ode(const std::string& name, dynsyslab::core::IntegrationMethod method) {
    return {name, "ode", [method] { (void)duffing(method, 20.0, 0.005); }};
}

const std::vector<Case>& suite() {
    using Method = dynsyslab::core::IntegrationMethod;
    static const std::vector<Case> cases = [] {
        std::vector<Case> result{ode("ode_rk4_duffing", Method::rk4),
            ode("ode_rk45_duffing", Method::dormandPrince45),
            ode("ode_dop853_duffing", Method::dormandPrince853),
            ode("ode_bdf_duffing", Method::bdf)};
        result.push_back({"ising_2d_32", "lattice", [] {
            dynsyslab::core::IsingRequest request;
            request.size = 32; request.temperature = 2.269; request.seed = 20260819;
            request.sweeps = 200;
            request.spins.assign(32 * 32, 1);
            (void)dynsyslab::core::runIsing2D(request);
        }});
        result.push_back({"edmd_duffing_degree4", "decomposition", [] {
            static const auto trajectory = duffing(Method::dormandPrince45, 20.0, 0.005);
            const std::size_t samples = trajectory.sampleCount;
            const std::vector<double>& states = trajectory.values;
            dynsyslab::core::EDMDOptions options;
            options.dictionary.polynomialDegree = 4;
            (void)dynsyslab::core::extendedDynamicModeDecomposition(states, samples, 2, 0.005, options);
        }});
        return result;
    }();
    return cases;
}

json::Value runtimeFingerprint() {
    const auto& build = dynsyslab::core::buildIdentity();
    return json::Value::Object{{"version", build.version}, {"source_commit", build.sourceCommit},
        {"source_modified", build.sourceModified}, {"compiler", build.compiler},
        {"build_type", build.buildType}, {"platform", build.platform},
        {"machine", dynsyslab::core::currentHost()},
        {"hardware_threads", static_cast<double>(std::thread::hardware_concurrency())}};
}

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    return values.size() % 2 ? values[middle] : 0.5 * (values[middle - 1] + values[middle]);
}

}  // namespace

std::vector<std::string> benchmarkCaseNames() {
    std::vector<std::string> names;
    for (const auto& item : suite()) names.push_back(item.name);
    return names;
}

json::Value runBenchmarks(const BenchmarkOptions& options) {
    if (!options.iterations) throw std::invalid_argument("benchmark iterations must be positive");
    std::vector<const Case*> selected;
    for (const auto& item : suite())
        if (options.cases.empty()
            || std::find(options.cases.begin(), options.cases.end(), item.name) != options.cases.end())
            selected.push_back(&item);
    for (const auto& name : options.cases)
        if (std::none_of(selected.begin(), selected.end(),
                         [&](const Case* item) { return item->name == name; }))
            throw std::invalid_argument("unknown benchmark case: " + name);
    json::Value::Array samples;
    for (const auto* item : selected) {
        for (std::size_t index = 0; index < options.warmup; ++index) item->run();
        std::vector<double> seconds;
        json::Value::Array raw;
        for (std::size_t index = 0; index < options.iterations; ++index) {
            const auto started = std::chrono::steady_clock::now();
            item->run();
            seconds.push_back(std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count());
            raw.emplace_back(seconds.back());
        }
        samples.emplace_back(json::Value::Object{{"name", item->name}, {"category", item->category},
            {"seconds", json::Value(std::move(raw))},
            {"minimum_seconds", *std::min_element(seconds.begin(), seconds.end())},
            {"median_seconds", median(seconds)},
            {"maximum_seconds", *std::max_element(seconds.begin(), seconds.end())}});
    }
    return json::Value::Object{{"schema", "dynsyslab.benchmark-report/v1"},
        {"created_utc", detail::timestamp()}, {"runtime", runtimeFingerprint()},
        {"iterations", static_cast<double>(options.iterations)},
        {"warmup", static_cast<double>(options.warmup)}, {"samples", json::Value(std::move(samples))}};
}

json::Value benchmarkGate(const json::Value& reference, const json::Value& candidate,
                          const json::Value& required, bool requireSameRuntime) {
    for (const auto* report : {&reference, &candidate})
        if (report->at("schema").string() != "dynsyslab.benchmark-report/v1")
            throw std::invalid_argument("benchmark gates need dynsyslab.benchmark-report/v1 reports");
    json::Value::Object differences;
    for (const char* key : {"machine", "platform", "compiler", "build_type"}) {
        const auto& first = reference.at("runtime").at(key);
        const auto& second = candidate.at("runtime").at(key);
        if (json::dump(first) != json::dump(second))
            differences[key] = json::Value(json::Value::Array{first, second});
    }
    const auto medians = [](const json::Value& report) {
        std::map<std::string, double> result;
        for (const auto& sample : report.at("samples").array())
            result[sample.at("name").string()] = sample.at("median_seconds").number();
        return result;
    };
    const auto before = medians(reference), after = medians(candidate);
    std::map<std::string, double> requirements;
    for (const auto& [name, value] : required.object()) {
        if (!(value.number() > 0.0)) throw std::invalid_argument("required speedups must be positive");
        if (name == "*") {
            for (const auto& [shared, ignored] : before)
                if (after.count(shared) && !required.object().count(shared))
                    requirements[shared] = value.number();
        } else requirements[name] = value.number();
    }
    if (requirements.empty()) throw std::invalid_argument("a benchmark gate needs at least one requirement");
    bool passed = differences.empty() || !requireSameRuntime;
    const bool comparable = passed;
    json::Value::Array results;
    for (const auto& [name, speedupRequired] : requirements) {
        if (!before.count(name) || !after.count(name))
            throw std::invalid_argument("benchmark " + name + " must be in both reports");
        const double reference_ = before.at(name), candidate_ = after.at(name);
        const double speedup = candidate_ == 0.0 ? std::numeric_limits<double>::max()
                                                 : reference_ / candidate_;
        const bool ok = speedup >= speedupRequired;
        passed = passed && ok;
        results.emplace_back(json::Value::Object{{"name", name}, {"reference_seconds", reference_},
            {"candidate_seconds", candidate_}, {"speedup", speedup},
            {"required_speedup", speedupRequired}, {"passed", ok}});
    }
    return json::Value::Object{{"schema", "dynsyslab.performance-gate/v1"},
        {"comparable_runtime", comparable}, {"runtime_differences", json::Value(std::move(differences))},
        {"passed", passed}, {"results", json::Value(std::move(results))}};
}

}  // namespace dynsyslab::workflow
