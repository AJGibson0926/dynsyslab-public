#include "CliArguments.hpp"
#include "CliUsage.hpp"
#include "CliHub.hpp"
#include "CliMcp.hpp"
#include "CliAgentSetup.hpp"
#include "CliQuickstart.hpp"
#include "CliDataPack.hpp"
#include "CliArchiveTools.hpp"
#include "CliArchive.hpp"
#include "CliExperimentReport.hpp"
#include "CliQueue.hpp"
#include <dynsyslab/workflow/JobGroup.hpp>
#include <dynsyslab/math/Query.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/Experiment.hpp>
#include <dynsyslab/workflow/ObservableQuery.hpp>
#include <dynsyslab/workflow/KoopmanSpectralQuery.hpp>
#include <dynsyslab/workflow/CoordinateClosureQuery.hpp>
#include <dynsyslab/workflow/CoordinateInverse.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/PeriodicQuery.hpp>
#include <dynsyslab/workflow/LieAlgebraQuery.hpp>
#include <dynsyslab/workflow/KromAuditQuery.hpp>
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/KromValidationPlan.hpp>
#include <dynsyslab/workflow/KromTrainingPlan.hpp>
#include <dynsyslab/workflow/Queue.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanInspection.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LegacyCsv.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/TimeSeriesArchive.hpp>
#include <dynsyslab/core/Version.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif
namespace {
namespace json = dynsyslab::core::json;
using dynsyslab::cli::Arguments;
using dynsyslab::cli::arguments;
using dynsyslab::cli::usage;
std::string readText(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to open job manifest: " + path.u8string());
    std::ostringstream out; out << stream.rdbuf(); return out.str();
}
dynsyslab::core::FourierWindow fourierWindow(const std::string& value) {
    if (value == "rectangular") return dynsyslab::core::FourierWindow::rectangular;
    if (value == "hann") return dynsyslab::core::FourierWindow::hann;
    if (value == "hamming") return dynsyslab::core::FourierWindow::hamming;
    if (value == "blackman") return dynsyslab::core::FourierWindow::blackman;
    throw std::invalid_argument("Fourier window must be rectangular, hann, hamming, or blackman");
}
void fourier(const Arguments& args) {
    if (args.positional.size() < 2)
        throw std::invalid_argument("fourier requires a trajectory NPZ archive");
    const auto archive = dynsyslab::core::readNpz(
        std::filesystem::path(args.positional[1]));
    const dynsyslab::core::TimeSeriesArchiveSelection selection{
        {"time", "controlled_time", "control_time", "free_time"},
        {"values", "states", "controlled_states", "free_training_values"}};
    const auto series = dynsyslab::core::loadTimeSeries(archive, selection);
    if (args.fourierComponent >= series.componentCount())
        throw std::out_of_range("Fourier component index exceeds the TimeSeries dimension");
    dynsyslab::core::FourierOptions options;
    options.oneSided=args.fourierOneSided;
    options.removeMean=args.fourierRemoveMean;
    options.window=fourierWindow(args.fourierWindow);
    const auto transform=series.fourier().transform(options);
    std::cout << std::setprecision(17)
        << "{\"schema\":\"dynsyslab.fourier/v1\",\"source\":"
        << json::dump(json::Value(
            std::filesystem::path(args.positional[1]).generic_u8string()))
        << ",\"component\":" << args.fourierComponent
        << ",\"sample_count\":" << series.sampleCount()
        << ",\"bin_count\":" << transform->binCount
        << ",\"delta_t\":" << transform->deltaTime
        << ",\"dominant_frequency\":"
        << transform->dominantFrequency(args.fourierComponent)
        << ",\"frequencies\":[";
    for (std::size_t bin=0;bin<transform->binCount;++bin) {
        if (bin) std::cout << ',';
        std::cout << transform->frequencies[bin];
    }
    std::cout << "],\"amplitudes\":[";
    for (std::size_t bin=0;bin<transform->binCount;++bin) {
        if (bin) std::cout << ',';
        std::cout << transform->amplitude(bin, args.fourierComponent);
    }
    std::cout << "]}\n";
}

void importData(const Arguments& args) {
    if(args.positional.size()<3||args.positional[1]!=L"legacy-csv")
        throw std::invalid_argument("import requires legacy-csv and a source file");
    dynsyslab::core::LegacyCsvImportOptions options;options.label=args.importLabel;
    options.strict=!args.allowUnknownRecords;options.inputScale=args.inputScale;
    options.coordinateSystem=args.importCoordinates;
    options.mappedSystem.identifier=args.importSystem;
    options.mappedSystem.preset=args.importPreset;
    for(const auto& item:json::parse(args.importParameters).object())
        options.mappedSystem.parameters[item.first]=item.second.number();
    const auto result=dynsyslab::core::importLegacyDynBatchCsv(
        std::filesystem::path(args.positional[2]),args.archiveRoot,options);
    (void)dynsyslab::core::indexArchiveFiles(args.archiveRoot,result.archives);
    std::cout<<json::dump(result.summary())<<'\n';
}

std::filesystem::path executablePath() {
#ifdef _WIN32
    std::wstring value(32768, L'\0');
    const DWORD length = GetModuleFileNameW(nullptr, value.data(),
                                            static_cast<DWORD>(value.size()));
    if (!length || length >= value.size())
        throw std::runtime_error("unable to resolve native executable path");
    value.resize(length);
    return std::filesystem::path(value);
#else
    return std::filesystem::canonical("/proc/self/exe");
#endif
}

void doctor(const Arguments& args) {
    const auto executable = executablePath();
    const auto root = executable.parent_path().parent_path();
    const auto start = root / L"resources" / L"docs" / L"START_HERE.md";
    const auto sdk = root / L"lib" / L"cmake" / L"DynSysLab"
        / L"DynSysLabExtensions.cmake";
    const auto computeUi = executable.parent_path() / L"DynSysLab.exe";
    const auto viewerUi = executable.parent_path() / L"DynSysLabViewer.exe";
    const auto queueModule = executable.parent_path() / L"dynsyslab_queue.dll";
    const auto icon = root / L"resources" / L"assets" / L"dynsyslab-app.ico";
    json::Value::Object paths;
    const auto describe = [](const std::filesystem::path& path) {
        json::Value::Object item;
        item["path"] = path.generic_u8string();
        item["exists"] = std::filesystem::exists(path);
        return json::Value(std::move(item));
    };
    paths["executable"] = describe(executable);
    paths["archive"] = describe(args.archiveRoot);
    paths["archive"].object()["source"] = args.archiveRootSource;
    paths["start_document"] = describe(start);
    paths["extension_sdk"] = describe(sdk);
    paths["solver_ui"] = describe(computeUi);
    paths["viewer_ui"] = describe(viewerUi);
    paths["queue_module"] = describe(queueModule);
    paths["application_icon"] = describe(icon);
    json::Value::Object result;
    result["schema"] = "dynsyslab.native-diagnostics/v1";
    result["version"] = std::string(dynsyslab::core::version);
    result["compiled"] = true;
    result["python_required"] = false;
    result["ok"] = std::filesystem::is_regular_file(executable)
        && std::filesystem::is_regular_file(start)
        && std::filesystem::is_regular_file(sdk)
        && std::filesystem::is_regular_file(computeUi)
        && std::filesystem::is_regular_file(viewerUi)
        && std::filesystem::is_regular_file(queueModule)
        && std::filesystem::is_regular_file(icon);
    result["paths"] = json::Value(std::move(paths));
    std::cout << json::dump(json::Value(std::move(result))) << '\n';
}

void docs() {
    const auto path = executablePath().parent_path().parent_path()
        / L"resources" / L"docs" / L"START_HERE.md";
    if (!std::filesystem::is_regular_file(path))
        throw std::runtime_error("installed start document is missing");
    std::wcout << path.wstring() << L'\n';
}

void requireJob(const json::Value& job) {
    dynsyslab::workflow::validateJob(job);
}

json::Value runJson(const dynsyslab::workflow::RunResult& run) {
    return dynsyslab::workflow::queueRunJson(run);
}

// --max-threads / --max-memory-mib cap the job's own `execution` limits.
json::Value execute(const json::Value& job,
                    const std::filesystem::path& archiveRoot,
                    const std::filesystem::path& progress,
                    const dynsyslab::core::ExecutionLimits& requested) {
    dynsyslab::core::ExecutionContext execution;
    execution.limits = dynsyslab::workflow::combinedLimits(
        requested, dynsyslab::workflow::jobExecutionLimits(job));
    if (execution.limits.maximumThreads)
        dynsyslab::core::setProcessThreadLimit(execution.limits.maximumThreads);
    return runJson(dynsyslab::workflow::executeJob(job, archiveRoot, progress, &execution));
}

void capabilities() {
    std::cout << "{\"schema\":\"dynsyslab.native-capabilities/v2\",\"version\":\""
              << dynsyslab::core::version
              << "\",\"compiled\":true,\"python_required\":false,\"controllers\":["
              << "\"waveform\",\"random_pulse\",\"pid\",\"lqr\",\"sdre\","
              << "\"mpc\",\"klqr\",\"kmpc\"],\"postprocessing\":["
              << "\"svd\",\"pca\",\"pod\",\"bagged_pca\",\"robust_pca\","
              << "\"hankel_pod\",\"hankel_svd\",\"hankel_matrix\",\"dmd\",\"dmdc\",\"edmdc\",\"hdmd\",\"sindy\",\"edmd\","
              << "\"koopman_spectral\",\"resdmd\",\"mpedmd\",\"hermitian_dmd\",\"tlsdmd\",\"multdmd\",\"stochastic_resdmd\",\"specrkhs\","
              << "\"havok\",\"fourier\",\"control_regression\",\"krom_extract\","
              << "\"koopman_validation\",\"boundary_eigenfunctions\","
              << "\"eigenvalue_optimization\","
              << "\"echo_history_koopman\",\"history_basis\",\"history_projection\","
              << "\"observability\",\"klqe\",\"linear_analysis\",\"linear_systems\",\"mode_selection\",\"controllability\",\"observable_reconstruction\"],\"workflows\":[\"krom_control_search\"],\"desktop_ui\":{"
              << "\"solver\":true,\"postprocessor\":true,\"viewer\":true,"
              << "\"native_jobs\":true,\"archive_tools\":true,"
              << "\"legacy_csv_import\":true,\"krom_control_search\":true,"
              << "\"initial_condition_preview\":true},"
              // CLI-only query commands (like observable, analytic, phase-space): no desktop workflow.
              << "\"queries\":[\"observable_closure\",\"lie_algebra\",\"krom_audit\",\"mathematics\"],\"initial_condition_modes\":[";
    const auto modes = dynsyslab::workflow::supportedOdeInitialConditionModes();
    for (std::size_t index = 0; index < modes.size(); ++index) {
        if (index) std::cout << ',';
        std::cout << '\"' << modes[index] << '\"';
    }
    std::cout << "],\"research_catalogue_actions\":[\"list\",\"add\"],"
              << "\"phase_space_modes\":[\"lookup\",\"compute\"],"
              << "\"analytic_solutions\":true,\"special_functions\":true}\n";
}

void systems() {
    std::cout << "{\"schema\":\"dynsyslab.system-catalog/v1\",\"systems\":[";
    const auto catalog = dynsyslab::core::availableSystems();
    for (std::size_t index = 0; index < catalog.size(); ++index) {
        if (index) std::cout << ',';
        const auto& item = catalog[index];
        std::cout << "{\"id\":\"" << item.identifier << "\",\"label\":\""
                  << item.label << "\",\"state_dimension\":" << item.stateDimension
                  << ",\"family\":\"" << item.family
                  << "\",\"coordinates\":\"" << item.preferredCoordinates << "\"}";
    }
    std::cout << "]}\n";
}

}  // namespace

int wmain(int argc, wchar_t** argv) {
    try {
        const Arguments args = arguments(argc, argv);
        if (!args.compression.empty())  // codec of every archive this run writes
            dynsyslab::core::setArchiveCompression(dynsyslab::core::parseArchiveCompression(args.compression));
        if (args.positional.empty() || args.help || args.positional.front() == L"help") { usage(); return 0; }
        const std::wstring command = args.positional.front();
        if (command == L"capabilities") capabilities();
        else if (command == L"doctor") doctor(args);
        else if (command == L"hub")
            return dynsyslab::cli::hub({args.archiveRoot, args.format, args.hook, args.peek});
        else if (command == L"mcp") return dynsyslab::cli::mcpServer();
        else if (command == L"agent-setup")
            return dynsyslab::cli::agentSetup({args.remove, args.agentCommands, args.home, args.workspace, args.format});
        else if (command == L"quickstart")
            return dynsyslab::cli::quickstart({dynsyslab::core::prepareOutputDirectory(args.archiveRoot, "--archive-root"),
                                               args.archiveRootSource, args.format});
        else if (command == L"datapack") return dynsyslab::cli::dataPackCommand(args);
        else if (command == L"docs") docs();
        else if (command == L"systems") systems();
        else if (command == L"krom-info") {
            if(args.positional.size()!=2)throw std::invalid_argument("krom-info requires a KROM archive");
            if(dynsyslab::core::readNpzIndex(args.positional[1]).contains("spectral_krom_contract_json"))
                std::cout<<json::dump(dynsyslab::workflow::evaluateKoopmanSpectralQuery(json::Value::Object{
                    {"action","inspect"},{"archive",dynsyslab::core::narrow(args.positional[1])},{"include_model",false}}))<<'\n';
            else std::cout<<json::dump(dynsyslab::core::inspectKoopmanModel(dynsyslab::core::loadKoopmanModel(args.positional[1])))<<'\n';
        } else if (command == L"krom-training-plan") {  // non-executing; see KromTrainingPlan.hpp
            if(args.positional.size()!=2)throw std::invalid_argument("krom-training-plan requires a request JSON file");
            std::cout << json::dump(dynsyslab::workflow::designKromTrainingPlan(
                json::parse(readText(args.positional[1])))) << std::endl;
        } else if (command == L"krom-validation-plan") {
            if(args.positional.size()!=2)throw std::invalid_argument("krom-validation-plan requires a KROM archive");
            std::cout << json::dump(dynsyslab::workflow::planAutonomousKoopmanValidation(
                args.positional[1],args.archiveRoot)) << '\n';
        }
        else if (command == L"research") {
            if (args.positional.size() != 2) throw std::invalid_argument("research requires a query JSON file");
            const auto query = json::parse(readText(args.positional[1]));
            std::cout << json::dump(dynsyslab::workflow::evaluateResearchCatalogueQuery(query, args.archiveRoot)) << '\n';
        }
        else if (command == L"phase-space") {
            if(args.positional.size()!=2)throw std::invalid_argument("phase-space requires a query JSON file");
            const auto query=json::parse(readText(args.positional[1]));
            std::cout << json::dump(dynsyslab::workflow::evaluatePhaseSpaceQuery(
                query, query.at("archive_root").isNull() ? args.archiveRoot : std::filesystem::path(dynsyslab::core::wide(query.at("archive_root").string())))) << '\n';
        }
        else if (command == L"analytic" || command == L"special-function") {
            if(args.positional.size()!=2)throw std::invalid_argument("query requires a JSON file");
            const auto query=json::parse(readText(args.positional[1]));
            std::cout << json::dump(command==L"analytic"
                ?dynsyslab::workflow::evaluateAnalyticQuery(query)
                :dynsyslab::workflow::evaluateSpecialFunctionQuery(query)) << '\n';
        }
        else if(command==L"koopman-spectral") {
            if(args.positional.size()!=2)throw std::invalid_argument("koopman-spectral requires a query JSON file");
            std::cout<<json::dump(dynsyslab::workflow::evaluateKoopmanSpectralQuery(
                json::parse(readText(args.positional[1]))))<<'\n';
        }
        else if (command == L"observable") {
            if(args.positional.size()!=2)throw std::invalid_argument("observable requires a query JSON file");
            std::cout << json::dump(dynsyslab::workflow::evaluateObservableQuery(
                json::parse(readText(args.positional[1])))) << '\n';
        }
        else if (command == L"krom-audit") {
            if(args.positional.size()!=2)throw std::invalid_argument("krom-audit requires a query JSON file");
            std::cout << json::dump(dynsyslab::workflow::evaluateKromAuditQuery(
                json::parse(readText(args.positional[1])))) << '\n';
        }
        else if (command == L"mathematics") {
            if(args.positional.size()!=2)throw std::invalid_argument("mathematics requires a query JSON file");
            std::cout << json::dump(dynsyslab::math::evaluateQuery(
                json::parse(readText(args.positional[1])))) << '\n';
        }
        else if (command == L"periodic-analysis") {
            if(args.positional.size()!=2)throw std::invalid_argument("periodic-analysis requires a query JSON file");
            std::cout << json::dump(dynsyslab::workflow::evaluatePeriodicQuery(
                json::parse(readText(args.positional[1])))) << '\n';
        }
        else if (command == L"lie-algebra") {
            if(args.positional.size()!=2)throw std::invalid_argument("lie-algebra requires a query JSON file");
            std::cout << json::dump(dynsyslab::workflow::evaluateLieAlgebraQuery(
                json::parse(readText(args.positional[1])))) << '\n';
        }
        else if (command == L"observable-closure") {
            if(args.positional.size()!=2)throw std::invalid_argument(
                "observable-closure requires a query JSON file");
            std::cout << json::dump(dynsyslab::workflow::evaluateCoordinateClosureQuery(
                json::parse(readText(args.positional[1])))) << '\n';
        }
        else if (command == L"materialize") {
            if(args.positional.size()!=3)throw std::invalid_argument(
                "materialize requires a recorded archive and an output NPZ");
            dynsyslab::cli::materializeArchive(args.positional[1], args.positional[2]);
        }
        else if (command == L"coordinate-inverse") {
            if(args.positional.size()!=4)throw std::invalid_argument(
                "coordinate-inverse requires a coordinate archive, an input and an output NPZ");
            const auto output = dynsyslab::core::prepareOutputFile(args.positional[3], "coordinate-inverse output");
            dynsyslab::workflow::invertCoordinates(args.positional[1],args.positional[2],output);
            std::cout << "{\"schema\":\"dynsyslab.coordinate-inverse/v1\",\"output\":"
                      << json::dump(json::Value(output.generic_u8string())) << "}\n";
        }
        else if (command == L"generate-ics") {
            if (args.positional.size() < 2)
                throw std::invalid_argument("generate-ics needs a solve job file");
            const auto job = json::parse(readText(args.positional[1]));
            if (job.at("kind").string() != "solve")
                throw std::invalid_argument("generate-ics requires a solve job");
            std::cout << json::dump(dynsyslab::workflow::odeInitialConditionSetJson(
                dynsyslab::workflow::resolveOdeInitialConditions(job))) << '\n';
        }
        else if (command == L"fourier") fourier(args);
        else if (command == L"archive") dynsyslab::cli::archiveCommand(args);
        else if (command == L"storage") dynsyslab::cli::storageCommand(args);
        else if (command == L"experiment") return dynsyslab::cli::experimentCommand(args);
        else if (command == L"import") {
            (void)dynsyslab::core::prepareOutputDirectory(args.archiveRoot, "--archive-root");
            importData(args);
        }
        else if (command == L"queue") {
            (void)dynsyslab::core::prepareOutputDirectory(args.archiveRoot, "--archive-root");
            return dynsyslab::cli::queueCommand(args);
        }
        else if (command == L"campaign") return dynsyslab::cli::campaignCommand(args);
        else if (command == L"benchmark") return dynsyslab::cli::benchmarkCommand(args);
        else if (command == L"template") return dynsyslab::cli::templateCommand(args);
        else if (command == L"validate" || command == L"plan" || command == L"run") {
            if (args.positional.size() < 2) throw std::invalid_argument("command needs a job file");
            auto job = dynsyslab::workflow::applyJobInputs(dynsyslab::workflow::bindJobFile(
                json::parse(readText(args.positional[1])), args.positional[1]), args.params, args.sets);
            dynsyslab::cli::applyJobGroupRuntime(job, args);
            requireJob(job);
            if (command == L"validate") {
                if (job.at("kind").string() == "job_group") (void)dynsyslab::workflow::planJob(job);
                std::cout << "valid\n";
            }
            else if (command == L"plan")
                std::cout << json::dump(dynsyslab::workflow::planJob(job)) << '\n';
            else {  // outputs are checked after validation and before any stage runs
                const auto root = dynsyslab::core::prepareOutputDirectory(args.archiveRoot, "--archive-root");
                const auto progress = args.progress.empty() ? std::filesystem::path{}
                    : dynsyslab::core::prepareOutputFile(args.progress, "--progress");
                std::cout << json::dump(execute(job, root, progress,
                    dynsyslab::cli::executionLimits(args))) << '\n';
            }
        } else { usage(); return 2; }
        return 0;
    } catch (const dynsyslab::core::OutputPathError& error) {
        std::cerr << json::dump(json::Value::Object{{"schema", "dynsyslab.cli-error/v1"}, {"error", error.what()},
            {"option", error.option()}, {"path", error.path().u8string()}}) << '\n';
        return 1;
    } catch (const std::exception& error) {
        std::cerr << json::dump(json::Value::Object{{"schema", "dynsyslab.cli-error/v1"}, {"error", error.what()}}) << '\n';
        return 1;
    }
}
