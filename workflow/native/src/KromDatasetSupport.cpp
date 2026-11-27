#include <set>
#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
std::string datasetDynamics(const SourceData& source) {
    if(std::any_of(source.inputs().begin(),source.inputs().end(),[](double u){return !std::isfinite(u)||u!=0.;}))
        return "forced";
    const auto meta=source.metadataJson.empty()?json::Value{}:json::parse(source.metadataJson);
    const auto forcing=meta.at("forcing_type").string();
    const auto controller=meta.at("controller_kind").string(meta.at("controller_configuration").at("kind").string());
    if((!controller.empty()&&controller!="none") || (!forcing.empty()&&forcing!="free"&&forcing!="unforced"&&forcing!="none"&&forcing!="unknown"))
        return "forced";
    if(!source.inputs().empty() || controller=="none" || forcing=="free" || forcing=="unforced")
        return "unforced";
    return "unknown";
}

std::string trainingDataRole(const std::string& method, const json::Value& settings) {
    static const std::set<std::string> autonomous{
        "dmd", "hdmd", "hankel_dmd", "sindy", "edmd", "havok", "koopman_spectral", "resdmd", "mpedmd",
        "hermitian_dmd", "tlsdmd", "multdmd", "stochastic_resdmd", "specrkhs", "boundary_eigenfunctions",
        "eigenvalue_optimization", "echo_history_koopman", "history_basis", "sindy_dynamics", "koopman_validation"};
    if (method == "dmdc" || method == "edmdc") return "forced_allowed";
    if (method == "control_regression" || method == "history_control_regression"
        || method == "bilinear_certificate") return "forced_required";
    if (method == "input_map") {
        const auto& map = settings.at("map");
        const auto representation = settings.at("representation").string(map.at("representation").string());
        return representation == "constant" || representation == "bilinear" ? "forced_required" : "any";
    }
    return autonomous.count(method) ? "autonomous" : "any";
}

void requireTrainingDataRole(const std::string& method, const SourceData& source, const json::Value& settings) {
    const auto role = trainingDataRole(method, settings);
    if (role == "any" || role == "forced_allowed") return;
    const auto dynamics = datasetDynamics(source);
    const auto where = " (" + source.path.filename().u8string() + ")";
    if (role == "forced_required" && dynamics != "forced")
        throw std::invalid_argument(method + " regresses an input map and needs forced training data with a "
                                    "recorded control history; this source is " + dynamics + where);
    if (role == "autonomous") {
        if (dynamics == "forced")
            throw std::invalid_argument(method + " extracts autonomous Koopman eigenfunctions/spectra and must be trained "
                "on unforced data; this source is forced" + where + ". Use dmdc/edmdc to identify from forced data, "
                "or input_map to regress B_psi from it.");
        // "unknown" (no forcing provenance: legacy imports, synthetic or measured records) is accepted;
        // settings.training_dynamics = "forced" lets a stage declare such a source forced and be refused.
        if (dynamics == "unknown" && settings.at("training_dynamics").string() == "forced")
            throw std::invalid_argument(method + " must be trained on unforced data; the stage declares this source forced" + where);
    }
}

void requireTrainingDataRole(const std::string& method, const std::filesystem::path& path, const json::Value& settings) {
    const auto role = trainingDataRole(method, settings);
    if (role == "any" || role == "forced_allowed") return;
    const auto archive = dynsyslab::core::readNpz(path);
    SourceData light; light.path = path;
    light.metadataJson = archive.metadataJson();
    for (const char* name : {"control_interval_values", "control_values"})
        if (archive.contains(name)) {  // only the input values are read: index them by row
            auto inputs = archive.at(name).doubles();
            std::vector<double> rows(inputs.size());
            for (std::size_t k = 0; k < rows.size(); ++k) rows[k] = double(k);
            light.record.inputs = dynsyslab::core::TimeSeries(std::move(rows), std::move(inputs), 1);
            break;
        }
    requireTrainingDataRole(method, light, settings);
}

void requireMatchingUnforcedDataset(const dynsyslab::core::KoopmanModel& model,const SourceData& source) {
    if(datasetDynamics(source)!="unforced")
        throw std::invalid_argument("Autonomous Koopman validation requires a confirmed unforced dataset; forced/unknown data cannot be scored with the free semigroup law");
    const auto actual=systemConfiguration(source);
    dynsyslab::core::requireCompatibleSystems(model.systemConfiguration(),actual);
    if(model.requiresExternalCoordinates())
        throw std::invalid_argument("History-coordinate KROM validation requires the frozen coordinate provider; raw dynamical states are not latent coordinates");
    const auto meta=json::parse(source.metadataJson);
    const auto coordinates=meta.at("archive_coordinate_system").string(meta.at("integration_form").string());
    const auto preferred=dynsyslab::core::createSystem(model.systemConfiguration())->descriptor().preferredCoordinates;
    if(!coordinates.empty() && coordinates!="automatic" && coordinates!=preferred)
        throw std::invalid_argument("Validation coordinates differ from the KROM system coordinates");
}

SourceData validationSamples(const SourceData& source,const dynsyslab::core::KoopmanModel& model,const json::Value& settings) {
    const double sourceDt=uniformDeltaTime(source);
    const double requested=settings.at("validation_delta_t").number(
        model.spectrum().sampleInterval()>0?model.spectrum().sampleInterval():sourceDt);
    const double ratio=requested/sourceDt;
    if(!std::isfinite(ratio)||ratio<1-1e-8||ratio>1e9||std::abs(ratio-std::round(ratio))>1e-7*std::max(1.,ratio))
        throw std::invalid_argument("Validation dt must be an integer multiple of the source dt; explicitly resample separately if necessary");
    const auto stride=static_cast<std::size_t>(std::llround(ratio));
    const auto& spectrum=model.spectrum();
    if(spectrum.domain()==dynsyslab::core::EigenvalueDomain::discreteTime) {
        const double multiple=requested/spectrum.sampleInterval();
        if(multiple<1-1e-8 || std::abs(multiple-std::round(multiple))>1e-7)
            throw std::invalid_argument("Sampled-origin eigenfunctions must be validated at their fitted dt or an integer multiple");
    }
    const double start=settings.at("start").number(source.time().front());
    const double stop=settings.at("stop").number(source.time().back());
    if(!std::isfinite(start)||!std::isfinite(stop)||!(stop>start))
        throw std::invalid_argument("Validation interval requires finite start < stop");
    std::vector<double> time,values;
    for(std::size_t k=0;k<source.samples();k+=stride) {
        if(source.time()[k]<start-1e-10 || source.time()[k]>stop+1e-10)continue;
        time.push_back(source.time()[k]);
        values.insert(values.end(),source.values().begin()+k*source.dimension(),
            source.values().begin()+(k+1)*source.dimension());
    }
    auto result=withSeries(source,std::move(time),std::move(values),source.dimension());
    result.record.inputs.reset();result.record.inputsOnIntervals=false;
    if(result.samples()<3)throw std::invalid_argument("Validation selection has fewer than three samples");
    return result;
}
}
