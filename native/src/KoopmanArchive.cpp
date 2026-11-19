#include <dynsyslab/core/KoopmanArchive.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include "KoopmanSpectrumArchive.hpp"
#include "KoopmanArchiveMetadata.hpp"
#include "KoopmanPolynomialContract.hpp"
#include "KoopmanBinding.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

const NpyArray& oneOf(const NpzRecord& record, const char* preferred, const char* fallback) {
    if (record.contains(preferred)) return record.at(preferred);
    return record.at(fallback);
}

std::vector<std::size_t> exponents(const NpyArray& source) {
    std::vector<std::size_t> result; result.reserve(source.elementCount());
    for (const double value : source.doubles()) {
        const double rounded = std::round(value);
        if (value < 0.0 || std::abs(value - rounded) > 1e-12)
            throw std::invalid_argument("polynomial exponents must be nonnegative integers");
        result.push_back(static_cast<std::size_t>(rounded));
    }
    return result;
}

std::vector<std::string> splitNames(const std::string& source) {
    std::vector<std::string> result; std::size_t begin=0;
    while (begin<source.size()) {
        const auto end=source.find('\n',begin);
        result.push_back(source.substr(begin,end==std::string::npos ? end : end-begin));
        if (end==std::string::npos) break; begin=end+1;
    }
    return result;
}

std::string joinNames(const std::vector<std::string>& names) {
    std::string result;
    for (std::size_t index=0;index<names.size();++index) {
        if (index) result.push_back('\n'); result+=names[index];
    }
    return result;
}

std::string safeKey(std::string value) {
    for (char& item : value)
        if (!std::isalnum(static_cast<unsigned char>(item)) && item!='_' && item!='-') item='_';
    return value;
}

std::vector<double> jsonNumbers(const json::Value& source) {
    std::vector<double> result;
    for (const auto& value : source.array()) result.push_back(value.number());
    return result;
}

std::vector<std::size_t> jsonIndices(const json::Value& source) {
    std::vector<std::size_t> result;
    for (const auto& value : source.array()) {
        const double number=value.number(-1.0), rounded=std::round(number);
        if (number<0.0 || std::abs(number-rounded)>1.0e-12)
            throw std::invalid_argument("legacy KROM retained indices are invalid");
        result.push_back(static_cast<std::size_t>(rounded));
    }
    return result;
}

KoopmanSINDyResult subset(KoopmanSINDyResult source,
                          const std::vector<std::size_t>& indices) {
    if (indices.empty()) return source;
    KoopmanSINDyResult result=source;
    result.continuousEigenvalues.clear(); result.coefficients.clear();
    result.trainingValues.clear(); result.residualRms.clear();
    result.residualMaximum.clear(); result.directFits=source.directFits?indices.size():0;
    result.logarithmBranches.clear();
    for (const auto index:indices) {
        if (index>=source.continuousEigenvalues.size())
            throw std::out_of_range("legacy KROM retained index exceeds spectrum");
        result.continuousEigenvalues.push_back(source.continuousEigenvalues[index]);
        if(!source.logarithmBranches.empty()) result.logarithmBranches.push_back(source.logarithmBranches.at(index));
        const auto begin=source.coefficients.begin()
            +static_cast<std::ptrdiff_t>(index*source.featureCount);
        result.coefficients.insert(result.coefficients.end(),begin,
            begin+static_cast<std::ptrdiff_t>(source.featureCount));
    }
    return result;
}

KoopmanSINDyResult nativeFit(const NpzRecord& archive) {
    const auto& coefficients=archive.at("coefficients");
    const auto& exponentArray=oneOf(archive,"polynomial_exponents","exponents");
    const auto& center=oneOf(archive,"state_center","center");
    const auto& scale=oneOf(archive,"state_scale","scale");
    const auto& eigenvalues=archive.at("continuous_eigenvalues");
    if (coefficients.shape.size()!=2 || exponentArray.shape.size()!=2
        || eigenvalues.shape.size()!=1 || center.shape.size()!=1 || scale.shape.size()!=1
        || coefficients.shape[0]!=eigenvalues.shape[0]
        || coefficients.shape[1]!=exponentArray.shape[0]
        || exponentArray.shape[1]!=center.shape[0] || scale.shape[0]!=center.shape[0])
        throw std::invalid_argument("polynomial KROM archive arrays have incompatible shapes");
    KoopmanSINDyResult fit;
    fit.continuousEigenvalues=eigenvalues.complexDoubles();
    fit.coefficients=coefficients.complexDoubles(); fit.exponents=exponents(exponentArray);
    fit.center=center.doubles(); fit.scale=scale.doubles();
    fit.stateDimension=center.shape[0]; fit.featureCount=coefficients.shape[1];
    fit.dynamicalStateDimension=archive.contains("dynamical_state_dimension")
        ? static_cast<std::size_t>(archive.at("dynamical_state_dimension").doubles().at(0))
        : fit.stateDimension;
    if (archive.contains("coordinate_observable_names"))
        fit.coordinateObservableNames=splitNames(archive.at("coordinate_observable_names").string());
    if (archive.contains("coordinate_observable_json"))
        fit.coordinateObservableDefinition=archive.at("coordinate_observable_json").string();
    if (archive.contains("coordinate_provider"))
        fit.coordinateProvider=archive.at("coordinate_provider").string();
    detail::readPolynomialSpectrumMetadata(archive,fit);
    if (archive.contains("training_residual_rms"))
        fit.residualRms=archive.at("training_residual_rms").doubles();
    if (archive.contains("training_residual_maximum"))
        fit.residualMaximum=archive.at("training_residual_maximum").doubles();
    return fit;
}

KoopmanSINDyResult legacyFit(const NpzRecord& archive, bool independent = true) {
    const auto& coefficients=archive.at("coefficients");
    const auto& eigenvalues=archive.at("eigenvalues");
    if (coefficients.shape.size()!=2 || eigenvalues.shape.size()!=1
        || coefficients.shape[0]!=eigenvalues.shape[0])
        throw std::invalid_argument("legacy polynomial KROM arrays have incompatible shapes");
    if (!archive.contains("metadata_json"))
        throw std::invalid_argument("legacy polynomial KROM requires metadata_json");
    const auto metadata=json::parse(archive.metadataJson());
    const auto& analysis=metadata.at("analysis_metadata");
    const auto& dictionary=analysis.at("dictionary");
    std::size_t dimension=0;
    if (archive.contains("training_states")
        && archive.at("training_states").shape.size()==2)
        dimension=archive.at("training_states").shape[1];
    if (!dimension)
        dimension=static_cast<std::size_t>(
            analysis.at("coordinate_observable").at("output_dimension").number());
    const auto degree=static_cast<std::size_t>(dictionary.at("degree").number());
    const bool includeConstant=dictionary.at("include_constant").boolean(true);
    KoopmanSINDyResult fit;
    const double deltaTime=analysis.at("krom").at("delta_t").number();
    if (!(deltaTime>0.0) || !std::isfinite(deltaTime))
        throw std::invalid_argument("legacy KROM requires a positive delta_t");
    fit.sampleTime=deltaTime; fit.eigenvalueOrigin="discrete-time";
    for (const auto discrete:eigenvalues.complexDoubles())
        fit.continuousEigenvalues.push_back(std::log(discrete)/deltaTime);
    fit.coefficients=coefficients.complexDoubles();
    fit.exponents=polynomialExponentTable(dimension,degree,includeConstant);
    fit.stateDimension=dimension; fit.dynamicalStateDimension=dimension;
    fit.featureCount=coefficients.shape[1];
    if (fit.exponents.size()!=fit.featureCount*dimension)
        throw std::invalid_argument("legacy KROM dictionary metadata does not match coefficients");
    const auto& preprocessing=analysis.at("state_preprocessing");
    fit.center=jsonNumbers(preprocessing.at("center"));
    fit.scale=jsonNumbers(preprocessing.at("scale"));
    const auto multiplier=jsonNumbers(preprocessing.at("multiplier"));
    if (fit.center.empty()) fit.center.assign(dimension,0.0);
    if (fit.scale.empty()) fit.scale.assign(dimension,1.0);
    if (fit.center.size()!=dimension || fit.scale.size()!=dimension)
        throw std::invalid_argument("legacy KROM state preprocessing is incomplete");
    if (!multiplier.empty()) {
        if (multiplier.size()!=dimension)
            throw std::invalid_argument("legacy KROM preprocessing multiplier is incomplete");
        for (std::size_t axis=0;axis<dimension;++axis) {
            if (!(multiplier[axis]>0.0) || !std::isfinite(multiplier[axis]))
                throw std::invalid_argument("legacy KROM preprocessing multiplier must be positive");
            fit.scale[axis]/=multiplier[axis];
        }
    }
    // Legacy archives stored explicit conjugate duals while their retained_indices
    // field identified the independent modes. Native KROMs store those independent
    // modes and form the split-real conjugate realization internally.
    if (!independent) return fit;
    return subset(std::move(fit),jsonIndices(
        analysis.at("krom").at("retained_indices")));
}

}  // namespace

KoopmanModel loadPolynomialKoopmanModel(const std::filesystem::path& path,
                                        const SystemConfiguration& system,
                                        const std::vector<std::size_t>& selectedModes,
                                        const ModelBindingOptions& binding) {
    const auto archive = readNpz(path);
    const auto resolved=detail::resolveArchiveBinding(archive,system,binding);
    KoopmanModel model=detail::hasPolynomialContract(archive)?detail::readPolynomialContract(archive,resolved.system):
        KoopmanModel(archive.contains("continuous_eigenvalues")?nativeFit(archive):legacyFit(archive),resolved.system);
    model.setBindingEvidenceJson(resolved.evidence);
    for(const auto& ref:koopmanDatasetReferencesFromMetadata(
            archive.metadataJson().empty()?json::Value{}:json::parse(archive.metadataJson())))
        model.addDatasetReference(ref);
    const std::size_t modes = model.liftedDimension(), inputs = model.inputDimension();
    if (archive.contains("control_constant_matrix")
        || archive.contains("constant_control_matrix")) {
        const auto& value = archive.contains("control_constant_matrix")
            ? archive.at("control_constant_matrix")
            : archive.at("constant_control_matrix");
        ConstantControlMap map{value.complexDoubles(),modes,inputs,
            "loaded from "+path.u8string()};
        if(archive.contains("control_constant_discrete_matrix"))
            map.discreteMatrix=archive.at("control_constant_discrete_matrix").complexDoubles();
        if(archive.contains("control_constant_sample_time"))
            map.sampleTime=archive.at("control_constant_sample_time").doubles().front();
        model.setControlMap("constant",std::move(map));
    }
    if (archive.contains("control_bilinear_operators")
        || archive.contains("bilinear_control_operators")) {
        BilinearControlMap map;
        const bool native=archive.contains("control_bilinear_operators");
        map.operators=archive.at(native?"control_bilinear_operators"
                                      :"bilinear_control_operators").complexDoubles();
        map.offsets=archive.at(native?"control_bilinear_offsets"
                                    :"bilinear_control_offsets").complexDoubles();
        const std::string conjugate=native?"control_bilinear_conjugate_operators"
            :"bilinear_control_conjugate_operators";
        if (archive.contains(conjugate))
            map.conjugateOperators=archive.at(conjugate).complexDoubles();
        map.modes = modes; map.inputs = inputs; map.provenance = "loaded from " + path.u8string();
        model.setControlMap("bilinear", std::move(map));
    }
    if (archive.contains("reconstruction_names"))
        for (const auto& name : splitNames(archive.at("reconstruction_names").string())) {
            const auto& value=archive.at("reconstruction__"+safeKey(name));
            if (value.shape.size()!=2 || value.shape[1]!=modes)
                throw std::invalid_argument("archived observable reconstruction has invalid shape");
            ReconstructionResult loaded{value.complexDoubles(),value.shape[0],modes,0.0,0.0,0.0};
            if(archive.contains("reconstruction_offset__"+safeKey(name)))
                loaded.offset=archive.at("reconstruction_offset__"+safeKey(name)).doubles();
            model.setReconstruction(name,std::move(loaded));
        }
    detail::readNamedInputMaps(archive,model);
    if(archive.contains("spectral_data_json"))model.setSpectralData(spectralDataFromJson(json::parse(archive.at("spectral_data_json").string())));
    return selectedModes.empty() ? model : model.extract(selectedModes);
}

KoopmanModel loadPolynomialKoopmanModel(
    const std::filesystem::path& path,
    const std::vector<std::size_t>& selectedModes) {
    const auto archive = readNpz(path);
    return loadPolynomialKoopmanModel(path, detail::systemFromMetadata(archive), selectedModes);
}

KoopmanModel loadKoopmanModel(const std::filesystem::path& path,
                              const SystemConfiguration& system,
                              const std::vector<std::size_t>& selectedModes,
                              const ModelBindingOptions& binding) {
    const auto archive=readNpz(path);
    if(archive.contains("coefficients"))
        return loadPolynomialKoopmanModel(path,system,selectedModes,binding);
    const auto resolved=detail::resolveArchiveBinding(archive,system,binding);
    const auto& states=archive.at("training_states");
    const auto& values=archive.contains("training_eigenfunction_values")
        ? archive.at("training_eigenfunction_values")
        : archive.at("eigenfunction_values");
    const auto& eigenvalues=archive.at("continuous_eigenvalues");
    if(states.shape.size()!=2||values.shape.size()!=2||eigenvalues.shape.size()!=1
        ||states.shape[0]!=values.shape[0]||values.shape[1]!=eigenvalues.shape[0])
        throw std::invalid_argument(
            "scattered KROM archive arrays have incompatible shapes");
    const auto metadata=archive.metadataJson().empty()?json::Value{}:
        json::parse(archive.metadataJson());
    ScatteredEigenfunctionData data;
    data.continuousEigenvalues=eigenvalues.complexDoubles();
    data.trainingStates=states.doubles();data.trainingValues=values.complexDoubles();
    data.samples=states.shape[0];data.stateDimension=states.shape[1];
    data.neighbors=static_cast<std::size_t>(
        metadata.at("interpolation_neighbors").number(8.0));
    data.finiteDifferenceRelativeStep=metadata.at(
        "gradient_relative_step").number(1.0e-6);
    data.provenance=metadata.at("eigenfunction_provenance").string(
        "loaded scattered-interpolation KROM from "+path.u8string());
    KoopmanModel model(std::move(data),resolved.system);
    model.setBindingEvidenceJson(resolved.evidence);
    for(const auto& ref:koopmanDatasetReferencesFromMetadata(metadata))
        model.addDatasetReference(ref);
    if(archive.contains("control_constant_matrix")
        ||archive.contains("constant_control_matrix")){
        const auto& value=archive.contains("control_constant_matrix")
            ?archive.at("control_constant_matrix"):archive.at("constant_control_matrix");
        ConstantControlMap map{value.complexDoubles(),model.liftedDimension(),
            model.inputDimension(),"loaded from "+path.u8string()};
        if(archive.contains("control_constant_discrete_matrix"))
            map.discreteMatrix=archive.at("control_constant_discrete_matrix").complexDoubles();
        if(archive.contains("control_constant_sample_time"))
            map.sampleTime=archive.at("control_constant_sample_time").doubles().front();
        model.setControlMap("constant",std::move(map));
    }
    if(archive.contains("control_bilinear_operators")
        ||archive.contains("bilinear_control_operators")){
        const bool native=archive.contains("control_bilinear_operators");
        BilinearControlMap map;
        map.operators=archive.at(native?"control_bilinear_operators":
            "bilinear_control_operators").complexDoubles();
        map.offsets=archive.at(native?"control_bilinear_offsets":
            "bilinear_control_offsets").complexDoubles();
        const std::string conjugate=native?"control_bilinear_conjugate_operators":
            "bilinear_control_conjugate_operators";
        if(archive.contains(conjugate))map.conjugateOperators=
            archive.at(conjugate).complexDoubles();
        map.modes=model.liftedDimension();map.inputs=model.inputDimension();
        map.provenance="loaded from "+path.u8string();
        model.setControlMap("bilinear",std::move(map));
    }
    if(archive.contains("reconstruction_names"))
        for(const auto& name:splitNames(archive.at("reconstruction_names").string())){
            const auto& value=archive.at("reconstruction__"+safeKey(name));
            ReconstructionResult loaded{value.complexDoubles(),value.shape.at(0),value.shape.at(1),0.0,0.0,0.0};
            if(archive.contains("reconstruction_offset__"+safeKey(name)))
                loaded.offset=archive.at("reconstruction_offset__"+safeKey(name)).doubles();
            model.setReconstruction(name,std::move(loaded));
        }
    detail::readNamedInputMaps(archive,model);
    if(archive.contains("spectral_data_json"))model.setSpectralData(spectralDataFromJson(json::parse(archive.at("spectral_data_json").string())));
    return selectedModes.empty()?model:model.extract(selectedModes);
}

KoopmanModel loadKoopmanModel(const std::filesystem::path& path,
                              const std::vector<std::size_t>& selectedModes) {
    const auto archive=readNpz(path);
    return loadKoopmanModel(path,detail::systemFromMetadata(archive),selectedModes);
}

KoopmanModel loadKoopmanEigenfunctionModel(const std::filesystem::path& path,
    const SystemConfiguration& system,const ModelBindingOptions& binding) {
    const auto archive = readNpz(path);
    if(detail::hasPolynomialContract(archive))return loadPolynomialKoopmanModel(path,system,{},binding);
    const auto resolved=detail::resolveArchiveBinding(archive,system,binding);
    if (archive.contains("coefficients")) {
        KoopmanModel model(archive.contains("continuous_eigenvalues") ? nativeFit(archive) : legacyFit(archive,false),resolved.system);
        model.setBindingEvidenceJson(resolved.evidence);
        for(const auto& ref:koopmanDatasetReferencesFromMetadata(
                archive.metadataJson().empty()?json::Value{}:json::parse(archive.metadataJson())))
            model.addDatasetReference(ref);
        return model;
    }
    return loadKoopmanModel(path,system,{},binding);
}

KoopmanModel loadKoopmanEigenfunctionModel(const std::filesystem::path& path) {
    return loadKoopmanEigenfunctionModel(path,detail::systemFromMetadata(readNpz(path)));
}

NpzRecord makeKoopmanArchive(const KoopmanModel& model,
                      const std::string& metadataJson,
                      std::map<std::string, NpyArray> extraArrays) {
    NpzRecord archive;archive.arrays=std::move(extraArrays);
    if(const auto* scattered=model.scatteredEigenfunctions()){
        archive.arrays["continuous_eigenvalues"]=NpyArray::fromComplex(
            {model.liftedDimension()},model.continuousEigenvalues());
        detail::writePolynomialSpectrumMetadata(archive,model.eigenfunctions());
        archive.arrays["training_states"]=NpyArray::fromDoubles(
            {scattered->samples,scattered->stateDimension},scattered->trainingStates);
        archive.arrays["eigenfunction_values"]=NpyArray::fromComplex(
            {scattered->samples,model.liftedDimension()},scattered->trainingValues);
    }else detail::writePolynomialContract(archive,model,metadataJson);
    for (const auto& item : model.controlMaps()) {
        if (const auto* constant = std::get_if<ConstantControlMap>(&item.second)) {
            archive.arrays["control_" + item.first + "_matrix"] = NpyArray::fromComplex(
                {constant->modes, constant->inputs}, constant->matrix);
            if(!constant->discreteMatrix.empty())
                archive.arrays["control_"+item.first+"_discrete_matrix"]=
                    NpyArray::fromComplex({constant->modes,constant->inputs},
                                          constant->discreteMatrix);
            if(constant->sampleTime>0.0)
                archive.arrays["control_"+item.first+"_sample_time"]=
                    NpyArray::fromDoubles({1},{constant->sampleTime});
        } else {
            const auto& bilinear = std::get<BilinearControlMap>(item.second);
            archive.arrays["control_" + item.first + "_operators"] = NpyArray::fromComplex(
                {bilinear.inputs, bilinear.modes, bilinear.modes}, bilinear.operators);
            archive.arrays["control_" + item.first + "_offsets"] = NpyArray::fromComplex(
                {bilinear.modes, bilinear.inputs}, bilinear.offsets);
            if (!bilinear.conjugateOperators.empty())
                archive.arrays["control_" + item.first + "_conjugate_operators"] = NpyArray::fromComplex(
                    {bilinear.inputs, bilinear.modes, bilinear.modes}, bilinear.conjugateOperators);
        }
    }
    detail::writeNamedInputMaps(archive,model);
    if (!model.reconstructions().empty()) {
        std::vector<std::string> names;
        for (const auto& item : model.reconstructions()) {
            names.push_back(item.first);
            archive.arrays["reconstruction__"+safeKey(item.first)]=NpyArray::fromComplex(
                {item.second.outputs,item.second.modes},item.second.matrix);
            if(!item.second.offset.empty())
                archive.arrays["reconstruction_offset__"+safeKey(item.first)]=NpyArray::fromDoubles(
                    {item.second.offset.size()},item.second.offset);
        }
        archive.arrays["reconstruction_names"]=NpyArray::fromString(joinNames(names));
    }
    json::Value metadata = metadataJson.empty()
        ? json::Value(json::Value::Object{}) : json::parse(metadataJson);
    if (!std::holds_alternative<json::Value::Object>(metadata.storage()))
        throw std::invalid_argument("KROM metadata must be a JSON object");
    auto references=model.datasetReferences();
    for(const auto& ref:koopmanDatasetReferencesFromMetadata(metadata))
        mergeKoopmanDatasetReference(references,ref);
    metadata.object()["dataset_references"]=koopmanDatasetReferencesJson(references);
    if(model.hasPhysicalSystem()) {
        metadata.object()["system_configuration"] = detail::koopmanSystemJson(model.systemConfiguration());
        metadata.object()["system"] = model.systemConfiguration().identifier;
        metadata.object()["system_preset"] = model.systemConfiguration().preset;
        metadata.object()["equations_version"] = model.systemConfiguration().equationsVersion;
    }
    metadata.object()["model_schema"]="dynsyslab.koopman-model/v2";
    metadata.object()["binding_evidence"]=json::parse(model.bindingEvidenceJson());
    metadata.object()["krom_mode_count"] = static_cast<double>(model.liftedDimension());
    metadata.object()["eigenfunction_representation"] =
        model.eigenfunctionRepresentation();
    metadata.object()["gradient_backend"] = model.scatteredEigenfunctions()
        ? "finite-difference" : "analytic";
    if(const auto* scattered=model.scatteredEigenfunctions()){
        metadata.object()["interpolation_neighbors"]=
            static_cast<double>(scattered->neighbors);
        metadata.object()["gradient_relative_step"]=
            scattered->finiteDifferenceRelativeStep;
        metadata.object()["eigenfunction_provenance"]=scattered->provenance;
    }
    archive.arrays["metadata_json"] = NpyArray::fromString(json::dump(metadata));
    archive.arrays["spectral_data_json"]=NpyArray::fromString(json::dump(spectralDataJson(model.spectralData())));
    return archive;
}

void saveKoopmanModel(const std::filesystem::path& path,const KoopmanModel& model,
    const std::string& metadataJson,std::map<std::string,NpyArray> extraArrays) {
    writeNpzAtomic(path,makeKoopmanArchive(model,metadataJson,std::move(extraArrays)));
}

void savePolynomialKoopmanModel(const std::filesystem::path& path,
                                const KoopmanModel& model,
                                const std::string& metadataJson,
                                std::map<std::string, NpyArray> extraArrays) {
    if(model.scatteredEigenfunctions())throw std::invalid_argument(
        "savePolynomialKoopmanModel received a scattered KROM");
    saveKoopmanModel(path,model,metadataJson,std::move(extraArrays));
}

}  // namespace dynsyslab::core
