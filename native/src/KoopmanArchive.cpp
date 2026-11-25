#include <optional>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <algorithm>

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

NpzRecord withReconstructions(NpzRecord archive, const std::filesystem::path& path, const ModelBindingOptions& binding) {
    for (const auto& reconstructions : binding.reconstructions)
        mergeReconstructions(archive, path, readNpz(reconstructions));
    return archive;
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
    const auto archive = withReconstructions(readNpz(path),path,binding);
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

namespace {
// Catalog name for an interpolation domain: one observable, a vector observable split into
// "name[k]" channels, or a comma-joined list of scalar observables.
std::string domainObservableName(const FunctionDomain& domain) {
    const auto label=domainName(domain);
    if(!domain.stateMapFromCatalog)throw std::invalid_argument("interpolation domain "+label
        +" is not catalog observables of the state; a K-ROM cannot evaluate it from a state");
    for(const auto& name:domain.names)
        if(name.empty())throw std::invalid_argument("interpolation domain "+label
            +" has channels of unknown provenance; it cannot be evaluated from a state");
    if(domain.names.empty())throw std::invalid_argument("interpolation domain "+label+" has no observables");
    if(domain.names.size()==1)return domain.names.front();
    const auto bracket=domain.names.front().find('[');
    const auto base=domain.names.front().substr(0,bracket);
    bool vector=bracket!=std::string::npos;
    for(std::size_t k=0;vector&&k<domain.names.size();++k)
        vector=domain.names[k]==base+"["+std::to_string(k)+"]";
    if(vector)return base;
    std::string joined;
    for(const auto& name:domain.names)joined+=(joined.empty()?"":",")+name;
    return joined;
}
}  // namespace

namespace {
// A spectral K-ROM in its own frozen dictionary: psi_j(q) = sum_k E_kj phi_k(q). The domain follows
// its Hankel window: one delay of a solve's full state is the physical state; one delay of named
// coordinate columns is those catalog observables; longer windows are data coordinates (the caller
// supplies q, or asks for the interpolation domain "state" instead).
KoopmanModel spectralExpansionModel(const NpzRecord& archive,const json::Value& metadata,
                                    const SystemConfiguration& system) {
    const auto krom=readSpectralKromArchive(archive);
    const std::size_t features=krom.dictionary.size(),modes=krom.multipliers.size();
    if(!features||!modes||krom.eigenfunctions.size()!=features*modes)
        throw std::invalid_argument("spectral K-ROM eigenfunction coefficients have the wrong shape");
    std::vector<std::complex<double>> coefficients(modes*features);
    for(std::size_t j=0;j<modes;++j)for(std::size_t k=0;k<features;++k)
        coefficients[j*features+k]=krom.eigenfunctions[k*modes+j];
    const auto rates=archive.contains("continuous_eigenvalues")
        ?archive.at("continuous_eigenvalues").complexDoubles():std::vector<std::complex<double>>{};
    auto spectrum=rates.size()==modes
        ?Spectrum::fromStoredRates(rates,krom.sampleInterval,EigenvalueDomain::discreteTime)
        :Spectrum::fromDiscreteMultipliers(krom.multipliers,krom.sampleInterval);
    FunctionDomain domain;domain.dimension=krom.dictionary.dimension;
    json::Value contract;
    try{contract=json::parse(krom.coordinateContract);}catch(const std::exception&){}
    const auto& channels=contract.at("channels");
    bool oneDelay=contract.at("delays").number(0)==1&&!channels.isNull()
        &&channels.array().size()==domain.dimension;
    std::vector<std::string> observables;
    for(std::size_t k=0;oneDelay&&k<channels.array().size();++k){
        const auto& channel=channels.array()[k];
        oneDelay=channel.at("lag_samples").number(0)==0&&channel.at("measurement_column").number(-1)==double(k);
        observables.push_back(channel.at("observable").string());
    }
    const bool fromSolve=metadata.at("coordinate_transform_kind").string().empty()
        &&metadata.at("value_space").string().empty();
    const auto physical=createSystem(system)->descriptor().stateDimension;
    if(oneDelay&&fromSolve&&physical==domain.dimension)domain=FunctionDomain::physicalState(domain.dimension);
    else if(oneDelay&&!fromSolve&&std::none_of(observables.begin(),observables.end(),
                [](const std::string& name){return name.empty();})){
        domain=FunctionDomain::fromCatalog(observables,domain.dimension);
    }
    FunctionExpansion functions(std::make_shared<const ExpansionDictionary>(ExpansionDictionary::fromSpectral(krom.dictionary)),
                                std::move(coefficients),std::move(domain));
    return KoopmanModel(std::move(functions),std::move(spectrum),system);
}
}  // namespace

KoopmanModel loadKoopmanModel(const std::filesystem::path& path,
                              const SystemConfiguration& system,
                              const std::vector<std::size_t>& selectedModes,
                              const ModelBindingOptions& binding) {
    const auto archive=withReconstructions(readNpz(path),path,binding);
    if(archive.contains("coefficients"))
        return loadPolynomialKoopmanModel(path,system,selectedModes,binding);
    const auto resolved=detail::resolveArchiveBinding(archive,system,binding);
    const auto metadata=archive.metadataJson().empty()?json::Value{}:
        json::parse(archive.metadataJson());
    // A spectral K-ROM is its frozen-dictionary expansion unless an interpolation domain is requested.
    std::optional<KoopmanModel> expansion;
    if(archive.contains("spectral_krom_contract_json")&&binding.interpolationDomain.empty())
        expansion.emplace(spectralExpansionModel(archive,metadata,resolved.system));
    // Eigenfunctions whose representation is an expansion (fitted by SINDy, DMD, EDMD, or given by
    // hand): the model in that expansion, in any dictionary.
    if(!expansion&&hasEigenfunctionObjects(archive)&&binding.interpolationDomain.empty()
        &&archive.contains("expansion_coefficients")){
        const auto set=readEigenfunctionSet(archive,{"expansion"});
        for(const auto& rep:set.representations())
            if(const auto* shared=dynamic_cast<const SharedRepresentation*>(rep.get()))
                if(const auto* fx=dynamic_cast<const FunctionExpansion*>(shared->function().get()))
                    expansion=resolved.system.identifier.empty()
                        ? KoopmanModel(EigenfunctionFitResult{*fx,set.spectrum(),{}})
                        : KoopmanModel(EigenfunctionFitResult{*fx,set.spectrum(),{}},resolved.system);
    }
    InterpolationData data;
    std::vector<std::complex<double>> rates;  // the interpolated eigenfunctions' eigenvalues
    std::string coordinateObservable;
    if(expansion){
    }else if(hasEigenfunctionObjects(archive)){
        // Eigenfunction objects: their interpolation, on the domain that is part of each eigenfunction.
        // The loader does not choose the domain; it only checks a requested one.
        const auto set=readEigenfunctionSet(archive,{"interpolation"});
        const auto& interpolation=interpolationOf(set);
        const auto& domain=interpolation.coordinates();
        if(!binding.interpolationDomain.empty()&&binding.interpolationDomain!=domainName(domain))
            throw std::invalid_argument("these eigenfunctions are interpolated on "+domainName(domain)+", not "
                +binding.interpolationDomain+"; the domain belongs to the eigenfunctions and is fixed when they are made");
        if(!domain.isPhysicalState())coordinateObservable=domainObservableName(domain);
        rates=set.spectrum().continuousRates();
        data=interpolation.data();
    }else{
        coordinateObservable=metadata.at("scattered_coordinate_observable").string("");
        const auto& states=archive.at(coordinateObservable.empty()?"training_states":"training_coordinates");
        const auto& values=archive.contains("training_eigenfunction_values")
            ? archive.at("training_eigenfunction_values")
            : archive.at("eigenfunction_values");
        const auto& eigenvalues=archive.at("continuous_eigenvalues");
        if(states.shape.size()!=2||values.shape.size()!=2||eigenvalues.shape.size()!=1
            ||states.shape[0]!=values.shape[0]||values.shape[1]!=eigenvalues.shape[0])
            throw std::invalid_argument(
                "scattered KROM archive arrays have incompatible shapes");
        rates=eigenvalues.complexDoubles();
        data.trainingStates=states.doubles();data.trainingValues=values.complexDoubles();
        data.samples=states.shape[0];data.stateDimension=states.shape[1];
        data.neighbors=static_cast<std::size_t>(
            metadata.at("interpolation_neighbors").number(8.0));
        if(!metadata.at("interpolation_method").isNull())
            data.method=interpolationMethodFromJson(metadata.at("interpolation_method"),&data.neighbors);
    }
    if(!expansion&&!binding.interpolationMethod.isNull())
        data.method=interpolationMethodFromJson(binding.interpolationMethod,&data.neighbors);
    data.coordinateObservable=coordinateObservable;
    data.finiteDifferenceRelativeStep=metadata.at(
        "gradient_relative_step").number(1.0e-6);
    data.provenance=metadata.at("eigenfunction_provenance").string(
        "loaded scattered-interpolation KROM from "+path.u8string());
    KoopmanModel model=expansion?std::move(*expansion):KoopmanModel(std::move(data),Spectrum::fromContinuousRates(rates),resolved.system);
    if(!expansion&&!coordinateObservable.empty()&&archive.contains("training_states")
        &&archive.contains("training_coordinates")&&!hasEigenfunctionObjects(archive)){
        // The stored coordinates must be this observable of the stored states.
        const auto& states=archive.at("training_coordinates");
        const auto physical=archive.at("training_states").doubles();
        const auto stored=states.doubles();
        const auto value=model.coordinateObservable().evaluateReal(physical.data());
        for(std::size_t axis=0;axis<value.size();++axis)
            if(std::abs(value[axis]-stored[axis])>1e-6*(1.0+std::abs(stored[axis])))
                throw std::invalid_argument("scattered KROM coordinates do not match "+coordinateObservable);
    }
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

KoopmanModel loadKoopmanModel(const std::filesystem::path& path,const ModelBindingOptions& binding) {
    return loadKoopmanModel(path,detail::systemFromMetadata(readNpz(path)),{},binding);
}

KoopmanModel loadKoopmanModelOnDomain(const std::filesystem::path& path,const std::string& domain,
                                      const json::Value& method) {
    ModelBindingOptions binding;
    binding.interpolationDomain=domain;
    binding.interpolationMethod=method;
    return loadKoopmanModel(path,detail::systemFromMetadata(readNpz(path)),{},binding);
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
        archive.arrays[scattered->coordinateObservable.empty()?"training_states":"training_coordinates"]=
            NpyArray::fromDoubles({scattered->samples,scattered->stateDimension},scattered->trainingStates);
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
        || !model.expansion().dictionary().exactGradient() ? "finite-difference" : "analytic";
    if(const auto* scattered=model.scatteredEigenfunctions()){
        metadata.object()["interpolation_neighbors"]=
            static_cast<double>(scattered->neighbors);
        metadata.object()["gradient_relative_step"]=
            scattered->finiteDifferenceRelativeStep;
        metadata.object()["eigenfunction_provenance"]=scattered->provenance;
        metadata.object()["interpolation_method"]=interpolationMethodJson(scattered->method,scattered->neighbors);
        if(!scattered->coordinateObservable.empty())
            metadata.object()["scattered_coordinate_observable"]=scattered->coordinateObservable;
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

json::Value reconstructionsBlock(const std::filesystem::path& eigenfunctions, std::size_t modes) {
    const auto absolute = std::filesystem::absolute(eigenfunctions).lexically_normal();
    return json::Value::Object{{"basis", json::Value::Object{{"path", absolute.generic_u8string()},
                                                             {"sha256", sha256File(absolute)}}},
                               {"modes", double(modes)}};
}

void mergeReconstructions(NpzRecord& eigenfunctions, const std::filesystem::path& eigenfunctionsPath,
                          const NpzRecord& reconstructions) {
    const auto meta = reconstructions.metadataJson().empty() ? json::Value{} : json::parse(reconstructions.metadataJson());
    const auto& block = meta.at("reconstructions");
    if (block.isNull()) throw std::invalid_argument("not a reconstructions object");
    if (block.at("basis").at("sha256").string() != sha256File(eigenfunctionsPath))
        throw std::invalid_argument("these reconstructions are written over " + block.at("basis").at("path").string()
                                    + ", not " + eigenfunctionsPath.u8string());
    auto names = eigenfunctions.contains("reconstruction_names")
        ? splitNames(eigenfunctions.at("reconstruction_names").string()) : std::vector<std::string>{};
    for (const auto& name : splitNames(reconstructions.at("reconstruction_names").string())) {
        if (std::find(names.begin(), names.end(), name) != names.end())
            throw std::invalid_argument("reconstruction " + name + " is given twice");
        names.push_back(name);
    }
    for (const auto& [key, value] : reconstructions.arrays)
        if (key.rfind("reconstruction", 0) == 0 && key != "reconstruction_names") eigenfunctions.arrays[key] = value;
    eigenfunctions.arrays["reconstruction_names"] = NpyArray::fromString(joinNames(names));
}

}  // namespace dynsyslab::core
