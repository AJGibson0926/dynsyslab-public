#include <dynsyslab/core/SpectralKromArchive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
// Bulk homogeneous arrays live in NPZ, not duplicated as decimal JSON strings.
json::Value pack(json::Value value,NpzRecord& archive,std::size_t& next) {
    if(std::holds_alternative<json::Value::Array>(value.storage())&&!value.array().empty()) {
        auto& a=value.array();
        bool flat=true;for(const auto& v:a)flat=flat&&std::holds_alternative<double>(v.storage());
        if(flat) {
            std::vector<double> data;for(const auto& v:a)data.push_back(v.number());
            auto key="spectral_array_"+std::to_string(next++);archive.arrays[key]=NpyArray::fromDoubles({data.size()},data);
            return json::Value::Object{{"$array",key}};
        }
        bool matrix=std::holds_alternative<json::Value::Array>(a.front().storage());
        std::size_t columns=matrix?static_cast<const json::Value&>(a.front()).array().size():0;
        for(const auto& row:a) {
            matrix=matrix&&row.array().size()==columns;
            for(const auto& item:row.array())matrix=matrix&&std::holds_alternative<double>(item.storage());
        }
        if(matrix&&columns) {
            std::vector<double> data;
            for(const auto& row:a)for(const auto& item:row.array())data.push_back(item.number());
            auto key="spectral_array_"+std::to_string(next++);archive.arrays[key]=NpyArray::fromDoubles({a.size(),columns},data);
            return json::Value::Object{{"$array",key}};
        }
        for(auto& v:a)v=pack(std::move(v),archive,next);
    }
    if(std::holds_alternative<json::Value::Object>(value.storage()))
        for(auto& v:value.object())v.second=pack(std::move(v.second),archive,next);
    return value;
}
json::Value unpack(json::Value value,const NpzRecord& archive) {
    if(std::holds_alternative<json::Value::Object>(value.storage())&&!static_cast<const json::Value&>(value).at("$array").isNull()) {
        const auto& array=archive.at(value.at("$array").string());
        if(array.shape.size()!=1&&array.shape.size()!=2)throw std::invalid_argument("Spectral archive array reference requires vector or matrix");
        const auto values=array.doubles();json::Value::Array a;
        if(array.shape.size()==1){for(double v:values)a.emplace_back(v);}
        else for(std::size_t i=0;i<array.shape[0];++i) {
            json::Value::Array row;for(std::size_t j=0;j<array.shape[1];++j)row.emplace_back(values[i*array.shape[1]+j]);a.emplace_back(row);
        }
        return a;
    }
    if(std::holds_alternative<json::Value::Array>(value.storage()))
        for(auto& v:value.array())v=unpack(std::move(v),archive);
    if(std::holds_alternative<json::Value::Object>(value.storage()))
        for(auto& v:value.object())v.second=unpack(std::move(v.second),archive);
    return value;
}
}
std::string spectralKromIdentity(const SpectralKROM& model) {
    // Include dictionaries, operators and readouts: changing a referenced decomposition changes identity.
    return sha256Text(json::dump(spectralKromJson(model)));
}
NpzRecord makeSpectralKromArchive(const SpectralKROM& model,const std::string& metadataJson) {
    NpzRecord a;std::size_t next=0;
    a.arrays["spectral_krom_contract_json"]=NpyArray::fromString(json::dump(pack(spectralKromJson(model),a,next)));
    auto metadata=json::parse(metadataJson);metadata.at("schema")="dynsyslab.spectral-krom/v1";
    metadata.at("model_identity")=spectralKromIdentity(model);metadata.at("method")=model.method;
    metadata.at("delta_t")=model.sampleInterval;metadata.at("eigenvalue_domain")="discrete-time";
    a.arrays["metadata_json"]=NpyArray::fromString(json::dump(metadata));
    a.arrays["discrete_eigenvalues"]=NpyArray::fromComplex({model.rank},model.multipliers);
    a.arrays["eigenvalues"]=a.arrays["discrete_eigenvalues"];
    return a;
}
SpectralKROM readSpectralKromArchive(const NpzRecord& a) {
    auto model=spectralKromFromJson(unpack(json::parse(a.at("spectral_krom_contract_json").string()),a));
    const auto expected=json::parse(a.metadataJson()).at("model_identity").string();
    if(expected.empty()||expected!=spectralKromIdentity(model))throw std::invalid_argument("Spectral KROM identity mismatch");
    return model;
}
SpectralKROM loadSpectralKrom(const std::filesystem::path& path){return readSpectralKromArchive(readNpz(path));}
Observable spectralKromObservable(SpectralKROM model,ObservableDecompositionReference reference,std::string name) {
    model.validate();
    if(reference.modelIdentity.empty()||reference.modelIdentity!=spectralKromIdentity(model))
        throw std::invalid_argument("Observable reference does not match spectral KROM identity");
    std::size_t outputs=1;
    if(reference.kind=="reconstruction")outputs=model.reconstructions.at(reference.component).outputs;
    else if(reference.kind=="eigenfunction") {
        if(reference.index>=model.rank)throw std::invalid_argument("Eigenfunction reference index out of range");
    } else if(reference.kind=="wave-packet") {
        if(reference.index>=model.spectral.wavePackets.size())throw std::invalid_argument("Wave packet reference index out of range");
        outputs=model.spectral.wavePackets[reference.index].grid.size();
    } else throw std::invalid_argument("Spectral measures are associated data, not pointwise observables");
    ObservableDescriptor descriptor;descriptor.name=std::move(name);descriptor.label=descriptor.name;
    descriptor.stateCoordinateSystem=model.coordinateContract;descriptor.complexValued=true;descriptor.decompositions={reference};
    const auto dimension=model.dictionary.dimension;
    return Observable(descriptor,dimension,outputs,[model=std::move(model),reference,dimension](const double* point) {
        std::vector<double> x(point,point+dimension);Observable::Values out;
        if(reference.kind=="reconstruction") {for(double v:model.reconstruct(reference.component,x))out.emplace_back(v);}
        else if(reference.kind=="eigenfunction")out.push_back(model.evaluateEigenfunctions(x)[reference.index]);
        else {
            const auto& p=model.spectral.wavePackets[reference.index];const auto values=model.dictionary.evaluate(x);
            out.resize(p.grid.size());
            for(std::size_t i=0;i<out.size();++i)for(std::size_t j=0;j<values.size();++j)out[i]+=p.coefficients[i*values.size()+j]*values[j];
        }
        return out;
    });
}
Observable existingKromObservable(const std::filesystem::path& path,ObservableDecompositionReference ref,std::string name) {
    const auto archive=readNpz(path);ref.artifact=path.generic_u8string();
    if(archive.contains("spectral_krom_contract_json"))return spectralKromObservable(readSpectralKromArchive(archive),ref,std::move(name));
    auto model=loadKoopmanEigenfunctionModel(path);
    if(ref.modelIdentity.empty()||ref.modelIdentity!=koopmanModelIdentity(model))throw std::invalid_argument("Observable reference does not match existing KROM identity");
    if(ref.kind!="reconstruction")throw std::invalid_argument("Legacy KROM observable binding currently requires named reconstruction");
    if(model.requiresExternalCoordinates())throw std::invalid_argument("Legacy KROM reference needs its external coordinate provider");
    ObservableDescriptor descriptor;descriptor.name=std::move(name);descriptor.label=descriptor.name;descriptor.complexValued=true;descriptor.decompositions={ref};
    return Observable(descriptor,model.stateDimension(),model.reconstruction(ref.component).outputs,
        [model=std::move(model),ref](const double* x){return model.reconstruct(ref.component,x);});
}
} // namespace dynsyslab::core
