#include "PostprocessInternal.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/Text.hpp>
#include <stdexcept>
namespace dynsyslab::workflow::postprocess_detail {
std::filesystem::path runObservableReconstruction(const std::filesystem::path& sourcePath,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    const auto krom=settings.at("krom").string();
    if(krom.empty())throw std::invalid_argument("observable reconstruction requires a KROM");
    const auto archive=dynsyslab::core::readNpz(sourcePath);
    SourceData source;source.path=sourcePath;source.metadataJson=archive.metadataJson();
    const auto array=settings.at("observable_array").string();
    const auto name=settings.at("observable").string(array);
    if(name.empty())throw std::invalid_argument("observable reconstruction requires a name");
    auto model=dynsyslab::core::loadKoopmanModel(dynsyslab::core::wide(krom));
    dynsyslab::core::ReconstructionResult fit;
    if(array.empty()) {
        source=loadSource(sourcePath);normalizeSourceCoordinates(source);
        // Loading with source physics checks system/parameter/coordinate compatibility.
        model=dynsyslab::core::loadKoopmanModel(dynsyslab::core::wide(krom),systemConfiguration(source));
        auto system=dynsyslab::core::createSystem(systemConfiguration(source));
        fit=dynsyslab::core::fitAndSetObservableReconstruction(model,name,source.values,
            source.samples,system->observableCatalog().at(name),settings.at("ridge").number());
    } else {
        // A common source archive supplies time-aligned measurements and Psi.
        // This route is usable with purely measured data and history coordinates.
        const auto& psi=archive.at(settings.at("eigenfunction_array").string("eigenfunction_values"));
        const auto& values=archive.at(array);
        if(psi.shape.size()!=2||psi.shape[1]!=model.liftedDimension()
            ||values.shape.empty()||values.shape.size()>2||values.shape[0]!=psi.shape[0])
            throw std::invalid_argument("measured reconstruction arrays are not sample aligned");
        source.samples=psi.shape[0];
        auto complexValues=[](const dynsyslab::core::NpyArray& array) {
            if(array.type==dynsyslab::core::ArrayType::complex128)return array.complexDoubles();
            const auto real=array.doubles();return std::vector<std::complex<double>>(real.begin(),real.end());
        };
        fit=dynsyslab::core::fitAndSetObservableReconstruction(model,name,complexValues(psi),
            source.samples,complexValues(values),values.shape.size()==2?values.shape[1]:1,
            settings.at("ridge").number());
    }
    const auto details="\"source_krom\":\""+escaped(krom)+"\",\"observable\":\""+escaped(name)
        +"\",\"reconstruction_fit\":\"least-squares\",\"ridge\":"+std::to_string(settings.at("ridge").number());
    auto record=dynsyslab::core::makeKoopmanArchive(model,metadata(source,"observable-reconstruction",label,details));
    record.arrays["reconstruction_fit_errors"]=dynsyslab::core::NpyArray::fromDoubles(
        {3},{fit.residualRms,fit.residualMaximum,fit.relativeResidual});
    return publish(root,source,"koopman-reconstruction",label,std::move(record));
}
} // namespace dynsyslab::workflow::postprocess_detail
