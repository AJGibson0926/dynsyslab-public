#include <dynsyslab/core/EdmdControlArchive.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
EdmdControlResult loadEdmdControlModel(const std::filesystem::path& path) {
    const auto archive=readNpz(path);
    const auto metadata=json::parse(archive.metadataJson());
    if(metadata.at("schema").string()!="dynsyslab.edmdc/v1")
        throw std::invalid_argument("expected a dynsyslab.edmdc/v1 model");
    auto size=[&](const char* name) {
        const double value=metadata.at(name).number(-1);
        if(!std::isfinite(value)||value<1||value!=std::floor(value)||value>1000000000)
            throw std::invalid_argument(std::string("invalid EDMDc ")+name);
        return static_cast<std::size_t>(value);
    };
    auto array=[&](const char* name,std::vector<std::size_t> shape)->const NpyArray& {
        const auto& value=archive.at(name);
        if(value.shape!=shape) throw std::invalid_argument(std::string("invalid EDMDc array shape: ")+name);
        return value;
    };
    EdmdControlResult result;
    auto& d=result.dictionary; auto& fit=result.lifted;
    d.dimension=size("snapshot_dimension"); d.degree=size("polynomial_degree");
    d.includeConstant=metadata.at("include_constant").boolean(true);
    fit.dimension=size("feature_count"); fit.rank=size("rank"); fit.inputDimension=size("input_dimension");
    if(fit.dimension>4096||fit.rank>fit.dimension||d.degree>32)
        throw std::invalid_argument("EDMDc dimensions exceed supported dictionary bounds");
    d.center=array("state_center",{d.dimension}).doubles();
    d.scale=array("state_scale",{d.dimension}).doubles();
    for(const auto value:array("polynomial_exponents",{fit.dimension,d.dimension}).doubles()) {
        if(!std::isfinite(value)||value<0||value>32||value!=std::floor(value))
            throw std::invalid_argument("invalid EDMDc monomial exponent");
        d.exponents.push_back(static_cast<std::size_t>(value));
    }
    fit.basis=array("projection_basis",{fit.dimension,fit.rank}).doubles();
    fit.operatorMatrix=array("operator",{fit.rank,fit.rank}).doubles();
    fit.controlMatrix=array("control_matrix",{fit.rank,fit.inputDimension}).doubles();
    fit.modes=array("modes",{fit.dimension,fit.rank}).complexDoubles();
    fit.exactModes=array("exact_modes",{fit.dimension,fit.rank}).complexDoubles();
    fit.deltaTime=metadata.at("delta_t").number();
    if(!std::isfinite(fit.deltaTime)||fit.deltaTime<=0)
        throw std::invalid_argument("EDMDc model requires positive dt");
    fit.spectrum=Spectrum::fromDiscreteMultipliers(array("discrete_eigenvalues",{fit.rank}).complexDoubles(),fit.deltaTime);
    fit.hasEigenfunctions=metadata.at("eigenfunctions_available").boolean();
    if(fit.hasEigenfunctions) {
        fit.eigenfunctionCoefficients=array("eigenfunction_coefficients",{fit.rank,fit.dimension}).complexDoubles();
        fit.modalInputMatrix=array("modal_input_matrix",{fit.rank,fit.inputDimension}).complexDoubles();
    }
    using Values=const std::vector<std::complex<double>>*;
    for(const Values values:std::initializer_list<Values>{&fit.spectrum.values(),&fit.modes,&fit.exactModes,
        &fit.eigenfunctionCoefficients,&fit.modalInputMatrix})
        if(!std::all_of(values->begin(),values->end(),[](auto z){return std::isfinite(z.real())&&std::isfinite(z.imag());}))
            throw std::invalid_argument("nonfinite EDMDc spectrum or coefficients");
    fit.augmentedRank=size("augmented_rank"); fit.augmentedNumericalRank=size("augmented_numerical_rank");
    if(!metadata.at("pairs").isNull()) fit.pairs=size("pairs");
    fit.fullRowRank=metadata.at("full_row_rank").boolean();
    fit.augmentedSingularValues=archive.at("augmented_singular_values").doubles();
    fit.outputSingularValues=archive.at("output_singular_values").doubles();
    fit.oneStepRelativeError=metadata.at("one_step_relative_error").number();
    fit.projectionRelativeError=metadata.at("projection_relative_error").number();
    fit.retainedConditionNumber=metadata.at("retained_condition_number").number();
    fit.eigenvectorConditionNumber=metadata.at("eigenvector_condition_number").number();
    result.snapshotOneStepRelativeError=metadata.at("snapshot_one_step_relative_error").number();
    // Shared APIs validate dictionary and real matrices; no source read or EOM.
    (void)predictEdmdControl(result,d.center,std::vector<double>(fit.inputDimension),1);
    return result;
}
} // namespace dynsyslab::core
