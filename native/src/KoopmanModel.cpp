#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

void finiteMapValues(const std::vector<std::complex<double>>& values) {
    for(auto z:values)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))
        throw std::invalid_argument("Control map coefficients must be finite");
}

bool realMode(std::complex<double> value) { return std::abs(value.imag()) <= 1e-12; }

void validate(const ConstantControlMap& value, std::size_t modes, std::size_t inputs) {
    finiteMapValues(value.matrix);finiteMapValues(value.discreteMatrix);
    if (value.modes != modes || value.inputs != inputs || value.matrix.size() != modes * inputs
        ||(!value.discreteMatrix.empty()&&value.discreteMatrix.size()!=modes*inputs)
        ||!std::isfinite(value.sampleTime)||value.sampleTime<0.0)
        throw std::invalid_argument("constant control map dimensions do not match KROM");
}

void validate(const BilinearControlMap& value, std::size_t modes, std::size_t inputs) {
    finiteMapValues(value.operators);finiteMapValues(value.offsets);finiteMapValues(value.conjugateOperators);
    if (value.modes != modes || value.inputs != inputs
        || value.operators.size() != inputs * modes * modes
        || value.offsets.size() != modes * inputs
        || (!value.conjugateOperators.empty()
            && value.conjugateOperators.size() != value.operators.size()))
        throw std::invalid_argument("bilinear control map dimensions do not match KROM");
}

std::vector<std::complex<double>> evaluateStored(
    const ApproximateControlMap& stored,
    const std::vector<std::complex<double>>& psi) {
    if (const auto* constant = std::get_if<ConstantControlMap>(&stored)) return constant->matrix;
    const auto& bilinear = std::get<BilinearControlMap>(stored);
    std::vector<std::complex<double>> result = bilinear.offsets;
    for (std::size_t input = 0; input < bilinear.inputs; ++input)
        for (std::size_t row = 0; row < bilinear.modes; ++row)
            for (std::size_t column = 0; column < bilinear.modes; ++column) {
                const auto at = (input * bilinear.modes + row) * bilinear.modes + column;
                result[row * bilinear.inputs + input] += bilinear.operators[at] * psi[column];
                if (!bilinear.conjugateOperators.empty())
                    result[row * bilinear.inputs + input] +=
                        bilinear.conjugateOperators[at] * std::conj(psi[column]);
            }
    return result;
}

}  // namespace

KoopmanModel::KoopmanModel(KoopmanSINDyResult eigenfunctions,
                           SystemConfiguration system)
    : KoopmanModel(polynomialEigenfunctionsFromLegacy(eigenfunctions),std::move(system),
                   fitDiagnosticsFromLegacy(eigenfunctions)) {
    eigenfunctions.dynamicalStateDimension=stateDimension();
    eigenfunctions_=std::move(eigenfunctions);
}

KoopmanModel::KoopmanModel(ScatteredEigenfunctionData eigenfunctions,
                           SystemConfiguration system)
    : systemConfiguration_(std::move(system)),
      scattered_(std::make_shared<ScatteredEigenfunctionSet>(
          std::move(eigenfunctions))) {
    system_=std::shared_ptr<ODESystem>(createSystem(systemConfiguration_).release());
    systemConfiguration_=system_->configuration();
    if(scattered_->data().stateDimension!=system_->descriptor().stateDimension)
        throw std::invalid_argument(
            "scattered KROM eigenfunctions and dynamical system differ in dimension");
    eigenfunctions_.continuousEigenvalues=scattered_->data().continuousEigenvalues;
    eigenfunctions_.stateDimension=scattered_->data().stateDimension;
    eigenfunctions_.dynamicalStateDimension=scattered_->data().stateDimension;
    scatteredSpectrum_=EigenfunctionSpectrum::fromContinuousRates(eigenfunctions_.continuousEigenvalues);
    scatteredCoordinates_.kind=EigenfunctionCoordinateKind::physicalState;
    scatteredCoordinates_.dimension=eigenfunctions_.stateDimension;
    scatteredCoordinates_.physicalStateDimension=eigenfunctions_.stateDimension;
    coordinateObservable_=std::make_shared<Observable>(
        system_->observableCatalog().at("state"));
}

std::size_t KoopmanModel::stateDimension() const {
    const auto dimension=physicalStateDimension();
    if(!dimension)throw std::logic_error("Physical state dimension is unknown; use coordinateDimension()");
    return *dimension;
}
std::size_t KoopmanModel::coordinateDimension() const noexcept {
    return coordinates().dimension;
}
std::size_t KoopmanModel::liftedDimension() const noexcept {
    return spectrum().size();
}
std::size_t KoopmanModel::splitDimension() const noexcept {
    std::size_t result = 0;
    for (std::size_t j=0;j<liftedDimension();++j) result += isRealMode(j) ? 1 : 2;
    return result;
}
std::size_t KoopmanModel::inputDimension() const {
    if(!system_)return suppliedInputDimension_;
    const auto initial = system_->initialState();
    const auto matrix = system_->inputMatrix(initial.data());
    if (matrix.empty() || matrix.size() % stateDimension())
        throw std::runtime_error("system input matrix has invalid dimensions");
    return matrix.size() / stateDimension();
}
const KoopmanSINDyResult& KoopmanModel::eigenfunctions() const {
    if(!spectrum().hasContinuousRates())throw std::logic_error("Legacy eigenfunctions() requires finite continuous rates; use polynomialEigenfunctions()");
    return eigenfunctions_;
}
const std::vector<std::complex<double>>&
KoopmanModel::continuousEigenvalues() const {
    return spectrum().continuousRates();
}
std::string KoopmanModel::eigenfunctionRepresentation() const {
    return scattered_ ? "scattered-interpolation" : "analytic-polynomial";
}
std::vector<KoopmanEigenfunctionDescriptor>
KoopmanModel::eigenfunctionDescriptors() const {
    std::vector<KoopmanEigenfunctionDescriptor> result;
    for(std::size_t index=0;index<liftedDimension();++index)
        result.push_back({index,spectrum().continuousRate(index).value_or(std::complex<double>{}),
            eigenfunctionRepresentation(),scattered_?"finite-difference":"analytic",
            spectrum().continuousRate(index).has_value()});
    return result;
}
const ScatteredEigenfunctionData*
KoopmanModel::scatteredEigenfunctions() const noexcept {
    return scattered_?&scattered_->data():nullptr;
}
const SystemConfiguration& KoopmanModel::systemConfiguration() const {
    if(!system_)throw std::logic_error("Coordinate-only KROM has no physical-system adapter");
    return systemConfiguration_;
}

std::vector<std::complex<double>> KoopmanModel::evaluate(const double* state) const {
    if(scattered_)return scattered_->evaluate(state);
    if(coordinates().kind==EigenfunctionCoordinateKind::physicalState)return evaluateCoordinates(state);
    if (!coordinateObservable_)
        throw std::logic_error(
            "causal-history KROM requires evaluateCoordinates after updating its history projector");
    const auto coordinates=coordinateObservable_->evaluateReal(state);
    return evaluateCoordinates(coordinates.data());
}

const Observable& KoopmanModel::coordinateObservable() const {
    if (!coordinateObservable_)
        throw std::logic_error(
            "causal-history KROM has no instantaneous state observable");
    return *coordinateObservable_;
}

bool KoopmanModel::requiresExternalCoordinates() const noexcept {
    return !scattered_ && coordinates().requiresExternalValues();
}

std::vector<std::complex<double>> KoopmanModel::evaluateCoordinates(const double* coordinates) const {
    if(scattered_)throw std::logic_error("direct coordinate evaluation requires polynomial eigenfunctions");
    return polynomial_->evaluate(coordinates);
}

std::vector<double> KoopmanModel::evaluateSplit(const double* state) const {
    return splitValues(evaluate(state));
}

std::vector<double> KoopmanModel::evaluateCoordinatesSplit(const double* coordinates) const {
    return splitValues(evaluateCoordinates(coordinates));
}

std::vector<double> KoopmanModel::splitValues(const std::vector<std::complex<double>>& values) const {
    if(values.size()!=eigenfunctionCount())throw std::invalid_argument("Eigenfunction value count differs from KROM");
    std::vector<double> result;result.reserve(splitDimension());
    for(std::size_t mode=0;mode<liftedDimension();++mode) {
        if(!std::isfinite(values[mode].real())||!std::isfinite(values[mode].imag()))
            throw std::invalid_argument("Eigenfunction values must be finite");
        result.push_back(values[mode].real());
        if(!isRealMode(mode))result.push_back(values[mode].imag());
    }
    return result;
}

std::vector<double> KoopmanModel::generatorSplit() const {
    const std::size_t dimension = splitDimension();
    std::vector<double> result(dimension * dimension, 0.0); std::size_t row = 0;
    for (const auto value : continuousEigenvalues()) {
        result[row * dimension + row] = value.real();
        if (realMode(value)) { ++row; continue; }
        result[row * dimension + row + 1] = -value.imag();
        result[(row + 1) * dimension + row] = value.imag();
        result[(row + 1) * dimension + row + 1] = value.real(); row += 2;
    }
    return result;
}

std::vector<std::complex<double>> KoopmanModel::exactControlMap(const double* state) const {
    if(!system_)throw std::logic_error("State-Jacobian control map needs a physical-system adapter");
    const auto original = system_->inputMatrix(state); const std::size_t inputs = inputDimension();
    if(scattered_){
        const auto gradients=scattered_->gradient(state);
        std::vector<std::complex<double>> result(liftedDimension()*inputs);
        for(std::size_t mode=0;mode<liftedDimension();++mode)
            for(std::size_t input=0;input<inputs;++input)
                for(std::size_t axis=0;axis<stateDimension();++axis)
                    result[mode*inputs+input]+=gradients[mode*stateDimension()+axis]
                        *original[axis*inputs+input];
        return result;
    }
    if (!coordinateObservable_)
        throw std::logic_error(
            "an exact state-Jacobian control map is unavailable for causal-history coordinates");
    const auto coordinates=coordinateObservable_->evaluateReal(state);
    const auto coordinateJacobian=coordinateObservable_->jacobian(state);
    std::vector<std::complex<double>> result(liftedDimension() * inputs);
    for (std::size_t mode = 0; mode < liftedDimension(); ++mode) {
        const auto gradient = polynomial_->gradient(mode,coordinates.data());
        for (std::size_t input = 0; input < inputs; ++input)
            for (std::size_t coordinate=0;coordinate<coordinateDimension();++coordinate)
                for (std::size_t axis = 0; axis < stateDimension(); ++axis)
                    result[mode * inputs + input] += gradient[coordinate]
                        * coordinateJacobian[coordinate*stateDimension()+axis]
                        * original[axis * inputs + input];
    }
    return result;
}

std::vector<std::complex<double>> KoopmanModel::controlMap(
    const double* state, const std::string& representation) const {
    if (representation == "exact") return exactControlMap(state);
    const auto found = controlMaps_.find(representation);
    if (found == controlMaps_.end()) throw std::out_of_range("unknown KROM control map");
    return evaluateStored(found->second, evaluate(state));
}

std::vector<std::complex<double>> KoopmanModel::controlMapCoordinates(
    const double* coordinates, const std::string& representation) const {
    return controlMapFromLifted(evaluateCoordinates(coordinates),representation);
}

std::vector<std::complex<double>> KoopmanModel::controlMapFromLifted(
    const std::vector<std::complex<double>>& values, const std::string& representation) const {
    if(values.size()!=eigenfunctionCount())throw std::invalid_argument("Input map requires mode-aligned values");
    if (representation == "exact")
        throw std::logic_error(
            "coordinate-supplied control evaluation requires an approximate control map");
    const auto found = controlMaps_.find(representation);
    if (found == controlMaps_.end()) throw std::out_of_range("unknown KROM control map");
    return evaluateStored(found->second, values);
}

std::vector<double> KoopmanModel::controlMapSplit(
    const double* state, const std::string& representation) const {
    return realControlMap(controlMap(state,representation));
}

std::vector<double> KoopmanModel::realControlMap(const std::vector<std::complex<double>>& complex) const {
    const std::size_t inputs=inputDimension();
    if(complex.size()!=eigenfunctionCount()*inputs)throw std::invalid_argument("Input map dimensions differ from KROM");
    for(auto z:complex)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))
        throw std::invalid_argument("Input map values must be finite");
    std::vector<double> result(splitDimension() * inputs); std::size_t split = 0;
    for (std::size_t mode = 0; mode < liftedDimension(); ++mode) {
        for (std::size_t input = 0; input < inputs; ++input)
            result[split * inputs + input] = complex[mode * inputs + input].real();
        ++split;
        if (!isRealMode(mode)) {
            for (std::size_t input = 0; input < inputs; ++input)
                result[split * inputs + input] = complex[mode * inputs + input].imag();
            ++split;
        }
    }
    return result;
}

void KoopmanModel::setControlMap(std::string name, ApproximateControlMap representation) {
    if (name.empty() || name == "exact") throw std::invalid_argument("invalid approximate map name");
    std::visit([&](const auto& value) { validate(value, liftedDimension(), inputDimension()); },
               representation);
    const auto identity = koopmanModelIdentity(*this);
    std::visit([&](auto& value) {
        if (!value.modelIdentity.empty() && value.modelIdentity != identity)
            throw std::invalid_argument("Control map belongs to different ordered eigenfunctions/coordinates/physics");
        value.modelIdentity = identity;
    }, representation);
    controlMaps_[std::move(name)] = std::move(representation);
}
bool KoopmanModel::hasControlMap(const std::string& name) const {
    return (name == "exact" && system_ && !requiresExternalCoordinates())
        || controlMaps_.find(name) != controlMaps_.end();
}
const ApproximateControlMap& KoopmanModel::storedControlMap(const std::string& name) const {
    return controlMaps_.at(name);
}
const std::map<std::string, ApproximateControlMap>& KoopmanModel::controlMaps() const noexcept {
    return controlMaps_;
}

void KoopmanModel::setReconstruction(std::string observable,
                                     ReconstructionResult reconstruction) {
    if (observable.empty() || !reconstruction.outputs
        || reconstruction.modes!=liftedDimension()
        || reconstruction.matrix.size()!=reconstruction.outputs*liftedDimension())
        throw std::invalid_argument("observable reconstruction dimensions do not match KROM");
    reconstructions_[std::move(observable)]=std::move(reconstruction);
}
bool KoopmanModel::hasReconstruction(const std::string& observable) const {
    return reconstructions_.count(observable)!=0;
}
const ReconstructionResult& KoopmanModel::reconstruction(
    const std::string& observable) const { return reconstructions_.at(observable); }
std::vector<std::complex<double>> KoopmanModel::reconstruct(
    const std::string& observable, const double* state) const {
    const auto& mapping=reconstruction(observable);
    const auto lifted=evaluate(state);
    if(mapping.modes!=lifted.size()
        || mapping.matrix.size()!=mapping.outputs*mapping.modes)
        throw std::runtime_error("stored observable reconstruction does not match KROM");
    std::vector<std::complex<double>> result(mapping.outputs);
    for(std::size_t output=0;output<mapping.outputs;++output)
        for(std::size_t mode=0;mode<mapping.modes;++mode)
            result[output]+=lifted[mode]*mapping.matrix[output*mapping.modes+mode];
    return result;
}
const std::map<std::string, ReconstructionResult>&
KoopmanModel::reconstructions() const noexcept { return reconstructions_; }

KoopmanModel KoopmanModel::extract(const std::vector<std::size_t>& indices) const {
    if (indices.empty()) throw std::invalid_argument("KROM extraction requires a mode");
    std::set<std::size_t> unique(indices.begin(), indices.end());
    if (unique.size() != indices.size() || *unique.rbegin() >= liftedDimension())
        throw std::out_of_range("KROM extraction indices are invalid");
    KoopmanModel result=scattered_
        ? KoopmanModel(scattered_->extract(indices)->data(),systemConfiguration_)
        : KoopmanModel(polynomial_->extract(indices),inputDimension(),fitDiagnostics_.extract(indices));
    if(!scattered_&&system_)result.attachPhysicalSystem(systemConfiguration_);
    result.datasetReferences_=datasetReferences_;
    result.bindingEvidenceJson_=bindingEvidenceJson_;
    const std::size_t inputs = inputDimension(), modes = indices.size();
    for(const auto& entry:discreteInputMaps_) {
        auto map=entry.second;map.matrix.clear();map.eigenfunctionCount=modes;map.modelIdentity.clear();
        for(auto j:indices)map.matrix.insert(map.matrix.end(),entry.second.matrix.begin()+j*inputs,
            entry.second.matrix.begin()+(j+1)*inputs);
        result.setDiscreteInputMap(entry.first,std::move(map));
    }
    for (const auto& item : controlMaps_) {
        if (const auto* constantSource = std::get_if<ConstantControlMap>(&item.second)) {
            ConstantControlMap map{{}, modes, inputs, constantSource->provenance + "; independent row subset"};
            map.sampleTime=constantSource->sampleTime;
            for (const auto row : indices)
                map.matrix.insert(map.matrix.end(), constantSource->matrix.begin() + row * inputs,
                                  constantSource->matrix.begin() + (row + 1) * inputs);
            if(!constantSource->discreteMatrix.empty())
                for(const auto row:indices)
                    map.discreteMatrix.insert(map.discreteMatrix.end(),
                        constantSource->discreteMatrix.begin()+row*inputs,
                        constantSource->discreteMatrix.begin()+(row+1)*inputs);
            result.setControlMap(item.first, std::move(map));
        } else {
            const auto& bilinearSource = std::get<BilinearControlMap>(item.second);
            BilinearControlMap map{{}, {}, {}, modes, inputs,
                bilinearSource.provenance + "; independent principal restriction"};
            map.operators.resize(inputs * modes * modes); map.offsets.resize(modes * inputs);
            if (!bilinearSource.conjugateOperators.empty()) map.conjugateOperators.resize(map.operators.size());
            for (std::size_t input = 0; input < inputs; ++input)
                for (std::size_t row = 0; row < modes; ++row) {
                    map.offsets[row * inputs + input] = bilinearSource.offsets[indices[row] * inputs + input];
                    for (std::size_t column = 0; column < modes; ++column) {
                        const auto target = (input * modes + row) * modes + column;
                        const auto origin = (input * bilinearSource.modes + indices[row])
                            * bilinearSource.modes + indices[column];
                        map.operators[target] = bilinearSource.operators[origin];
                        if (!map.conjugateOperators.empty())
                            map.conjugateOperators[target] = bilinearSource.conjugateOperators[origin];
                    }
                }
            result.setControlMap(item.first, std::move(map));
        }
    }
    for (const auto& item : reconstructions_) {
        ReconstructionResult projection=item.second;
        projection.modes=modes; projection.matrix.assign(projection.outputs*modes,{});
        for (std::size_t output=0;output<projection.outputs;++output)
            for (std::size_t column=0;column<modes;++column)
                projection.matrix[output*modes+column]=
                    item.second.matrix[output*liftedDimension()+indices[column]];
        result.setReconstruction(item.first,std::move(projection));
    }
    return result;
}

const std::vector<KoopmanDatasetReference>& KoopmanModel::datasetReferences() const noexcept {
    return datasetReferences_;
}
void KoopmanModel::addDatasetReference(KoopmanDatasetReference reference) {
    mergeKoopmanDatasetReference(datasetReferences_,std::move(reference));
}

void KoopmanModel::setBindingEvidenceJson(std::string evidence) {
    if (!std::holds_alternative<json::Value::Object>(json::parse(evidence).storage()))
        throw std::invalid_argument("Binding evidence must be a JSON object");
    bindingEvidenceJson_ = std::move(evidence);
}

}  // namespace dynsyslab::core
