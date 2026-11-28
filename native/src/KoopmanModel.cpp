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

// A K-ROM need not have an input map B: it is fitted later (input_map, control_regression, ...).
std::string missingControlMap(const std::string& name) {
    return "this K-ROM has no input map B named \"" + name + "\"; fit one first (input_map, control_regression, "
           "history_control_regression or DMDc/EDMDc) and use the krom that step publishes";
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

// An interpolation's own samples: its training values at its base points, in order.
EigenfunctionSamples interpolationSamples(const InterpolationData& data) {
    EigenfunctionSamples samples;
    samples.modes = data.outputs();
    samples.values = data.trainingValues;
    for (std::size_t k = 0; k < data.samples; ++k) { samples.sampleIndices.push_back(k); samples.time.push_back(double(k)); }
    return samples;
}

}  // namespace

KoopmanModel::KoopmanModel(InterpolationData eigenfunctions, Spectrum spectrum,
                           SystemConfiguration system)
    : KoopmanModel(std::move(spectrum), interpolationSamples(eigenfunctions), DataOrigin{},
                   {std::make_shared<const FunctionInterpolation>(eigenfunctions)}) {
    bindSystem(std::move(system));
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
void KoopmanModel::setAdditiveInputMatrix(std::vector<double> matrix) {
    const std::size_t n = stateDimension();
    if (!matrix.empty() && (n == 0 || matrix.size() % n))
        throw std::invalid_argument("additive input matrix must be state dimension x inputs");
    for (const double v : matrix)
        if (!std::isfinite(v)) throw std::invalid_argument("additive input matrix must be finite");
    additiveInput_ = std::move(matrix);
}
std::size_t KoopmanModel::inputDimension() const {
    if(!additiveInput_.empty())return additiveInput_.size()/stateDimension();
    if(!system_)return suppliedInputDimension_;
    const auto initial = system_->initialState();
    const auto matrix = system_->inputMatrix(initial.data());
    const std::size_t states = system_->descriptor().stateDimension;  // the system's, whatever the domain
    if (matrix.empty() || matrix.size() % states)
        throw std::runtime_error("system input matrix has invalid dimensions");
    return matrix.size() / states;
}
const std::vector<std::complex<double>>&
KoopmanModel::continuousEigenvalues() const {
    return spectrum().continuousRates();
}
std::string KoopmanModel::eigenfunctionRepresentation() const {
    if(scattered_)return "scattered-interpolation";
    return polynomial_->monomial()?"analytic-polynomial":"expansion";
}
std::vector<KoopmanEigenfunctionDescriptor>
KoopmanModel::eigenfunctionDescriptors() const {
    std::vector<KoopmanEigenfunctionDescriptor> result;
    for(std::size_t index=0;index<liftedDimension();++index)
        result.push_back({index,spectrum().continuousRate(index).value_or(std::complex<double>{}),
            eigenfunctionRepresentation(),
            scattered_||!polynomial_->dictionary().exactGradient()?"finite-difference":"analytic",
            spectrum().continuousRate(index).has_value()});
    return result;
}
const InterpolationData*
KoopmanModel::scatteredEigenfunctions() const noexcept {
    return scattered_?&scattered_->data():nullptr;
}
const SystemConfiguration& KoopmanModel::systemConfiguration() const {
    if(!system_)throw std::logic_error("Coordinate-only KROM has no physical-system adapter");
    return systemConfiguration_;
}

std::vector<std::complex<double>> KoopmanModel::evaluate(const double* state) const {
    const auto& q=coordinates();
    if(q.isPhysicalState())return scattered_?scattered_->evaluate(state):evaluateCoordinates(state);
    if (!q.stateMap)
        throw std::logic_error(
            "causal-history KROM requires evaluateCoordinates after updating its history projector");
    const auto values=q.stateMap->evaluateReal(state);
    return scattered_?scattered_->evaluate(values.data()):evaluateCoordinates(values.data());
}

const Observable& KoopmanModel::coordinateObservable() const {
    if (!coordinates().stateMap)
        throw std::logic_error(
            "causal-history KROM has no instantaneous state observable");
    return *coordinates().stateMap;
}

bool KoopmanModel::requiresExternalCoordinates() const noexcept {
    return !scattered_ && coordinates().requiresExternalValues();
}

std::vector<std::complex<double>> KoopmanModel::evaluateCoordinates(const double* coordinates) const {
    if(scattered_){
        if(!scattered_->data().coordinateObservable.empty())return scattered_->evaluate(coordinates);
        throw std::logic_error("direct coordinate evaluation requires polynomial eigenfunctions");
    }
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

std::vector<std::complex<double>> KoopmanModel::valuesFromSplit(
    const std::vector<double>& values) const {
    if(values.size()!=splitDimension())
        throw std::invalid_argument("Real lifted value count differs from KROM");
    for(double value:values)if(!std::isfinite(value))
        throw std::invalid_argument("Real lifted values must be finite");
    std::vector<std::complex<double>> result(liftedDimension());
    std::size_t index=0;
    for(std::size_t mode=0;mode<liftedDimension();++mode) {
        const double real=values[index++];
        const double imaginary=isRealMode(mode)?0.0:values[index++];
        result[mode]={real,imaginary};
    }
    return result;
}

std::vector<std::complex<double>> KoopmanModel::exactControlMap(const double* state) const {
    if(!system_&&additiveInput_.empty())throw std::logic_error("State-Jacobian control map needs a physical-system adapter");
    const auto original = additiveInput_.empty() ? system_->inputMatrix(state) : additiveInput_;
    const std::size_t inputs = inputDimension();
    if(scattered_){
        if(!scattered_->data().coordinateObservable.empty())throw std::logic_error(
            "an observable-bound scattered KROM needs a fitted input map");
        const auto gradients=scattered_->gradient(state);
        std::vector<std::complex<double>> result(liftedDimension()*inputs);
        for(std::size_t mode=0;mode<liftedDimension();++mode)
            for(std::size_t input=0;input<inputs;++input)
                for(std::size_t axis=0;axis<stateDimension();++axis)
                    result[mode*inputs+input]+=gradients[mode*stateDimension()+axis]
                        *original[axis*inputs+input];
        return result;
    }
    const auto& map=this->coordinates().stateMap;
    if (!map)
        throw std::logic_error(
            "an exact state-Jacobian control map is unavailable for causal-history coordinates");
    const auto coordinates=map->evaluateReal(state);
    const auto coordinateJacobian=map->jacobian(state);
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

std::vector<std::complex<double>> KoopmanModel::controlMapFromCoordinateSensitivity(
    const double* coordinates, const std::vector<double>& sensitivity) const {
    if (!polynomial_)
        throw std::logic_error("coordinate-sensitivity control maps need polynomial eigenfunctions");
    const std::size_t inputs = inputDimension(), dimension = coordinateDimension();
    if (sensitivity.size() != dimension * inputs)
        throw std::invalid_argument("coordinate sensitivity must be coordinates x inputs");
    std::vector<std::complex<double>> result(liftedDimension() * inputs);
    for (std::size_t mode = 0; mode < liftedDimension(); ++mode) {
        const auto gradient = polynomial_->gradient(mode, coordinates);
        for (std::size_t input = 0; input < inputs; ++input)
            for (std::size_t c = 0; c < dimension; ++c)
                result[mode * inputs + input] += gradient[c] * sensitivity[c * inputs + input];
    }
    return result;
}

std::vector<std::complex<double>> KoopmanModel::controlMap(
    const double* state, const std::string& representation) const {
    if (representation == "exact") return exactControlMap(state);
    const auto found = controlMaps_.find(representation);
    if (found == controlMaps_.end()) throw std::out_of_range(missingControlMap(representation));
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
    if (found == controlMaps_.end()) throw std::out_of_range(missingControlMap(representation));
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
    inherited_.controlMaps.erase(name);
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
    inherited_.readouts.erase(observable);
    spectrum_->setKoopmanModes(std::move(observable),std::move(reconstruction));
}
bool KoopmanModel::hasReconstruction(const std::string& observable) const {
    return spectrum_->koopmanModes(observable)!=nullptr;
}
const ReconstructionResult& KoopmanModel::reconstruction(
    const std::string& observable) const {
    const auto* modes=spectrum_->koopmanModes(observable);
    if(!modes)throw std::out_of_range("this KROM has no reconstruction of "+observable);
    return *modes;
}
std::vector<std::complex<double>> KoopmanModel::reconstruct(
    const std::string& observable, const double* state) const {
    const auto& mapping=reconstruction(observable);
    const auto lifted=evaluate(state);
    if(mapping.modes!=lifted.size()
        || mapping.matrix.size()!=mapping.outputs*mapping.modes)
        throw std::runtime_error("stored observable reconstruction does not match KROM");
    std::vector<std::complex<double>> result(mapping.outputs);
    for(std::size_t output=0;output<mapping.outputs;++output) {
        if(!mapping.offset.empty())result[output]=mapping.offset[output];
        for(std::size_t mode=0;mode<mapping.modes;++mode)
            result[output]+=lifted[mode]*mapping.matrix[output*mapping.modes+mode];
    }
    return result;
}
const std::map<std::string, ReconstructionResult>&
KoopmanModel::reconstructions() const noexcept { return spectrum_->koopmanModes(); }

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
