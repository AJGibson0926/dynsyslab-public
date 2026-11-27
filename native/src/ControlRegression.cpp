#include <dynsyslab/core/ControlRegression.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include "FlatMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using CMatrix=detail::FlatMatrix<std::complex<double>>;
using RMatrix=detail::FlatMatrix<double>;

CMatrix product(const CMatrix& left,const CMatrix& right) {
    return {numerics::multiply(left.values,left.rows,left.columns,
        right.values,right.columns),left.rows,right.columns};
}

RMatrix product(const RMatrix& left,const RMatrix& right) {
    return {numerics::multiply(left.values,left.rows,left.columns,
        right.values,right.columns),left.rows,right.columns};
}

// What a regression reads from a record: states (sample x state) and inputs (sample x input, or
// interval x input when the record holds inputs on intervals).
struct Controlled {
    const std::vector<double> &time, &states, &inputs;
    std::size_t samples{}, stateDimension{}, inputDimension{};
    bool onIntervals{};
};

Controlled controlled(const TimeSeriesRecord& r) {
    static const std::vector<double> none;
    return {r.series.timeVector(), r.series.valueVector(), r.inputs ? r.inputs->valueVector() : none,
            r.series.sampleCount(), r.series.componentCount(), r.inputs ? r.inputs->componentCount() : 0,
            r.inputsOnIntervals};
}

std::vector<Controlled> controlled(const std::vector<TimeSeriesRecord>& records) {
    std::vector<Controlled> result;
    for (const auto& record : records) result.push_back(controlled(record));
    return result;
}

void validateStates(const KoopmanModel& model, const std::vector<double>& states,
                    std::size_t samples) {
    if (samples < 3 || states.size() != samples * model.stateDimension())
        throw std::invalid_argument("control regression state array is inconsistent");
}

CMatrix lifts(const KoopmanModel& model, const std::vector<double>& states,
              std::size_t samples) {
    validateStates(model, states, samples);
    CMatrix result(samples,model.liftedDimension());
    for (std::size_t row = 0; row < samples; ++row) {
        const auto value = model.evaluate(states.data() + row * model.stateDimension());
        for (std::size_t column = 0; column < value.size(); ++column)
            result(row,column)=value[column];
    }
    return result;
}

CMatrix exactMaps(const KoopmanModel& model, const std::vector<double>& states,
                  std::size_t samples) {
    const std::size_t modes = model.liftedDimension(), inputs = model.inputDimension();
    CMatrix result(samples, modes * inputs);
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const auto value = model.exactControlMap(states.data() + sample * model.stateDimension());
        for (std::size_t index = 0; index < value.size(); ++index)
            result(sample,index)=value[index];
    }
    return result;
}

CMatrix ridgeSolve(const CMatrix& design, const CMatrix& target, double ridge) {
    if (ridge < 0.0 || !std::isfinite(ridge)) throw std::invalid_argument("ridge must be nonnegative");
    numerics::LeastSquaresOptions options;
    options.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    options.ridge=ridge;
    return {numerics::solveLeastSquares(design.values,design.rows,design.columns,
        target.values,target.columns,options),design.columns,target.columns};
}

ControlRegressionDiagnostics diagnostics(const char* method, const CMatrix& design,
                                         const CMatrix& target, const CMatrix& predicted) {
    if(predicted.rows!=target.rows||predicted.columns!=target.columns)
        throw std::invalid_argument("regression diagnostic shapes do not match");
    const auto singular=numerics::singularValues(
        design.values,design.rows,design.columns);
    const double maximum=singular.empty()?0.0:singular.front();
    std::size_t rank = 0; double minimum = std::numeric_limits<double>::infinity();
    for (const double value:singular)
        if (value > std::max(1.0, maximum) * 1e-12) {
            ++rank; minimum = std::min(minimum,value);
        }
    double residualSquared=0.0,targetSquared=0.0;
    for(std::size_t index=0;index<target.values.size();++index) {
        residualSquared+=std::norm(predicted.values[index]-target.values[index]);
        targetSquared+=std::norm(target.values[index]);
    }
    return {method,design.rows,design.columns,rank,
        minimum < std::numeric_limits<double>::infinity() ? maximum / minimum
                                                          : std::numeric_limits<double>::infinity(),
        std::sqrt(residualSquared/std::max<std::size_t>(1,target.values.size())),
        std::sqrt(residualSquared)/std::max(std::sqrt(targetSquared),
            std::numeric_limits<double>::min())};
}

BilinearControlFit bilinearFit(const KoopmanModel& model, const CMatrix& lift,
                               const CMatrix& target, double ridge,
                               bool conjugateComplete, const char* method) {
    const std::size_t modes = model.liftedDimension(), inputs = model.inputDimension();
    const std::size_t features = modes * (conjugateComplete ? 2 : 1) + 1;
    CMatrix design(lift.rows,features);
    for(std::size_t row=0;row<lift.rows;++row)
        for(std::size_t mode=0;mode<modes;++mode)design(row,mode)=lift(row,mode);
    std::size_t cursor = modes;
    if (conjugateComplete) {
        for(std::size_t row=0;row<lift.rows;++row)
            for(std::size_t mode=0;mode<modes;++mode)
                design(row,cursor+mode)=std::conj(lift(row,mode));
        cursor += modes;
    }
    for(std::size_t row=0;row<lift.rows;++row)design(row,cursor)=1.0;
    BilinearControlMap map{{}, {}, {}, modes, inputs, method};
    map.operators.resize(inputs * modes * modes); map.offsets.resize(modes * inputs);
    if (conjugateComplete) map.conjugateOperators.resize(map.operators.size());
    CMatrix predicted(target.rows,target.columns);
    for (std::size_t input = 0; input < inputs; ++input) {
        CMatrix block(target.rows,modes);
        for (std::size_t mode = 0; mode < modes; ++mode)
            for(std::size_t row=0;row<target.rows;++row)
                block(row,mode)=target(row,mode*inputs+input);
        const CMatrix coefficient = ridgeSolve(design, block, ridge);
        const CMatrix prediction=product(design,coefficient);
        for (std::size_t row = 0; row < modes; ++row) {
            map.offsets[row*inputs+input]=coefficient(cursor,row);
            for (std::size_t column = 0; column < modes; ++column) {
                const auto at = (input * modes + row) * modes + column;
                map.operators[at]=coefficient(column,row);
                if (conjugateComplete) map.conjugateOperators[at] =
                    coefficient(modes+column,row);
            }
            for(std::size_t sample=0;sample<target.rows;++sample)
                predicted(sample,row*inputs+input)=prediction(sample,row);
        }
    }
    return {std::move(map), diagnostics(method, design, target, predicted)};
}

void validateControlled(const KoopmanModel& model, const Controlled& data) {
    validateStates(model, data.states, data.samples);
    if (data.stateDimension != model.stateDimension() || data.inputDimension != model.inputDimension()
        || data.time.size() != data.samples
        // Inputs are held per sample or per interval; either count is accepted, as before the
        // record carried the alignment flag (the regressions index both the same way).
        || (data.inputs.size() != data.samples * data.inputDimension
            && data.inputs.size() != (data.samples - 1) * data.inputDimension))
        throw std::invalid_argument("controlled dataset dimensions do not match KROM");
    for (std::size_t index = 1; index < data.time.size(); ++index)
        if (!(data.time[index] > data.time[index - 1]))
            throw std::invalid_argument("controlled dataset time must increase");
}

struct ResidualData { CMatrix lift; CMatrix input; CMatrix residual; };

ResidualData residualData(const KoopmanModel& model, const Controlled& data) {
    validateControlled(model, data); const CMatrix all = lifts(model, data.states, data.samples);
    const std::size_t rows = data.samples - 2, modes = model.liftedDimension();
    CMatrix active(rows, modes), input(rows, data.inputDimension), residual(rows, modes);
    for (std::size_t sample = 1; sample + 1 < data.samples; ++sample) {
        const double dt = data.time[sample + 1] - data.time[sample - 1];
        const std::size_t row=sample-1;
        for (std::size_t control = 0; control < data.inputDimension; ++control)
            input(row,control)=
                data.inputs[std::min(sample, data.samples - 2) * data.inputDimension + control];
        for (std::size_t mode = 0; mode < modes; ++mode) {
            active(row,mode)=all(sample,mode);
            const auto derivative=(all(sample+1,mode)-all(sample-1,mode))/dt;
            residual(row,mode)=derivative-all(sample,mode)
                * model.continuousEigenvalues()[mode];
        }
    }
    return {std::move(active), std::move(input), std::move(residual)};
}

ResidualData residualData(const KoopmanModel& model,
                          const std::vector<Controlled>& datasets) {
    if(datasets.empty())throw std::invalid_argument("control regression ensemble is empty");
    std::size_t rows=0;for(const auto& data:datasets)rows+=data.samples-2;
    ResidualData result{CMatrix(rows,model.liftedDimension()),
        CMatrix(rows,model.inputDimension()),CMatrix(rows,model.liftedDimension())};
    std::size_t offset=0;
    for(const auto& data:datasets){const auto block=residualData(model,data);
        for(std::size_t row=0;row<block.lift.rows;++row) {
            for(std::size_t column=0;column<block.lift.columns;++column)
                result.lift(offset+row,column)=block.lift(row,column);
            for(std::size_t column=0;column<block.input.columns;++column)
                result.input(offset+row,column)=block.input(row,column);
            for(std::size_t column=0;column<block.residual.columns;++column)
                result.residual(offset+row,column)=block.residual(row,column);
        }
        offset+=block.lift.rows;}
    return result;
}

double uniformStep(const Controlled& data) {
    const double step=data.time[1]-data.time[0];
    if (!(step>0.0) || !std::isfinite(step))
        throw std::invalid_argument("multistep control regression requires positive sample time");
    for(std::size_t index=2;index<data.time.size();++index)
        if(std::abs((data.time[index]-data.time[index-1])-step)
            >1.0e-10*std::max(1.0,std::abs(step)))
            throw std::invalid_argument("multistep control regression requires uniform sample time");
    return step;
}

bool realMode(std::complex<double> value) {
    return std::abs(value.imag())<=1.0e-12;
}

std::pair<RMatrix,RMatrix> heldKoopmanPair(const KoopmanModel& model,double step) {
    const std::size_t dimension=model.splitDimension();
    RMatrix transition(dimension,dimension),integral(dimension,dimension);
    std::size_t row=0;
    for(const auto value:model.continuousEigenvalues()){
        const double a=value.real(),b=value.imag();
        if(realMode(value)){
            const double discrete=std::exp(a*step);
            transition(row,row)=discrete;
            integral(row,row)=std::abs(a)>1.0e-14?(discrete-1.0)/a:step;
            ++row;continue;
        }
        const double scale=std::exp(a*step),angle=b*step;
        const double cosine=scale*std::cos(angle),sine=scale*std::sin(angle);
        transition(row,row)=cosine;transition(row,row+1)=-sine;
        transition(row+1,row)=sine;transition(row+1,row+1)=cosine;
        const double denominator=a*a+b*b;
        double cosineIntegral=step,sineIntegral=0.0;
        if(denominator>1.0e-28){
            cosineIntegral=(scale*(a*std::cos(angle)+b*std::sin(angle))-a)/denominator;
            sineIntegral=(scale*(a*std::sin(angle)-b*std::cos(angle))+b)/denominator;
        }
        integral(row,row)=cosineIntegral;integral(row,row+1)=-sineIntegral;
        integral(row+1,row)=sineIntegral;integral(row+1,row+1)=cosineIntegral;
        row+=2;
    }
    return {std::move(transition),std::move(integral)};
}

RMatrix splitReconstruction(const KoopmanModel& model,
                            const ReconstructionResult& reconstruction) {
    if(reconstruction.modes!=model.liftedDimension()
        ||reconstruction.matrix.size()!=reconstruction.outputs*reconstruction.modes)
        throw std::invalid_argument("multistep reconstruction does not match KROM");
    RMatrix result(reconstruction.outputs,model.splitDimension());
    std::size_t split=0;
    for(std::size_t mode=0;mode<model.liftedDimension();++mode){
        for(std::size_t output=0;output<reconstruction.outputs;++output){
            const auto coefficient=reconstruction.matrix[output*reconstruction.modes+mode];
            result(output,split)=realMode(model.continuousEigenvalues()[mode])
                ?coefficient.real():2.0*coefficient.real();
            if(!realMode(model.continuousEigenvalues()[mode]))
                result(output,split+1)=-2.0*coefficient.imag();
        }
        split+=realMode(model.continuousEigenvalues()[mode])?1:2;
    }
    return result;
}

double intervalInput(const Controlled& data,std::size_t interval,
                     std::size_t input) {
    return data.inputs[interval*data.inputDimension+input];
}

}  // namespace

ConstantControlFit fitConstantExactControlMap(const KoopmanModel& model,
                                               const std::vector<double>& states,
                                               std::size_t samples) {
    validateStates(model, states, samples); const CMatrix target = exactMaps(model, states, samples);
    CMatrix design(samples,1);std::fill(design.values.begin(),design.values.end(),1.0);
    const CMatrix coefficient=ridgeSolve(design,target,0.0);
    ConstantControlMap map{{}, model.liftedDimension(), model.inputDimension(), "pointwise exact-map mean"};
    map.matrix=coefficient.values;
    return {std::move(map), diagnostics("pointwise-constant-exact-map/v1",
        design,target,product(design,coefficient))};
}

BilinearControlFit fitBilinearExactControlMap(const KoopmanModel& model,
                                               const std::vector<double>& states,
                                               std::size_t samples, double ridge,
                                               bool conjugateComplete) {
    const auto lift = lifts(model, states, samples); const auto target = exactMaps(model, states, samples);
    return bilinearFit(model, lift, target, ridge, conjugateComplete,
        conjugateComplete ? "pointwise-widely-linear-exact-map/v1" : "pointwise-bilinear-exact-map/v1");
}

ConstantControlFit fitConstantControlResidual(const KoopmanModel& model,
                                               const TimeSeriesRecord& data, double ridge) {
    return fitConstantControlResidual(model,std::vector<TimeSeriesRecord>{data},ridge);
}

ConstantControlFit fitConstantControlResidual(const KoopmanModel& model,
    const std::vector<TimeSeriesRecord>& data,double ridge) {
    const auto source = residualData(model, controlled(data));
    const CMatrix coefficient = ridgeSolve(source.input, source.residual, ridge);
    ConstantControlMap map{{}, model.liftedDimension(), model.inputDimension(), "continuous controlled residual"};
    map.matrix.resize(model.liftedDimension() * model.inputDimension());
    for (std::size_t mode = 0; mode < model.liftedDimension(); ++mode)
        for (std::size_t input = 0; input < model.inputDimension(); ++input)
            map.matrix[mode * model.inputDimension() + input] =
                coefficient(input,mode);
    return {std::move(map), diagnostics("continuous-constant-lifted-residual/v1", source.input,
        source.residual,product(source.input,coefficient))};
}

ConstantControlFit fitConstantControlMultistepOutput(
    const KoopmanModel& model,const TimeSeriesRecord& data,
    const Observable& observable,const ReconstructionResult& reconstruction,
    double ridge) {
    return fitConstantControlMultistepOutput(model,
        std::vector<TimeSeriesRecord>{data},observable,reconstruction,ridge);
}

ConstantControlFit fitConstantControlMultistepOutput(
    const KoopmanModel& model,const std::vector<TimeSeriesRecord>& records,
    const Observable& observable,const ReconstructionResult& reconstruction,
    double ridge) {
    const auto datasets = controlled(records);
    if(datasets.empty())throw std::invalid_argument("multistep regression ensemble is empty");
    for(const auto& data:datasets)validateControlled(model,data);
    if(observable.stateDimension()!=model.stateDimension()
        ||observable.outputDimension()!=reconstruction.outputs)
        throw std::invalid_argument("multistep observable and reconstruction do not match");
    const double step=uniformStep(datasets.front());
    for(const auto& data:datasets)if(std::abs(uniformStep(data)-step)
        >1.0e-10*std::max(1.0,std::abs(step)))
        throw std::invalid_argument("multistep trajectories must share one sample time");
    const auto dynamics=heldKoopmanPair(model,step);
    const RMatrix C=splitReconstruction(model,reconstruction);
    const std::size_t lifted=model.splitDimension(),inputs=model.inputDimension();
    std::size_t rows=0;for(const auto& data:datasets)
        rows+=(data.samples-1)*reconstruction.outputs;
    CMatrix design(rows,lifted*inputs),target(rows,1);
    std::size_t outputRow=0;
    for(const auto& data:datasets){
        const auto observed=observable.evaluateRealBatch(data.states,data.samples);
        std::vector<RMatrix> response;
        response.reserve(inputs);
        for(std::size_t input=0;input<inputs;++input)
            response.emplace_back(reconstruction.outputs,lifted);
        auto freeState=model.evaluateSplit(data.states.data());
        for(std::size_t sample=1;sample<data.samples;++sample){
            for(std::size_t input=0;input<inputs;++input) {
                response[input]=product(response[input],dynamics.first);
                const double held=intervalInput(data,sample-1,input);
                for(std::size_t index=0;index<response[input].values.size();++index)
                    response[input].values[index]+=C.values[index]*held;
            }
            std::vector<double> nextState(lifted);
            for(std::size_t row=0;row<lifted;++row)
                for(std::size_t column=0;column<lifted;++column)
                    nextState[row]+=dynamics.first(row,column)*freeState[column];
            freeState=std::move(nextState);
            std::vector<double> freeOutput(reconstruction.outputs);
            for(std::size_t output=0;output<reconstruction.outputs;++output)
                for(std::size_t state=0;state<lifted;++state)
                    freeOutput[output]+=C(output,state)*freeState[state];
            for(std::size_t output=0;output<reconstruction.outputs;++output,++outputRow){
                for(std::size_t input=0;input<inputs;++input)
                    for(std::size_t state=0;state<lifted;++state)
                        design(outputRow,input*lifted+state)=response[input](output,state);
                target(outputRow,0)=observed[sample*reconstruction.outputs+output]
                    -freeOutput[output];
            }
        }
    }
    const CMatrix coefficient=ridgeSolve(design,target,ridge);
    RMatrix discrete(lifted,inputs);
    for(std::size_t input=0;input<inputs;++input)
        for(std::size_t state=0;state<lifted;++state)
            discrete(state,input)=coefficient(input*lifted+state,0).real();
    numerics::LeastSquaresOptions solveOptions;
    solveOptions.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    const RMatrix continuous{numerics::solveLeastSquares(
        dynamics.second.values,lifted,lifted,discrete.values,inputs,solveOptions),
        lifted,inputs};
    ConstantControlMap map{{},model.liftedDimension(),inputs,
        "Korda-Mezic multistep observable-output regression",{},step};
    std::size_t split=0;
    for(const auto value:model.continuousEigenvalues()){
        for(std::size_t input=0;input<inputs;++input){
            map.matrix.push_back({continuous(split,input),
                realMode(value)?0.0:continuous(split+1,input)});
            map.discreteMatrix.push_back({discrete(split,input),
                realMode(value)?0.0:discrete(split+1,input)});
        }
        split+=realMode(value)?1:2;
    }
    return {std::move(map),diagnostics("korda-mezic-multistep-output/v1",
        design,target,product(design,coefficient))};
}

BilinearControlFit fitBilinearControlResidual(const KoopmanModel& model,
                                               const TimeSeriesRecord& data, double ridge,
                                               bool conjugateComplete) {
    return fitBilinearControlResidual(model,std::vector<TimeSeriesRecord>{data},
                                      ridge,conjugateComplete);
}

BilinearControlFit fitBilinearControlResidual(const KoopmanModel& model,
    const std::vector<TimeSeriesRecord>& data,double ridge,bool conjugateComplete) {
    const auto source = residualData(model, controlled(data)); const std::size_t modes = model.liftedDimension();
    const std::size_t inputs=model.inputDimension();
    CMatrix features(source.lift.rows,
        inputs*(modes*(conjugateComplete?2:1)+1));
    std::size_t column = 0;
    for (std::size_t input = 0; input < inputs; ++input) {
        for (std::size_t mode = 0; mode < modes; ++mode)
        {
            for(std::size_t row=0;row<features.rows;++row)
                features(row,column)=source.input(row,input)*source.lift(row,mode);
            ++column;
        }
        if (conjugateComplete) for (std::size_t mode = 0; mode < modes; ++mode)
        {
            for(std::size_t row=0;row<features.rows;++row)
                features(row,column)=source.input(row,input)*std::conj(source.lift(row,mode));
            ++column;
        }
        for(std::size_t row=0;row<features.rows;++row)
            features(row,column)=source.input(row,input);
        ++column;
    }
    const CMatrix coefficient = ridgeSolve(features, source.residual, ridge);
    BilinearControlMap map{{}, {}, {}, modes, inputs, "continuous controlled residual"};
    map.operators.resize(inputs * modes * modes); map.offsets.resize(modes * inputs);
    if (conjugateComplete) map.conjugateOperators.resize(map.operators.size());
    column = 0;
    for (std::size_t input = 0; input < inputs; ++input) {
        for (std::size_t predictor = 0; predictor < modes; ++predictor, ++column)
            for (std::size_t row = 0; row < modes; ++row)
                map.operators[(input * modes + row) * modes + predictor] =
                    coefficient(column,row);
        if (conjugateComplete) for (std::size_t predictor = 0; predictor < modes; ++predictor, ++column)
            for (std::size_t row = 0; row < modes; ++row)
                map.conjugateOperators[(input * modes + row) * modes + predictor] =
                    coefficient(column,row);
        for (std::size_t row = 0; row < modes; ++row)
                map.offsets[row * inputs + input] =
                coefficient(column,row);
        ++column;
    }
    return {std::move(map), diagnostics("continuous-bilinear-lifted-residual/v1", features,
        source.residual,product(features,coefficient))};
}

ControlMapErrorGrid evaluateControlMapError2D(
    const KoopmanModel& model, const std::string& approximation,
    double axis0Minimum, double axis0Maximum, double axis1Minimum, double axis1Maximum,
    std::size_t columns, std::size_t rows, const std::vector<double>& fixedState, double relativeFloor) {
    if (model.stateDimension() < 2 || columns < 2 || rows < 2
        || !(axis0Minimum < axis0Maximum) || !(axis1Minimum < axis1Maximum))
        throw std::invalid_argument("invalid two-dimensional control-map error grid");
    std::vector<double> state = fixedState;
    if (state.empty()) state = createSystem(model.systemConfiguration())->initialState();
    if (state.size() != model.stateDimension()) throw std::invalid_argument("fixed grid state dimension mismatch");
    ControlMapErrorGrid result; result.columns = columns; result.rows = rows;
    result.axis0.resize(columns); result.axis1.resize(rows);
    result.absoluteError.resize(rows * columns); result.relativeError.resize(rows * columns);
    for (std::size_t column = 0; column < columns; ++column)
        result.axis0[column] = axis0Minimum + (axis0Maximum-axis0Minimum)*column/(columns-1);
    for (std::size_t row = 0; row < rows; ++row)
        result.axis1[row] = axis1Minimum + (axis1Maximum-axis1Minimum)*row/(rows-1);
    double absoluteSum = 0.0, relativeSum = 0.0;
    for (std::size_t row = 0; row < rows; ++row) for (std::size_t column = 0; column < columns; ++column) {
        state[0] = result.axis0[column]; state[1] = result.axis1[row];
        const auto exact = model.controlMap(state.data(), "exact");
        const auto fitted = model.controlMap(state.data(), approximation);
        double exactNorm = 0.0, errorNorm = 0.0;
        for (std::size_t index = 0; index < exact.size(); ++index) {
            exactNorm += std::norm(exact[index]); errorNorm += std::norm(fitted[index] - exact[index]);
        }
        const double absolute = std::sqrt(errorNorm);
        const double relative = absolute / std::max({std::sqrt(exactNorm), relativeFloor,
                                                     std::numeric_limits<double>::min()});
        const auto at = row * columns + column; result.absoluteError[at] = absolute;
        result.relativeError[at] = relative; absoluteSum += absolute; relativeSum += relative;
        result.absoluteMaximum = std::max(result.absoluteMaximum, absolute);
        result.relativeMaximum = std::max(result.relativeMaximum, relative);
    }
    result.absoluteMean = absoluteSum / (rows * columns);
    result.relativeMean = relativeSum / (rows * columns);
    const auto percentile = [](std::vector<double> values, double q) {
        std::sort(values.begin(), values.end());
        const auto rank = static_cast<std::size_t>(std::ceil(q * double(values.size())));
        return values[std::min(values.size() - 1, rank ? rank - 1 : 0)];
    };
    result.absoluteMedian = percentile(result.absoluteError, 0.5);
    result.absoluteP95 = percentile(result.absoluteError, 0.95);
    result.relativeMedian = percentile(result.relativeError, 0.5);
    result.relativeP95 = percentile(result.relativeError, 0.95);
    return result;
}

}  // namespace dynsyslab::core
