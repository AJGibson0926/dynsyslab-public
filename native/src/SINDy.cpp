#include <dynsyslab/core/SINDy.hpp>
#include "PolynomialEvaluation.hpp"
#include <dynsyslab/numerics/Nullspace.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {


double power(double value, std::size_t exponent) {
    double result = 1.0;
    for (std::size_t index = 0; index < exponent; ++index) result *= value;
    return result;
}

std::vector<double> featureMatrix(
    const std::vector<double>& states,std::size_t samples,std::size_t dimension,
    const std::vector<std::size_t>& exponents,
    const std::vector<double>& center,const std::vector<double>& scale) {
    const std::size_t features = exponents.size() / dimension;
    std::vector<double> result(samples*features);
    for (std::size_t row = 0; row < samples; ++row)
        for (std::size_t feature = 0; feature < features; ++feature) {
            double value = 1.0;
            for (std::size_t axis = 0; axis < dimension; ++axis)
                value *= power((states[row*dimension+axis] - center[axis]) / scale[axis],
                               exponents[feature * dimension + axis]);
            result[row*features+feature] = value;
        }
    return result;
}

std::vector<double> finiteDifference(
    const std::vector<double>& values,std::size_t samples,std::size_t columns,
    const std::vector<double>& time) {
    std::vector<double> result(values.size());
    const std::size_t last = samples - 1;
    const double firstStep = time[1] - time[0], lastStep = time[time.size()-1] - time[time.size()-2];
    for(std::size_t column=0;column<columns;++column) {
        result[column]=(-3.0*values[column]+4.0*values[columns+column]
            -values[2*columns+column])/(2.0*firstStep);
        result[last*columns+column]=(3.0*values[last*columns+column]
            -4.0*values[(last-1)*columns+column]
            +values[(last-2)*columns+column])/(2.0*lastStep);
    }
    for (std::size_t row = 1; row < last; ++row) {
        const double backward = time[row] - time[row-1];
        const double forward = time[row+1] - time[row];
        for(std::size_t column=0;column<columns;++column)
            result[row*columns+column]=
                (-forward/(backward*(backward+forward)))*values[(row-1)*columns+column]
                +((forward-backward)/(backward*forward))*values[row*columns+column]
                +(backward/(forward*(backward+forward)))*values[(row+1)*columns+column];
    }
    return result;
}

std::vector<double> featureDerivative(
    const std::vector<double>& states,const std::vector<double>& velocities,
    std::size_t samples,std::size_t dimension,
    const std::vector<std::size_t>& exponents,
    const std::vector<double>& center,const std::vector<double>& scale) {
    const std::size_t features = exponents.size() / dimension;
    std::vector<double> result(samples*features);
    for (std::size_t row = 0; row < samples; ++row)
        for (std::size_t feature = 0; feature < features; ++feature)
            for (std::size_t differentiated = 0; differentiated < dimension; ++differentiated) {
                const std::size_t count = exponents[feature * dimension + differentiated];
                if (!count) continue;
                double derivative = static_cast<double>(count) / scale[differentiated];
                for (std::size_t axis = 0; axis < dimension; ++axis) {
                    std::size_t exponent = exponents[feature * dimension + axis];
                    if (axis == differentiated) --exponent;
                    derivative *= power((states[row*dimension+axis] - center[axis])
                                        / scale[axis], exponent);
                }
                result[row*features+feature] += derivative
                    * velocities[row*dimension+differentiated];
            }
    return result;
}

std::vector<std::complex<double>> sparseNullVector(
    const std::vector<std::complex<double>>& source,std::size_t rows,
    std::size_t columns,double threshold,std::size_t maximumIterations) {
    std::vector<double> scales(columns);
    for (std::size_t column = 0; column < columns; ++column) {
        for(std::size_t row=0;row<rows;++row)
            scales[column]+=std::norm(source[row*columns+column]);
        scales[column]=std::sqrt(scales[column]);
        if (scales[column] <= 1e-15) scales[column] = 1.0;
    }
    std::vector<std::size_t> active(columns);
    std::iota(active.begin(), active.end(), std::size_t{0});
    std::vector<std::complex<double>> vector;
    auto solve = [&]() {
        std::vector<std::complex<double>> block(rows*active.size());
        for(std::size_t row=0;row<rows;++row)
            for(std::size_t index=0;index<active.size();++index)
                block[row*active.size()+index]=source[row*columns+active[index]]
                    /scales[active[index]];
        return numerics::smallestRightSingularVector(
            block,rows,active.size());
    };
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        vector = solve(); double largest=0.0;
        for(const auto value:vector)largest=std::max(largest,std::abs(value));
        const double cutoff = threshold * largest;
        std::vector<std::size_t> retained;
        for (std::size_t index = 0; index < active.size(); ++index)
            if (std::abs(vector[index]) >= cutoff) retained.push_back(active[index]);
        if (retained.empty()) {
            const auto largestIndex=static_cast<std::size_t>(std::distance(vector.begin(),
                std::max_element(vector.begin(),vector.end(),[](const auto left,const auto right) {
                    return std::abs(left)<std::abs(right);
                })));
            retained.push_back(active[largestIndex]);
        }
        if (retained == active) break; active = std::move(retained);
    }
    vector = solve(); std::vector<std::complex<double>> result(columns);
    for (std::size_t index = 0; index < active.size(); ++index)
        result[active[index]] = vector[index] / scales[active[index]];
    double norm=0.0;for(const auto value:result)norm+=std::norm(value);
    norm=std::sqrt(norm);
    if (norm == 0.0) throw std::runtime_error("SINDy null vector vanished");
    for(auto& value:result)value/=norm;
    return result;
}

bool conjugatePartner(std::complex<double> left, std::complex<double> right) {
    return std::abs(right - std::conj(left)) <= 1e-9 + 1e-7 * std::abs(left);
}

}  // namespace


KoopmanSINDyResult fitKoopmanEigenfunctions(
    const std::vector<double>& stateValues, std::size_t samples, std::size_t stateDimension,
    const std::vector<double>& time, const std::vector<std::complex<double>>& eigenvalues,
    const std::vector<double>& velocityValues, const KoopmanSINDyOptions& options) {
    if (samples < 3 || !stateDimension || stateValues.size() != samples * stateDimension
        || time.size() != samples || eigenvalues.empty() || !options.polynomialDegree
        || !options.maximumIterations || options.threshold < 0.0)
        throw std::invalid_argument("invalid Koopman-with-SINDy inputs");
    for (std::size_t index = 1; index < time.size(); ++index)
        if (!(time[index] > time[index-1])) throw std::invalid_argument("SINDy time must increase");
    std::vector<double> center(stateDimension, 0.0), scale(stateDimension, 1.0);
    if (options.normalizeStates) for (std::size_t axis = 0; axis < stateDimension; ++axis) {
        for(std::size_t sample=0;sample<samples;++sample)
            center[axis]+=stateValues[sample*stateDimension+axis];
        center[axis]/=static_cast<double>(samples);
        double squared=0.0;
        for(std::size_t sample=0;sample<samples;++sample) {
            const double value=stateValues[sample*stateDimension+axis]-center[axis];
            squared+=value*value;
        }
        scale[axis]=std::sqrt(squared/static_cast<double>(samples));
        if (scale[axis] <= 1e-15) scale[axis] = 1.0;
    }
    const auto exponents = polynomialExponentTable(
        stateDimension, options.polynomialDegree, options.includeConstant);
    const auto theta = featureMatrix(
        stateValues,samples,stateDimension,exponents,center,scale);
    std::vector<double> velocities;
    if (velocityValues.empty()) velocities = finiteDifference(
        stateValues,samples,stateDimension,time);
    else {
        if (velocityValues.size() != stateValues.size()) throw std::invalid_argument("SINDy velocity shape mismatch");
        velocities = velocityValues;
    }
    const auto gamma = featureDerivative(stateValues,velocities,samples,
        stateDimension,exponents,center,scale);
    KoopmanSINDyResult result;
    result.continuousEigenvalues = eigenvalues; result.samples = samples;
    result.stateDimension = stateDimension; result.dynamicalStateDimension=stateDimension;
    result.featureCount = exponents.size() / stateDimension;
    result.exponents = exponents; result.center = center; result.scale = scale;
    const double dt=time[1]-time[0];
    bool uniform=true;
    for(std::size_t k=1;k<time.size();++k)
        if(std::abs((time[k]-time[k-1])-dt)>1e-7*dt) {uniform=false;break;}
    if(uniform) result.sampleTime=dt;
    result.coefficients.resize(eigenvalues.size() * result.featureCount);
    result.residualRms.resize(eigenvalues.size()); result.residualMaximum.resize(eigenvalues.size());
    std::vector<bool> assigned(eigenvalues.size(), false);
    for (std::size_t mode = 0; mode < eigenvalues.size(); ++mode) {
        if (assigned[mode]) continue;
        std::vector<std::complex<double>> generator(samples*result.featureCount);
        for(std::size_t index=0;index<generator.size();++index)
            generator[index]=gamma[index]-eigenvalues[mode]*theta[index];
        const auto coefficient = sparseNullVector(generator,samples,
            result.featureCount,options.threshold,options.maximumIterations);
        std::copy(coefficient.begin(), coefficient.end(),
                  result.coefficients.begin() + static_cast<std::ptrdiff_t>(mode * result.featureCount));
        assigned[mode] = true; ++result.directFits;
        if (options.reuseConjugateFits && std::abs(eigenvalues[mode].imag()) > 1e-12)
            for (std::size_t partner = mode + 1; partner < eigenvalues.size(); ++partner)
                if (!assigned[partner] && conjugatePartner(eigenvalues[mode], eigenvalues[partner])) {
                    for (std::size_t feature = 0; feature < result.featureCount; ++feature)
                        result.coefficients[partner * result.featureCount + feature]
                            = std::conj(coefficient[feature]);
                    assigned[partner] = true; break;
                }
    }
    if (options.materializeTrainingValues) result.trainingValues.resize(samples * eigenvalues.size());
    for (std::size_t mode = 0; mode < eigenvalues.size(); ++mode) {
        double residualSquared=0.0;
        for (std::size_t sample = 0; sample < samples; ++sample) {
            std::complex<double> residual{},value{};
            for(std::size_t feature=0;feature<result.featureCount;++feature) {
                const auto coefficient=
                    result.coefficients[mode*result.featureCount+feature];
                value+=theta[sample*result.featureCount+feature]*coefficient;
                residual+=(gamma[sample*result.featureCount+feature]
                    -eigenvalues[mode]*theta[sample*result.featureCount+feature])
                    *coefficient;
            }
            const double magnitude=std::abs(residual);
            residualSquared+=magnitude*magnitude;
            result.residualMaximum[mode]=std::max(
                result.residualMaximum[mode],magnitude);
            if (options.materializeTrainingValues)
                result.trainingValues[sample * eigenvalues.size() + mode] = value;
        }
        result.residualRms[mode] = std::sqrt(residualSquared / samples);
    }
    return result;
}

KoopmanSINDyResult fitKoopmanEigenfunctionsInCoordinates(
    const std::vector<double>& dynamicalStates, std::size_t samples,
    std::size_t dynamicalStateDimension, const std::vector<double>& time,
    const std::vector<std::complex<double>>& eigenvalues,
    const Observable& coordinates, const std::vector<double>& dynamicalVelocities,
    const KoopmanSINDyOptions& options) {
    if (coordinates.stateDimension()!=dynamicalStateDimension)
        throw std::invalid_argument("SINDy coordinate observable is not state aligned");
    auto values=coordinates.evaluateRealBatch(dynamicalStates,samples);
    const auto velocities=dynamicalVelocities.empty()?std::vector<double>{}
        :coordinates.pushForwardRealBatch(dynamicalStates,dynamicalVelocities,samples);
    auto result=fitKoopmanEigenfunctions(values,samples,coordinates.outputDimension(),time,
                                         eigenvalues,velocities,options);
    result.dynamicalStateDimension=dynamicalStateDimension;
    if (!coordinates.expressions().empty())
        result.coordinateObservableDefinition=json::dump(coordinates.expressionDefinition());
    result.coordinateObservableNames=coordinates.descriptor().outputNames;
    if (result.coordinateObservableNames.empty() && coordinates.outputDimension()==1)
        result.coordinateObservableNames.push_back(coordinates.descriptor().name);
    return result;
}


}  // namespace dynsyslab::core
