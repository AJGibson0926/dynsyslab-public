#include "KromEvidence.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

dynsyslab::core::NpyArray selectedColumns(
    const dynsyslab::core::NpyArray& source,
    const std::vector<std::size_t>& indices) {
    if (indices.empty()) return source;
    if (source.shape.size()!=2)
        throw std::invalid_argument("selected KROM evidence must be a matrix");
    const auto values=source.complexDoubles();
    std::vector<std::complex<double>> selected;
    selected.reserve(source.shape[0]*indices.size());
    for (std::size_t row=0;row<source.shape[0];++row)
        for (const auto column:indices) {
            if (column>=source.shape[1])
                throw std::out_of_range("selected KROM mode exceeds archived values");
            selected.push_back(values[row*source.shape[1]+column]);
        }
    return dynsyslab::core::NpyArray::fromComplex(
        {source.shape[0],indices.size()},selected);
}

dynsyslab::core::NpyArray selectedVector(
    const dynsyslab::core::NpyArray& source,
    const std::vector<std::size_t>& indices) {
    if (indices.empty()) return source;
    if (source.shape.size()!=1)
        throw std::invalid_argument("selected KROM spectrum must be a vector");
    const auto values=source.complexDoubles();
    std::vector<std::complex<double>> selected;
    for (const auto index:indices) {
        if (index>=values.size())
            throw std::out_of_range("selected KROM mode exceeds archived spectrum");
        selected.push_back(values[index]);
    }
    return dynsyslab::core::NpyArray::fromComplex({selected.size()},selected);
}

}  // namespace

std::vector<std::size_t> inheritedSourceModes(
    const dynsyslab::core::NpzRecord& parent,
    const std::vector<std::size_t>& selectedModes,
    std::size_t outputModes) {
    std::vector<std::size_t> parentModes;
    if (parent.contains("source_mode_indices"))
        for (const auto value:parent.at("source_mode_indices").doubles())
            parentModes.push_back(static_cast<std::size_t>(std::llround(value)));
    if (parentModes.empty()) {
        std::size_t count=outputModes;
        if (!selectedModes.empty())
            count=*std::max_element(selectedModes.begin(),selectedModes.end())+1;
        parentModes.resize(count);
        for (std::size_t index=0;index<count;++index) parentModes[index]=index;
    }
    if (selectedModes.empty()) {
        if (parentModes.size()!=outputModes)
            throw std::invalid_argument(
                "archived source-mode identity map does not match KROM dimension");
        return parentModes;
    }
    std::vector<std::size_t> result;
    for (const auto index:selectedModes) {
        if (index>=parentModes.size())
            throw std::out_of_range("selected KROM mode exceeds source identity map");
        result.push_back(parentModes[index]);
    }
    return result;
}

void inheritKromEvidence(
    const dynsyslab::core::NpzRecord& parent,
    const std::vector<std::size_t>& selectedModes,
    std::map<std::string,dynsyslab::core::NpyArray>& output) {
    for (const auto* name:{"training_states","training_coordinates","training_time",
            "training_coordinate_velocities","training_coordinate_jacobians"})
        if (parent.contains(name)) output[name]=parent.at(name);
    if (parent.contains("eigenfunction_values"))
        output["eigenfunction_values"]=selectedColumns(
            parent.at("eigenfunction_values"),selectedModes);
    const auto* discrete=parent.contains("discrete_eigenvalues")
        ? &parent.at("discrete_eigenvalues")
        : (parent.contains("eigenvalues")?&parent.at("eigenvalues"):nullptr);
    if (discrete) {
        output["discrete_eigenvalues"]=selectedVector(*discrete,selectedModes);
        output["eigenvalues"]=output["discrete_eigenvalues"];
    }
}

}  // namespace dynsyslab::workflow::postprocess_detail
