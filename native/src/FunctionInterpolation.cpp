#include <dynsyslab/core/FunctionInterpolation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {

struct FunctionInterpolation::Index {
    struct Node { std::size_t sample{}, axis{}; int left{-1}, right{-1}; };
    const InterpolationData* data{};
    std::vector<Node> nodes;
    int root{-1};

    explicit Index(const InterpolationData& source) : data(&source) {
        std::vector<std::size_t> samples(source.samples);
        for (std::size_t index=0;index<samples.size();++index)samples[index]=index;
        root=build(samples,0,samples.size(),0);
    }

    double coordinate(std::size_t sample,std::size_t axis) const {
        return data->trainingStates[sample*data->stateDimension+axis];
    }

    int build(std::vector<std::size_t>& samples,std::size_t first,
              std::size_t last,std::size_t depth) {
        if(first==last)return -1;const std::size_t axis=depth%data->stateDimension;
        const std::size_t middle=first+(last-first)/2;
        std::nth_element(samples.begin()+static_cast<std::ptrdiff_t>(first),
            samples.begin()+static_cast<std::ptrdiff_t>(middle),
            samples.begin()+static_cast<std::ptrdiff_t>(last),
            [&](std::size_t left,std::size_t right){
                return coordinate(left,axis)<coordinate(right,axis);});
        const int node=static_cast<int>(nodes.size());nodes.push_back({samples[middle],axis});
        nodes[static_cast<std::size_t>(node)].left=build(samples,first,middle,depth+1);
        nodes[static_cast<std::size_t>(node)].right=build(samples,middle+1,last,depth+1);
        return node;
    }

    double distanceSquared(std::size_t sample,const double* point) const {
        double result=0.0;for(std::size_t axis=0;axis<data->stateDimension;++axis){
            const double delta=coordinate(sample,axis)-point[axis];result+=delta*delta;}
        return result;
    }

    using Neighbor=std::pair<double,std::size_t>;
    void queryNode(int node,const double* point,std::size_t count,
                   std::priority_queue<Neighbor>& heap) const {
        if(node<0)return;const auto& item=nodes[static_cast<std::size_t>(node)];
        const double distance=distanceSquared(item.sample,point);
        if(heap.size()<count)heap.push({distance,item.sample});
        else if(distance<heap.top().first){heap.pop();heap.push({distance,item.sample});}
        const double split=point[item.axis]-coordinate(item.sample,item.axis);
        const int near=split<0.0?item.left:item.right;
        const int far=split<0.0?item.right:item.left;
        queryNode(near,point,count,heap);
        if(heap.size()<count||split*split<heap.top().first)
            queryNode(far,point,count,heap);
    }

    std::vector<Neighbor> query(const double* point,std::size_t count) const {
        std::priority_queue<Neighbor> heap;queryNode(root,point,count,heap);
        std::vector<Neighbor> result(heap.size());
        for(auto position=result.rbegin();position!=result.rend();++position){
            *position=heap.top();heap.pop();}
        return result;
    }
};

FunctionInterpolation::FunctionInterpolation(InterpolationData data)
    : data_(std::move(data)) {
    const std::size_t modes=data_.outputs();
    if(!data_.samples||!data_.stateDimension||!modes
        ||data_.trainingStates.size()!=data_.samples*data_.stateDimension
        ||data_.trainingValues.size()!=data_.samples*modes
        ||!data_.neighbors||!(data_.finiteDifferenceRelativeStep>0.0))
        throw std::invalid_argument("scattered eigenfunction data are incomplete");
    (void)interpolationMethodFromJson(interpolationMethodJson(data_.method,data_.neighbors));  // validates the method
    data_.neighbors=std::min(data_.neighbors,data_.samples);
    for(const double value:data_.trainingStates)if(!std::isfinite(value))
        throw std::invalid_argument("scattered eigenfunction states must be finite");
    for(const auto value:data_.trainingValues)
        if(!std::isfinite(value.real())||!std::isfinite(value.imag()))
            throw std::invalid_argument("scattered eigenfunction values must be finite");
    index_=std::make_shared<Index>(data_);
    if(data_.coordinateObservable.empty()){
        coordinates_=FunctionDomain::physicalState(data_.stateDimension);
    }else{  // catalog observables, e.g. "echo,echo_dot"
        std::vector<std::string> names;const auto& bound=data_.coordinateObservable;std::size_t start=0;
        for(std::size_t comma;(comma=bound.find(',',start))!=std::string::npos;start=comma+1)
            names.push_back(bound.substr(start,comma-start));
        names.push_back(bound.substr(start));
        coordinates_=FunctionDomain::fromCatalog(std::move(names),data_.stateDimension);
    }
}

FunctionInterpolation::FunctionInterpolation(InterpolationData data, FunctionDomain coordinates)
    : FunctionInterpolation(std::move(data)) {
    if(coordinates.dimension!=data_.stateDimension)
        throw std::invalid_argument("interpolation domain dimension differs from its points");
    coordinates_=std::move(coordinates);
}

std::shared_ptr<const FunctionInterpolation> FunctionInterpolation::withMethod(InterpolationMethod method,
                                                                             std::size_t neighbors) const {
    InterpolationData data=data_;
    data.method=std::move(method);
    data.neighbors=neighbors;
    return std::make_shared<FunctionInterpolation>(std::move(data),coordinates_);
}

std::shared_ptr<const FunctionInterpolation> interpolateSamples(const std::vector<double>& states,
    std::size_t samples, const FunctionDomain& domain, std::vector<std::complex<double>> values,
    InterpolationMethod method, std::size_t neighbors, std::string provenance) {
    if(!samples||states.size()%samples)throw std::invalid_argument("interpolateSamples: states are not samples x dimension");
    if(domain.needsCatalog())throw std::invalid_argument("interpolateSamples: bind the domain's catalog map first");
    const std::size_t width=states.size()/samples;
    InterpolationData data;
    data.samples=samples;
    data.stateDimension=domain.dimension;
    data.trainingValues=std::move(values);
    data.method=std::move(method);
    data.neighbors=neighbors;
    data.provenance=std::move(provenance);
    if(!domain.stateMap){  // supplied coordinates: the states are the base points
        if(width!=domain.dimension)throw std::invalid_argument("interpolateSamples: supplied points do not match the domain");
        data.trainingStates=states;
    }else{
        const auto& map=*domain.stateMap;
        if(map.stateDimension()!=width||map.outputDimension()!=domain.dimension)
            throw std::invalid_argument("interpolateSamples: the domain map does not fit the states");
        data.trainingStates=map.evaluateRealBatch(states,samples);
    }
    return std::make_shared<FunctionInterpolation>(std::move(data),domain);
}

std::shared_ptr<const FunctionRepresentation> FunctionInterpolation::onDomain(FunctionDomain coordinates) const {
    if(coordinates.dimension!=coordinates_.dimension)throw std::invalid_argument("interpolation domain dimension differs");
    // A fresh object (and index): the neighbour index points into its own sample arrays.
    auto out=std::make_shared<FunctionInterpolation>(data_);
    out->coordinates_=std::move(coordinates);
    return out;
}

const InterpolationData& FunctionInterpolation::data() const noexcept {
    return data_;
}

std::vector<std::complex<double>> FunctionInterpolation::evaluate(
    const double* state) const {
    for(std::size_t axis=0;axis<data_.stateDimension;++axis)
        if(!std::isfinite(state[axis]))throw std::invalid_argument(
            "scattered eigenfunction state must be finite");
    const auto neighbors=index_->query(state,data_.neighbors);
    if(data_.method.name!="inverse_distance"||data_.method.power!=2.0)
        return interpolateScatteredNeighbours(data_,neighbors,state);
    std::vector<std::complex<double>> result(modes());
    if(neighbors.front().first<=1.0e-26){
        const auto sample=neighbors.front().second;
        std::copy_n(data_.trainingValues.begin()+static_cast<std::ptrdiff_t>(sample*modes()),
            modes(),result.begin());return result;
    }
    double total=0.0;for(const auto& neighbor:neighbors){
        const double weight=1.0/std::max(neighbor.first,1.0e-26);total+=weight;
        for(std::size_t mode=0;mode<modes();++mode)
            result[mode]+=weight*data_.trainingValues[neighbor.second*modes()+mode];
    }
    for(auto& value:result)value/=total;return result;
}

std::vector<std::complex<double>> FunctionInterpolation::jacobian(
    const double* state) const {
    std::vector<std::complex<double>> result(modes()*data_.stateDimension);
    std::vector<double> plus(state,state+data_.stateDimension),minus(plus);
    for(std::size_t axis=0;axis<data_.stateDimension;++axis){
        const double step=data_.finiteDifferenceRelativeStep
            *std::max(1.0,std::abs(state[axis]));
        plus[axis]+=step;minus[axis]-=step;
        const auto high=evaluate(plus.data()),low=evaluate(minus.data());
        for(std::size_t mode=0;mode<modes();++mode)
            result[mode*data_.stateDimension+axis]=(high[mode]-low[mode])/(2.0*step);
        plus[axis]=state[axis];minus[axis]=state[axis];
    }
    return result;
}

std::shared_ptr<const FunctionInterpolation> FunctionInterpolation::select(
    const std::vector<std::size_t>& selected) const {
    if(selected.empty())throw std::invalid_argument("scattered KROM extraction needs modes");
    InterpolationData subset=data_;
    subset.trainingValues.assign(data_.samples*selected.size(),{});
    for(std::size_t output=0;output<selected.size();++output){
        if(selected[output]>=modes())throw std::out_of_range("scattered KROM mode index");
        for(std::size_t sample=0;sample<data_.samples;++sample)
            subset.trainingValues[sample*selected.size()+output]=
                data_.trainingValues[sample*modes()+selected[output]];
    }
    subset.provenance+=subset.provenance.empty()?"independent mode subset":
        "; independent mode subset";
    return std::make_shared<FunctionInterpolation>(std::move(subset),coordinates_);
}

}  // namespace dynsyslab::core
