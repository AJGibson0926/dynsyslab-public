#include "lookup.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace receivers {
namespace {
bool closer(const Neighbor& a, const Neighbor& b) {
    return a.squaredDistance == b.squaredDistance ? a.row < b.row : a.squaredDistance < b.squaredDistance;
}
}
Lookup::Lookup(std::vector<double> features, std::vector<double> states,
               std::vector<double> scales, std::size_t stateDimension)
    : points_(std::move(features)), states_(std::move(states)), scales_(std::move(scales)),
      stateDimension_(stateDimension) {
    const auto d = scales_.size();
    if (!d || !stateDimension_ || points_.empty() || points_.size()%d ||
        states_.size()/stateDimension_ != points_.size()/d || states_.size()%stateDimension_)
        throw std::invalid_argument("Training feature and state matrix dimensions disagree");
    for (double s : scales_) if (!std::isfinite(s) || s <= 0) throw std::invalid_argument("Scales must be finite and positive");
    for (double s : states_) if (!std::isfinite(s)) throw std::invalid_argument("Nonfinite training label");
    for (std::size_t i = 0; i < points_.size(); ++i) {
        points_[i] = std::asinh(points_[i]/scales_[i%d]);
        if (!std::isfinite(points_[i])) throw std::invalid_argument("Nonfinite scaled training feature");
    }
    order_.resize(points_.size()/d);
    std::iota(order_.begin(),order_.end(),std::size_t{0});
    build(0,order_.size(),0);
}
void Lookup::build(std::size_t first, std::size_t last, std::size_t axis) {
    if (last-first < 2) return;
    const auto mid=first+(last-first)/2, d=scales_.size();
    std::nth_element(order_.begin()+first,order_.begin()+mid,order_.begin()+last,[&](auto a,auto b) {
        const double x=points_[a*d+axis], y=points_[b*d+axis];
        return x==y ? a<b : x<y;
    });
    build(first,mid,(axis+1)%d);
    build(mid+1,last,(axis+1)%d);
}
void Lookup::search(std::size_t first, std::size_t last, std::size_t axis,
                    const std::vector<double>& query, std::size_t k, std::vector<Neighbor>& hits) const {
    if (first==last) return;
    const auto mid=first+(last-first)/2, row=order_[mid], d=scales_.size();
    const double* point=points_.data()+row*d;
    double square=0;
    for (std::size_t j=0;j<d;++j) {const double diff=query[j]-point[j]; square+=diff*diff;}
    const Neighbor hit{square,row};
    if (hits.size()<k || closer(hit,hits.back())) {
        if (hits.size()==k) hits.pop_back();
        hits.insert(std::lower_bound(hits.begin(),hits.end(),hit,closer),hit);
    }
    const double delta=query[axis]-point[axis];
    const auto next=(axis+1)%d;
    if (delta<=0) {
        search(first,mid,next,query,k,hits);
        if (hits.size()<k || delta*delta<=hits.back().squaredDistance) search(mid+1,last,next,query,k,hits);
    } else {
        search(mid+1,last,next,query,k,hits);
        if (hits.size()<k || delta*delta<=hits.back().squaredDistance) search(first,mid,next,query,k,hits);
    }
}
Estimate Lookup::predict(const double* features, std::size_t k) const {
    if (!features || !k || k>order_.size()) throw std::invalid_argument("Invalid query or neighbor count");
    std::vector<double> query(scales_.size());
    for (std::size_t j=0;j<query.size();++j) {
        query[j]=std::asinh(features[j]/scales_[j]);
        if (!std::isfinite(query[j])) throw std::invalid_argument("Nonfinite scaled query feature");
    }
    Estimate result;
    result.neighbors.reserve(k);
    search(0,order_.size(),0,query,k,result.neighbors);
    const double nearest=std::max(1e-24,result.neighbors.front().squaredDistance);
    double sum=0;
    for (const auto& hit : result.neighbors) {
        const double weight=nearest/std::max(1e-24,hit.squaredDistance);
        result.weights.push_back(weight);sum+=weight;
    }
    result.values.assign(stateDimension_,0.);
    for (std::size_t i=0;i<k;++i) {
        result.weights[i]/=sum;
        for (std::size_t j=0;j<stateDimension_;++j)
            result.values[j]+=result.weights[i]*states_[result.neighbors[i].row*stateDimension_+j];
    }
    return result;
}
}
