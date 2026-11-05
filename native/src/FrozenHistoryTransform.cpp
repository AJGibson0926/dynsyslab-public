#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
namespace dynsyslab::core {
FrozenHistoryTransform::FrozenHistoryTransform(PooledHankelBasis basis,double dt,
    std::vector<std::size_t> coordinates,std::string measurement,
    std::string timeUnits,std::string observationUnits)
    :basis_(std::make_shared<const PooledHankelBasis>(std::move(basis))),sampleInterval_(dt),
     coordinates_(std::move(coordinates)),measurement_(std::move(measurement)),
     timeUnits_(std::move(timeUnits)),observationUnits_(std::move(observationUnits)) {
    const auto& b=*basis_;
    if(b.delays<2||!b.rank||b.rank>b.delays||b.rank>std::numeric_limits<std::size_t>::max()/b.delays
        ||b.leftSingularVectors.size()!=b.rank*b.delays||b.singularValues.size()!=b.rank
        ||!std::isfinite(dt)||dt<0||!std::isfinite(b.observationScale)||b.observationScale<=0)
        throw std::invalid_argument("Invalid frozen history transform dimensions/scales/clock");
    for(double v:b.leftSingularVectors)if(!std::isfinite(v))throw std::invalid_argument("History basis must be finite");
    for(double v:b.singularValues)if(!std::isfinite(v)||v<=std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("History singular values must be positive and finite");
    if(coordinates_.empty()){coordinates_.resize(b.rank);std::iota(coordinates_.begin(),coordinates_.end(),0);}
    const std::set<std::size_t> unique(coordinates_.begin(),coordinates_.end());
    if(unique.size()!=coordinates_.size()||*unique.rbegin()>=b.rank)
        throw std::invalid_argument("History coordinate selection is invalid");
    auto numbers=[](const auto& values){json::Value::Array a;for(auto v:values)a.emplace_back(double(v));return a;};
    identity_="sha256:"+sha256Text(json::dump(json::Value::Object{
        {"schema","dynsyslab.history-transform/v1"},{"order","oldest-to-newest"},
        {"projection","U-transpose-h-over-sigma"},{"delays",double(b.delays)},{"rank",double(b.rank)},
        {"basis",numbers(b.leftSingularVectors)},{"singular_values",numbers(b.singularValues)},
        {"fitting_observation_scale",b.observationScale},{"sample_interval",dt},
        {"retained_coordinates",numbers(coordinates_)},{"measurement",measurement_},
        {"time_units",timeUnits_},{"observation_units",observationUnits_}}));
}
}
