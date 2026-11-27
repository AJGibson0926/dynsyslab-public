#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
void writeFrozenHistoryTransform(NpzRecord& a,const FrozenHistoryTransform& transform) {
    if(!(transform.sampleInterval()>0))throw std::invalid_argument("Published history transform needs a sample interval");
    const auto& b=transform.basis();
    a.arrays["left_singular_vectors"]=NpyArray::fromDoubles({b.rows,b.rank},b.leftSingularVectors);
    a.arrays["singular_values"]=NpyArray::fromDoubles({b.rank},b.singularValues);
    if(!b.cumulativeEnergy.empty())a.arrays["cumulative_energy"]=NpyArray::fromDoubles({b.cumulativeEnergy.size()},b.cumulativeEnergy);
    a.arrays["observation_scale"]=NpyArray::fromDoubles({1},{b.observationScale});
    a.arrays["history_sample_time"]=NpyArray::fromDoubles({1},{transform.sampleInterval()});
    a.arrays["history_delay_samples"]=NpyArray::fromDoubles({1},{double(b.rows)});
    json::Value::Array coordinates;for(auto n:transform.retainedCoordinates())coordinates.emplace_back(double(n));
    a.arrays["history_transform_json"]=NpyArray::fromString(json::dump(json::Value::Object{
        {"schema","dynsyslab.history-transform/v1"},{"identity",transform.identity()},
        {"sample_interval",transform.sampleInterval()},{"retained_coordinates",coordinates},
        {"measurement",transform.measurement()},{"time_units",transform.timeUnits()},
        {"observation_units",transform.observationUnits()},
        {"training_trajectories",double(b.trajectories)},{"training_windows",double(b.columns)},
        {"total_hankel_energy",b.totalEnergy}}));
}
FrozenHistoryTransform readFrozenHistoryTransform(const NpzRecord& a) {
    const auto& u=a.at("left_singular_vectors");
    if(u.shape.size()!=2)throw std::invalid_argument("History U must be delay x rank");
    HankelLeftBasis b;b.rows=u.shape[0];b.rank=u.shape[1];b.leftSingularVectors=u.doubles();
    b.singularValues=a.at("singular_values").doubles();
    if(a.contains("cumulative_energy"))b.cumulativeEnergy=a.at("cumulative_energy").doubles();
    if(a.contains("observation_scale"))b.observationScale=a.at("observation_scale").doubles().at(0);
    const auto meta=a.metadataJson().empty()?json::Value{}:json::parse(a.metadataJson());
    std::vector<std::size_t> coordinates;
    auto size=[](const json::Value& v) {
        const double n=v.number(-1);
        if(!std::isfinite(n)||n<0||n!=std::floor(n)||n>double(std::numeric_limits<std::size_t>::max()/2))
            throw std::invalid_argument("Invalid history size/index");
        return static_cast<std::size_t>(n);
    };
    if(a.contains("history_transform_json")) {
        const auto spec=json::parse(a.at("history_transform_json").string());
        if(spec.at("schema").string()!="dynsyslab.history-transform/v1")throw std::invalid_argument("Unsupported history transform schema");
        for(const auto& n:spec.at("retained_coordinates").array())coordinates.push_back(size(n));
        b.trajectories=size(spec.at("training_trajectories"));b.columns=size(spec.at("training_windows"));
        b.totalEnergy=spec.at("total_hankel_energy").number();
        FrozenHistoryTransform transform(std::move(b),spec.at("sample_interval").number(-1),std::move(coordinates),
            spec.at("measurement").string(),spec.at("time_units").string(),spec.at("observation_units").string());
        if(transform.sampleInterval()<=0||transform.identity()!=spec.at("identity").string())
            throw std::invalid_argument("History transform content/clock differs from its recorded identity");
        return transform;
    }
    // Existing native echo archives already store the complete frozen basis.
    // Its content identity is computed during import, never claimed from a path.
    const auto rank=meta.at("history_coordinate_rank").isNull()?b.rank:size(meta.at("history_coordinate_rank"));
    if(!rank||rank>b.rank)throw std::invalid_argument("Invalid archived history coordinate rank");
    for(std::size_t j=0;j<rank;++j)coordinates.push_back(j);
    const double dt=a.at("history_sample_time").doubles().at(0);
    if(!std::isfinite(dt)||dt<=0)throw std::invalid_argument("Archived history interval must be positive");
    return FrozenHistoryTransform(std::move(b),dt,std::move(coordinates),
        meta.at("history_observation").string(),meta.at("history_time_units").string("model-time"));
}
FrozenHistoryTransform loadFrozenHistoryTransform(const std::filesystem::path& path) {
    return readFrozenHistoryTransform(readNpz(path));
}
void saveFrozenHistoryTransform(const std::filesystem::path& path,
    const FrozenHistoryTransform& transform,const std::string& metadata) {
    auto meta=json::parse(metadata);
    if(!std::holds_alternative<json::Value::Object>(meta.storage()))throw std::invalid_argument("Transform metadata must be an object");
    meta.object()["schema"]="dynsyslab.derived/v1";meta.object()["method"]="history-transform";
    meta.object()["transform_identity"]=transform.identity();
    NpzRecord a;writeFrozenHistoryTransform(a,transform);
    a.arrays["metadata_json"]=NpyArray::fromString(json::dump(meta));writeNpzAtomic(path,a);
}
}
