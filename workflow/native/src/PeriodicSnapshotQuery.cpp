#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/PeriodicSnapshots.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double num(const J&v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");return v.number();}
const A&array(const J&v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>1000000)throw std::invalid_argument("Bounded array required");return v.array();}
std::vector<double>numbers(const J&v){std::vector<double>r;for(const auto&e:array(v))r.push_back(num(e));return r;}
}
core::json::Value evaluatePeriodicSnapshotQuery(const J&q){
    std::vector<math::PeriodicSnapshotOrbit>orbits;A frequencies,masses;
    for(const auto&o:array(q.at("orbits"))){math::PeriodicSnapshotOrbit orbit;orbit.weight=num(o.at("weight"));
        for(const auto&channel:array(o.at("channels"))){math::FourierSeries s{num(channel.at("angular_frequency")),{}};
            for(const auto&row:array(channel.at("coefficients"))){const auto&c=array(row);if(c.size()!=3)throw std::invalid_argument("Fourier row needs [n,re,im]");
                double n=num(c[0]);if(n!=std::floor(n)||std::abs(n)>1000000)throw std::invalid_argument("Invalid harmonic");
                if(!s.coefficients.emplace(int(n),std::complex<double>(num(c[1]),num(c[2]))).second)throw std::invalid_argument("Duplicate harmonic");}
            orbit.channels.push_back(std::move(s));}
        if(orbit.channels.empty())throw std::invalid_argument("Empty orbit channels");
        frequencies.emplace_back(orbit.channels[0].angularFrequency);masses.emplace_back(orbit.weight);orbits.push_back(std::move(orbit));}
    const double count=num(q.at("phase_count"));if(count<3||count>1000000||count!=std::floor(count))throw std::invalid_argument("Invalid phase_count");
    const double dt=num(q.at("sample_interval")),phase=q.at("phase_origin").isNull()?0:num(q.at("phase_origin"));
    if(!q.at("split_complex").isNull()&&!std::holds_alternative<bool>(q.at("split_complex").storage()))throw std::invalid_argument("split_complex must be boolean");
    const bool split=q.at("split_complex").boolean();const auto offsets=numbers(q.at("time_offsets"));
    const auto pairs=math::periodicSnapshotPairs(orbits,std::size_t(count),offsets,dt,split,phase);
    const auto path=q.at("output_archive").string();if(path.empty())throw std::invalid_argument("output_archive is required");
    O metadata{{"schema","dynsyslab.koopman-snapshot-pairs/v1"},{"generator","periodic-fourier-product-quadrature/v1"},
        {"sample_interval",dt},{"time_offsets",q.at("time_offsets")},{"phase_count",count},{"phase_origin",phase},
        {"orbit_frequencies",frequencies},{"orbit_masses",masses},{"split_complex",split},
        {"layout","rows orbit/phase; columns offset/channel/real-imag"},{"source",q.at("source")},
        {"measure",q.at("measure").string("Explicit orbit masses times normalized uniform phase; masses not automatically normalized")}};
    core::NpzRecord archive;archive.arrays["x"]=core::NpyArray::fromDoubles({pairs.samples,pairs.dimension},pairs.x);
    archive.arrays["y"]=core::NpyArray::fromDoubles({pairs.samples,pairs.dimension},pairs.y);
    archive.arrays["weights"]=core::NpyArray::fromDoubles({pairs.samples},pairs.weights);
    archive.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(metadata));
    core::writeNpzAtomicNew(std::filesystem::u8path(path),archive);
    return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-snapshot-pairs"},
        {"archive",path},{"samples",double(pairs.samples)},{"dimension",double(pairs.dimension)},{"metadata",metadata}};
}
}
