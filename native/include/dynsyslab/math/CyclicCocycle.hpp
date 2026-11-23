#pragma once
#include <dynsyslab/core/Json.hpp>
#include <vector>
#include <cstddef>
namespace dynsyslab::math {
// Finite periodic matrix cocycle, phase-major row-major matrices. Quotient
// projection is Haar conditional expectation, not a physical renormalization.
core::json::Value cyclicCocycle(const std::vector<double>& matrices,std::size_t dimension,
    std::size_t quotientPeriod,const std::vector<double>& initial,std::size_t steps,std::size_t phase=0);
// Sampled nested residue partitions; cannot certify an infinite inverse limit.
core::json::Value cyclicPartitionGeometry(const std::vector<double>& points,std::size_t dimension,
    std::size_t base=2,const std::vector<double>& images={});
// Explicit divisibility tower, including mixed-prime and odd-base doubling
// hierarchies. Every period divides the sample cycle; infinity is not inferred.
core::json::Value cyclicTowerPartitionGeometry(const std::vector<double>& points,std::size_t dimension,
    const std::vector<std::size_t>& periods,const std::vector<double>& images={});
core::json::Value evaluateCyclicCocycleQuery(const core::json::Value&);
}
