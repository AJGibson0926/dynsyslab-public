#pragma once
#include "JetPolynomial.hpp"

namespace dynsyslab::core::jet {
std::string fiberProblem(const std::vector<Polynomial>& jet,const std::vector<Polynomial>& different,
    std::size_t dimension,const std::set<std::size_t>& poles,const std::vector<std::size_t>& positiveAxes,
    const std::vector<Polynomial>& graph={});
json::Value decideFiber(const std::string& problem,std::size_t timeoutMs);
} // namespace dynsyslab::core::jet
