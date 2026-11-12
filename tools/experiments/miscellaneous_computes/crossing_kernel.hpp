#pragma once
#include <cstddef>
struct mc_crossing_v1 {
    std::size_t first,second;
    double u,v; // fractions along the two time-ordered segments
};
// Two-column row-major coordinates. Reports isolated nonparallel segment
// intersections. Pass nullptr for a counting pass. Status 0 success, 1 bad
// input, 2 insufficient output capacity, 3 allocation failure.
extern "C" int mc_crossings_v1(const double* xy,std::size_t count,
    mc_crossing_v1* out,std::size_t capacity,std::size_t* found);
