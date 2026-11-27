#pragma once
// Geometry of sampled orbits in a plane of coordinates, for injectivity audits of embeddings
// (a delay or observable map is injective on an attractor only if distinct states never share a
// point; Takens 1981, Lecture Notes in Math. 898:366-381; Sauer, Yorke and Casdagli 1991,
// J. Stat. Phys. 65:579-616). Two orbits whose projections cross at different physical states
// witness a non-injective coordinate map.
//
//   findPeaks             local maxima with a minimum prominence and spacing (topographic
//                         prominence: height above the higher of the two lowest points that
//                         separate the peak from higher ground on either side)
//   representativeCycle   one late orbit period between consecutive peaks, resampled uniformly
//   peakClosure           median |z(p_{k+1}) - z(p_k)| over peaks, relative to the RMS orbit radius
//                         (a periodic orbit closes exactly; a drifting one does not)
//   findCurveIntersections  proper crossings of two polylines by the parametric segment test
//                         (O'Rourke, Computational Geometry in C, 2nd ed. 1998, section 7.2)
//   polygonArea           |shoelace sum| / 2 of a closed sampled curve
//   nearestSeparation     smallest distance between the samples of two curves
//   refineVisits          the two visit times of a crossing witness located on an independent
//                         (finer) record: z(t1) = z(t2) solved by Newton's method on its linear
//                         interpolant inside the witness's brackets
//   matchVisits           one-to-one matching of visit pairs of two records within time
//                         tolerances, maximum cardinality then least cost (Hungarian method,
//                         Kuhn 1955, Naval Res. Logist. Q. 2:83-97; Jonker and Volgenant 1987,
//                         Computing 38:325-340)
#include <cstddef>
#include <utility>
#include <vector>

namespace dynsyslab::core {

// Indices of local maxima of y whose prominence is at least `prominence` and that are at least
// `distance` samples apart (lower peaks within `distance` of a higher one are dropped). Plateaus
// count once, at their first sample.
[[nodiscard]] std::vector<std::size_t> findPeaks(const std::vector<double>& y, std::size_t distance, double prominence);

// The segment of a (samples x width) row-major path between two consecutive peaks, linearly
// resampled to `samples` points (rows x width). `which` counts cycles back from the last complete
// one (0: the last). Throws when the peaks give no such cycle.
[[nodiscard]] std::vector<double> representativeCycle(const std::vector<double>& path, std::size_t width,
                                                      const std::vector<std::size_t>& peaks, std::size_t samples,
                                                      std::size_t which = 0);

// Median distance between the path at consecutive peaks over its RMS radius about its mean;
// NaN with fewer than three peaks.
[[nodiscard]] double peakClosure(const std::vector<double>& path, std::size_t width, const std::vector<std::size_t>& peaks);

struct CurveIntersection {
    std::size_t first{}, second{};  // segment indices on curve a and curve b
    double u{}, v{};                // fractions along them, strictly inside (0, 1)
    double x{}, y{};                // the crossing point
};

// Proper, nonparallel crossings of polylines a and b given as (x, y) sample pairs. Crossings
// closer than `merge` (in the plane) to an earlier one are merged.
[[nodiscard]] std::vector<CurveIntersection> findCurveIntersections(const std::vector<double>& ax, const std::vector<double>& ay,
                                                                    const std::vector<double>& bx, const std::vector<double>& by,
                                                                    double merge = 0);

[[nodiscard]] double polygonArea(const std::vector<double>& x, const std::vector<double>& y);

// min_{i,j} |a_i - b_j| over the first `width` columns of two row-major paths.
[[nodiscard]] double nearestSeparation(const std::vector<double>& a, const std::vector<double>& b, std::size_t width);

struct RefinedVisit {
    double first{}, second{};  // visit times
    double gap{};              // |z(first) - z(second)| on the record, in scaled coordinates
    bool matched{};            // converged inside the brackets with gap <= tolerance
};

// time (samples), x and y (samples): the record; brackets [a0, a1] and [b0, b1] with a1 <= b0.
[[nodiscard]] RefinedVisit refineVisits(const std::vector<double>& time, const std::vector<double>& x, const std::vector<double>& y,
                                        double a0, double a1, double b0, double b1, double scaleX = 1, double scaleY = 1,
                                        double tolerance = 1e-9);

struct VisitMatching {
    std::vector<std::pair<std::size_t, std::size_t>> matches;
    std::vector<std::size_t> unmatchedFirst, unmatchedSecond;
};

// Visit pairs (t1, t2) of two records; a match needs |dt1| <= tolerance1 and |dt2| <= tolerance2.
[[nodiscard]] VisitMatching matchVisits(const std::vector<std::pair<double, double>>& first,
                                        const std::vector<std::pair<double, double>>& second, double tolerance1, double tolerance2);

}  // namespace dynsyslab::core
