#include <algorithm>
#include <cmath>
#include <dynsyslab/core/AnalyticSolution.hpp>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
AnalyticComparison compareAnalyticSolution(const std::vector<AnalyticSolutionFamily> &families,
                                           const AnalyticState &initial, double t0,
                                           const std::vector<double> &time,
                                           const std::vector<double> &numerical,
                                           const std::string &selection) {
    AnalyticComparison result;
    const auto d = initial.size();
    if (!d || numerical.size() != time.size() * d)
        throw std::invalid_argument("Analytic comparison shape mismatch");
    if (selection == "off") {
        result.status = "disabled";
        result.reason = "Analytic comparison disabled by solve settings";
        return result;
    }
    AnalyticTrajectory trajectory;
    bool found = false;
    for (const auto &family : families) {
        if (selection != "auto" && selection != family.descriptor.identifier)
            continue;
        found = true;
        AnalyticTrajectory candidate;
        try {
            candidate = family.instantiate(initial, t0);
        } catch (const std::exception &e) {
            candidate.reason = std::string("Evaluator initialization failed: ") + e.what();
        }
        if (candidate) {
            trajectory = std::move(candidate);
            result.descriptor = family.descriptor;
            break;
        }
        if (!result.reason.empty())
            result.reason += "; ";
        result.reason += family.descriptor.identifier + ": " + candidate.reason;
    }
    if (!found) {
        result.reason = families.empty()
                            ? "This system exposes no analytic solution family"
                            : "Requested analytic family is not a property of this system";
    }
    if (!trajectory)
        return result;
    result.values.assign(numerical.size(), std::numeric_limits<double>::quiet_NaN());
    result.valid.assign(time.size(), 0);
    result.maximumAbsoluteError.assign(d, 0);
    result.status = "complete";
    result.reason.clear();
    bool failed = false;
    for (std::size_t i = 0; i < time.size(); ++i) {
        if (failed)
            continue;
        try {
            const auto x = trajectory.evaluate(time[i]);
            if (x.size() != d)
                throw std::runtime_error("Analytic evaluator changed state dimension");
            for (auto q : x)
                if (!std::isfinite(q))
                    throw std::domain_error("Analytic solution is outside its finite domain");
            result.valid[i] = 1;
            for (std::size_t j = 0; j < d; ++j) {
                result.values[i * d + j] = x[j];
                result.maximumAbsoluteError[j] =
                    std::max(result.maximumAbsoluteError[j], std::abs(x[j] - numerical[i * d + j]));
            }
        } catch (const std::exception &e) {
            result.status = "partial";
            result.reason = e.what();
            failed = true;
        }
    }
    return result;
}
} // namespace dynsyslab::core
