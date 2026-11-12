#include "rollout_kernel.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

extern "C" int mc_rk4_v1(mc_rhs_v1 rhs, void* ctx, std::size_t n,
    const double* initial, const double* time, std::size_t count, double maxStep,
    std::size_t budget, const double* limits, double* out, std::size_t* written,
    std::size_t* evaluations) {
    if (!written || !evaluations) return 4;
    *written = *evaluations = 0;
    if (!rhs || !initial || !time || !out || !n || !count || !budget ||
        n > std::numeric_limits<std::size_t>::max()/count ||
        !std::isfinite(maxStep) || maxStep <= 0) return 4;
    for (std::size_t i=0; i<count; ++i)
        if (!std::isfinite(time[i]) || (i && time[i]<=time[i-1])) return 4;
    for (std::size_t j=0; j<n; ++j)
        if (!std::isfinite(initial[j]) || (limits && !(limits[j]>0))) return 4;
    auto status = [&](const double* x) {
        for (std::size_t j=0; j<n; ++j) {
            if (!std::isfinite(x[j])) return 2;
            if (limits && std::abs(x[j]) >= limits[j]) return 1;
        }
        return 0;
    };
    if (auto err=status(initial)) return err;
    try {
        std::vector<double> x(initial,initial+n), stage(n), k(n), total(n);
        std::copy(x.begin(),x.end(),out);
        *written=1;
        std::size_t steps=0;
        for (std::size_t i=1; i<count; ++i) {
            const double duration=time[i]-time[i-1];
            const double required=std::max(1.,std::ceil(duration/maxStep));
            if (required>=static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
                required>static_cast<double>(budget-steps)) return 3;
            const auto sub=static_cast<std::size_t>(required);
            const double h=duration/sub;
            for (std::size_t s=0; s<sub; ++s) {
                stage=x;
                for (int pass=0; pass<4; ++pass) {
                    if (auto err=status(stage.data())) return err;
                    ++*evaluations;
                    if (rhs(ctx,stage.data(),k.data())) return 2;
                    for (std::size_t j=0; j<n; ++j) {
                        if (!std::isfinite(k[j])) return 2;
                        if (pass==0) total[j]=k[j];
                        else total[j]+=(pass==3?1.:2.)*k[j];
                        stage[j]=x[j]+(pass<2?.5:1.)*h*k[j];
                    }
                }
                for (std::size_t j=0; j<n; ++j) x[j]+=h*total[j]/6.;
                if (auto err=status(x.data())) return err;
            }
            steps+=sub;
            std::copy(x.begin(),x.end(),out+i*n);
            ++*written;
        }
        return 0;
    } catch (...) { return 5; }
}
