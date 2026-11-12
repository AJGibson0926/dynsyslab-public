#include "polynomial_kernel.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
namespace {
struct Polynomial {
    std::vector<std::vector<std::pair<std::size_t,int>>> powers;
    std::vector<std::vector<std::pair<std::size_t,double>>> rows;
    std::vector<double> values;
};
}
extern "C" void* mc_poly_new_v1(std::size_t n,std::size_t outputs,std::size_t f,
                                const int* e,const double* weights) {
    const auto lim=std::numeric_limits<std::size_t>::max();
    if (!n || !outputs || !f || !e || !weights || f>lim/n || f>lim/outputs) return nullptr;
    try {
        auto p=std::make_unique<Polynomial>();
        p->powers.resize(f); p->rows.resize(outputs); p->values.resize(f);
        for (std::size_t k=0;k<f;++k) for (std::size_t j=0;j<n;++j) {
            int power=e[k*n+j];
            if (power<0) return nullptr;
            if (power) p->powers[k].emplace_back(j,power);
        }
        for (std::size_t j=0;j<outputs;++j) for(std::size_t k=0;k<f;++k) {
            double c=weights[j*f+k];
            if (!std::isfinite(c)) return nullptr;
            if (c!=0.) p->rows[j].emplace_back(k,c);
        }
        return p.release();
    } catch (...) { return nullptr; }
}
extern "C" int mc_poly_eval_v1(void* field,const double* x,double* out) {
    if (!field || !x || !out) return 1;
    auto& p=*static_cast<Polynomial*>(field);
    for (std::size_t k=0;k<p.powers.size();++k) {
        double v=1.;
        for(auto [j,degree]:p.powers[k]) v*=degree==1?x[j]:std::pow(x[j],degree);
        p.values[k]=v;
    }
    for (std::size_t j=0;j<p.rows.size();++j) {
        double v=0.;
        for(auto [k,c]:p.rows[j]) v+=c*p.values[k];
        if (!std::isfinite(v)) return 1;
        out[j]=v;
    }
    return 0;
}
extern "C" void mc_poly_drop_v1(void* field) {delete static_cast<Polynomial*>(field);}
