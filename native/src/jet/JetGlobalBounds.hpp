#pragma once
#include "JetPolynomial.hpp"
#include <numeric>

namespace dynsyslab::core::jet {

// If S flips selected unrestricted axes, f(Sx)=S f(x) and h(Sx)=h(x)
// imply equality of all measured derivatives along paired trajectories.
inline json::Value signSymmetry(const std::vector<Polynomial>& f,const std::vector<Polynomial>& h,
    const std::vector<std::size_t>& positive) {
    const auto n=f.size();
    if (n>8) return {};
    for (unsigned mask=1;mask<(1u<<n);++mask) {
        auto flipped=[&](std::size_t i) {return bool(mask & (1u<<i));};
        if (std::any_of(positive.begin(),positive.end(),flipped)) continue;
        auto parityMatches=[&](const Polynomial& p,bool odd) {
            for (const auto& term:p) {
                int parity=0;
                for (std::size_t i=0;i<n;++i) if (flipped(i)) parity+=term.first[i];
                if ((parity%2!=0)!=odd) return false;
            }
            return true;
        };
        if (!std::all_of(h.begin(),h.end(),[&](const auto& p){return parityMatches(p,false);})) continue;
        bool equivariant=true;
        for (std::size_t i=0;i<n;++i) equivariant=equivariant && parityMatches(f[i],flipped(i));
        if (equivariant) {
            json::Value::Array axes;
            for (std::size_t i=0;i<n;++i) if (flipped(i)) axes.emplace_back(double(i));
            return json::Value::Object{{"method","output_preserving_sign_symmetry"},
                {"flip_axes",axes},{"state_recovery","impossible_at_every_order"}};
        }
    }
    return {};
}

// Multiply BOTH copies of f by the same nonzero monomial and clear each
// observation difference separately. Neither changes order of contact at a
// regular point. The doubled field has 2*n variables.
inline int pairedDegree(const std::vector<Polynomial>& polys, std::size_t n) {
    Powers clear(n);
    int maxTotal=0;
    for (const auto& p:polys) for (const auto& term:p) {
        for (std::size_t i=0;i<n;++i) clear[i]=std::max(clear[i],-term.first[i]);
        maxTotal=std::max(maxTotal,std::accumulate(term.first.begin(),term.first.end(),0));
    }
    return maxTotal+2*std::accumulate(clear.begin(),clear.end(),0);
}

inline json::Value finiteOrderBound(const std::vector<Polynomial>& f,const std::vector<Polynomial>& h,
    const std::set<std::size_t>& poles) {
    const auto n=f.size();
    const auto m=int(2*n);
    int delta=std::max(1,pairedDegree(f,n));
    int d=std::max(m,pairedDegree(h,n));
    auto base=[m](int d,int delta){return std::int64_t(d)+std::int64_t(m-1)*(delta-1);};
    json::Value::Array reciprocal;
    // u=1/x_i is a global chart on both signs of an excluded pole.
    // The transformed i-th rate is -u^2 f_i. Degrees alone suffice here.
    for (auto axis:poles) {
        auto transformed=[axis](std::vector<Polynomial> rows,bool field) {
            for (std::size_t i=0;i<rows.size();++i) {
                Polynomial p;
                for (const auto& [powers,c]:rows[i]) {
                    auto key=powers;
                    const bool rate=field && i==axis;
                    key[axis]=(rate?2:0)-key[axis];
                    addTerm(p,key,rate?-c:c);
                }
                rows[i]=std::move(p);
            }
            return rows;
        };
        const int candidateDelta=std::max(1,pairedDegree(transformed(f,true),n));
        const int candidateD=std::max(m,pairedDegree(transformed(h,false),n));
        if (base(candidateD,candidateDelta)<base(d,delta)) {
            d=candidateD; delta=candidateDelta; reciprocal={json::Value(double(axis))};
        }
    }
    boost::multiprecision::cpp_int bound=boost::multiprecision::cpp_int(1)<<(m+1);
    for (int i=0;i<m;++i) bound*=base(d,delta);
    return json::Value::Object{{"theorem","Binyamini (2016), Corollary 1"},
        {"reciprocal_chart_axes",std::move(reciprocal)},
        {"pair_dimension",double(m)},{"field_degree_bound",double(delta)},
        {"output_degree_bound",double(d)},{"sufficient_order",bound.str()}};
}

// Recognize by coefficient identities, independent of catalog name. The
// convergence argument applies on the whole positive-radius component.
inline json::Value volumeRecovery(const std::vector<Polynomial>& f,const std::vector<Polynomial>& h,
    const std::set<std::size_t>& poles) {
    if (f.size()!=2 || h.size()!=1 || poles.count(1) || f[0]!=Polynomial{{{0,1},1}}) return {};
    Polynomial rest=f[1];
    auto take=[&](Powers key) {Rational c=rest[key]; rest.erase(key); return c;};
    const Rational g=take({-5,0}), b=take({-1,0}), c=-take({-2,0}), d=-take({-2,1});
    const Rational inertia=take({-1,2});
    if (!rest.empty() || g<=0 || b>=0 || c<0 || d<=0 || inertia!=-Rational(3)/2) return {};
    const Polynomial volume{{{3,0},1}};
    const auto acceleration=lie(lie(volume,f),f);
    if (h.front().empty()) return {};
    const auto& lead=*acceleration.begin();
    const auto found=h.front().find(lead.first);
    if (found==h.front().end()) return {};
    const Rational scale=found->second/lead.second;
    auto residual=h.front();
    addScaled(residual,acceleration,-scale);
    if (!residual.empty()) return {};
    const Polynomial energy{{{3,2},Rational(1)/2},{{-1,0},g},{{3,0},-b/3},{{2,0},c/2}};
    const Polynomial dissipation{{{1,2},-d}};
    if (lie(energy,f)!=dissipation) return {};
    return json::Value::Object{{"method","convergence_and_second_derivative_of_volume"},
        {"domain","x0 > 0, x1 real"},{"finite_injective_order_exists",true},
        {"observable_to_volume_acceleration_scale",scale.str()},
        {"energy",encode(energy)},{"energy_derivative",encode(dissipation)}};
}
} // namespace dynsyslab::core::jet
