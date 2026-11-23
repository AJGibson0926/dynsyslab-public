#include "Exact.hpp"

namespace dynsyslab::core::lie::detail {
// Exact circle lift: (s,c,v)=(sin(q),cos(q),q_dot), rho=s^2+c^2=1.
// F=(c*v,-s*v,A*s+B*c+C+D*v), E=d_v, R=c*d_s-s*d_c.
// [E,F]=R+D*E. Angular derivatives of F give Y_s and Y_c when
// A^2+B^2 != 0. For Y_h=h(q)*d_v:
// [Y_h,[F,Y_k]]=Y_(d_q(h*k)).
// Hence Y_(n+1)=(2[Y_s,[F,Y_n]]+(n-1)Y_(n-1))/(n+1),
// Y_n=sin(n*q)*d_v, n>=1, Y_0=0. The fields are independent on every
// open physical patch, not merely off the embedded circle.
// Project proof: docs/case-studies/PENDULUM_DERIVATIVE_CLOSURE.md.
// Finite embedding obstruction: Katayama/Susuki, arXiv:2611.20438, Thm 1.
bool sinusoidalMechanicalCascade(const std::vector<F>& generators,const Span& basis) {
  F input(3);input[2][{0,0,0}]=1;
  if(!basis.reduce(flatten(input)).second.empty())return false;
  const P radius{{{2,0,0},Q(1)},{{0,2,0},Q(1)}};
  // All supplied directions must restrict to the physical circle.
  // This conservative recognizer requires exact preservation of every radius.
  for(const auto& g:generators)
    if(g.size()!=3||!action(g,radius).empty())return false;
  for(const auto& f:generators) {
    if(f[0]!=P{{{0,1,1},Q(1)}}||f[1]!=P{{{1,0,1},Q(-1)}})continue;
    bool affine=true,restoring=false;
    for(const auto& [k,value]:f[2]) {
      if(k[0]+k[1]+k[2]>1)affine=false;
      if((k==Key{1,0,0}||k==Key{0,1,0})&&value!=0)restoring=true;
    }
    if(affine&&restoring)return true;
  }
  return false;
}
} // namespace dynsyslab::core::lie::detail
