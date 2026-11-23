#include <dynsyslab/core/LieAlgebra.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core::lie;
namespace {
void check(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
Polynomial p(std::initializer_list<Term> t){return {3,t};}
Polynomial sum(Polynomial a,const Polynomial& b,const std::string& factor="1") {
  const auto scaled=multiply(b,p({{{0,0,0},factor}}));
  a.terms.insert(a.terms.end(),scaled.terms.begin(),scaled.terms.end());
  return multiply(a,p({{{0,0,0},"1"}}));
}
bool equal(const Polynomial& a,const Polynomial& b){return sum(a,b,"-1").terms.empty();}
void run() {
  const auto zero=p({}),one=p({{{0,0,0},"1"}}),s=p({{{1,0,0},"1"}}),
             c=p({{{0,1,0},"1"}}),v=p({{{0,0,1},"1"}});
  const auto radius=sum(multiply(s,s),multiply(c,c));
  VectorField f{multiply(c,v),multiply(sum(zero,s,"-1"),v),
                p({{{1,0,0},"-196133/20000"},{{0,0,1},"-1/50"}})};
  VectorField e{zero,zero,one},r{c,sum(zero,s,"-1"),zero};
  check(generatedAlgebra({f,e}).status=="proven_infinite","Damped torque-input certificate");
  check(generatedAlgebra({f,r,e}).status=="proven_infinite","Full physical gradient certificate");
  auto undamped=f;undamped[2]=p({{{1,0,0},"-1"}});
  check(generatedAlgebra({undamped,e}).status=="proven_infinite","Undamped certificate");
  auto shifted=f;shifted[2]=p({{{1,0,0},"2"},{{0,1,0},"-3"},{{0,0,0},"4"},{{0,0,1},"5"}});
  check(generatedAlgebra({shifted,e}).status=="proven_infinite","Phase-shifted restoring force and bias");
  check(generatedAlgebra({f}).status=="closed","Do not infer input closure from drift alone");
  auto rotor=f;rotor[2]=p({{{0,0,1},"-1/50"}});
  check(generatedAlgebra({rotor,e}).status=="closed","Nonzero sinusoid hypothesis matters");
  auto quadraticDamping=f;quadraticDamping[2]=sum(f[2],multiply(v,v));
  check(generatedAlgebra({quadraticDamping,e},Limits{2,10,8}).status=="limit_reached",
        "Unsupported velocity dependence must not receive this certificate");
  check(generatedAlgebra({f,e},Limits{1,10,8}).status=="limit_reached","Budget status preserved");
  const auto gradient=invariantModule({one,s,c,v},{r,e},false);
  check(gradient.status=="closed","Physical derivative closure of simple trig bank");
  check(invariantModule({one,s,c,v},{f,r,e},false).status=="not_invariant","Drift breaks same bank");
  check(invariantModule({one},{f,r,e},false).status=="closed","Constants exception");
  check(lieDerivative(f,radius).terms.empty()&&lieDerivative(e,radius).terms.empty(),"Physical surface invariant");

  // Independent polynomial checks of the analytic all-n recurrence.
  // S_n=Im(c+i*s)^n, C_n=Re(c+i*s)^n. No sampled rank test is used.
  Polynomial previous=zero,sn=s,cn=c;
  const VectorField y1{zero,zero,s};
  for(unsigned n=1;n<=12;++n) {
    const auto next=sum(multiply(c,sn),multiply(s,cn));
    const auto nextCos=sum(multiply(c,cn),multiply(s,sn),"-1");
    const auto twice=bracket(y1,bracket(f,{zero,zero,sn}));
    check(twice[0].terms.empty()&&twice[1].terms.empty(),"Vertical harmonic bracket");
    auto lhs=sum(zero,twice[2],"2");
    auto rhs=sum(sum(zero,next,std::to_string(n+1)),multiply(radius,previous),"-"+std::to_string(n-1));
    check(equal(lhs,rhs),"Harmonic recurrence with radius factor");
    // On rho=1, these are precisely sin(n*q), including near unstable equilibria.
    previous=sn;sn=next;cn=nextCos;
  }
}
}
int main() {
  try{run();std::cout<<"Sinusoidal closure certificates passed\n";return 0;}
  catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
