#include <dynsyslab/core/LieAlgebra.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core::lie;
namespace {
void check(bool ok, const char *why) {
  if (!ok)
    throw std::runtime_error(why);
}
Polynomial p(std::initializer_list<Term> t, std::size_t n = 2) {
  return {n, t};
}
bool equal(const Polynomial &a, const Polynomial &b) {
  // Canonicalization by multiplication with one.
  auto x = multiply(
      a, Polynomial{a.dimension, {{std::vector<unsigned>(a.dimension), "1"}}});
  auto y = multiply(
      b, Polynomial{b.dimension, {{std::vector<unsigned>(b.dimension), "1"}}});
  if (x.terms.size() != y.terms.size())
    return false;
  for (std::size_t i = 0; i < x.terms.size(); ++i)
    if (x.terms[i].powers != y.terms[i].powers ||
        x.terms[i].coefficient != y.terms[i].coefficient)
      return false;
  return true;
}
bool zero(const VectorField &f) {
  for (const auto &a : f)
    if (!a.terms.empty())
      return false;
  return true;
}
template <class F> void rejects(F f) {
  bool caught = false;
  try {
    f();
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  check(caught, "Invalid input not rejected");
}
void run() {
  auto one = p({{{0, 0}, "1"}}), x = p({{{1, 0}, "1"}}), y = p({{{0, 1}, "1"}}),
       x2 = p({{{2, 0}, "1"}}), z = p({});
  VectorField e1{one, z}, e2{z, one}, g{z, x};
  VectorField f{p({{{1, 0}, "-1"}}), p({{{0, 1}, "-3"}, {{2, 0}, "1"}})};
  check(equal(bracket(e1, g)[1], one), "Heisenberg sign");
  check(zero(bracket(e1, e2)), "Translations commute");
  auto h = generatedAlgebra({e1, e2, g});
  check(h.status == "closed" && h.basis.size() == 3 && h.nilpotent &&
            h.solvable && h.center.size() == 1,
        "Heisenberg structure");
  check(h.lowerCentralDimensions == std::vector<std::size_t>({3, 1, 0}),
        "Nilpotence class");
  auto t = generatedAlgebra({f, e1, e2});
  check(t.status == "closed" && t.basis.size() == 4 && t.solvable &&
            !t.nilpotent && t.center.empty(),
        "Triangular algebra");
  check(t.derivedDimensions == std::vector<std::size_t>({4, 3, 1, 0}),
        "Derived series");
  auto ideal = generatedIdeal(t, {e1});
  check(ideal.status == "closed" && ideal.basis.size() == 3 && ideal.nilpotent,
        "Constructed Heisenberg ideal");
  auto abelian = generatedAlgebra({e1, e2, e1});
  check(abelian.basis.size() == 2 && abelian.center.size() == 2 &&
            abelian.nilpotent,
        "Dependent generators");
  auto sl2 = generatedAlgebra(
      {{p({{{0}, "1"}}, 1)}, {p({{{1}, "1"}}, 1)}, {p({{{2}, "1"}}, 1)}});
  check(sl2.status == "closed" && sl2.basis.size() == 3 && !sl2.solvable &&
            !sl2.nilpotent && sl2.center.empty(),
        "sl2 structure");
  check(sl2.killingForm == std::vector<std::string>(
                               {"0", "0", "-4", "0", "2", "0", "-4", "0", "0"}),
        "sl2 Killing form");
  check(rankAt({e1, g}, {"0", "0"}) == 1 && rankAt({e1, g}, {"1", "0"}) == 2,
        "Point rank varies");
  auto eta = p({{{0, 1}, "1"}, {{2, 0}, "-1"}});
  auto m = invariantModule({one, x, x2, eta}, {f, e1, e2});
  check(m.status == "closed" && m.basis.size() == 4 &&
            m.linearStateReconstruction,
        "Triangular module");
  check(m.matrices[0] == std::vector<std::string>(
                             {"0", "0", "0", "0", "0", "-1", "0", "0", "0", "0",
                              "-2", "0", "0", "0", "0", "-3"}),
        "Exact eigenvalues");
  auto minimal = invariantModule({one, x, eta}, {f, e1, e2});
  check(minimal.basis.size() == 3 && !minimal.linearStateReconstruction,
        "Nonlinear readout not falsely rejected");
  auto completed = invariantModule({one, x, y}, {f, e1, e2});
  check(completed.status == "closed" && completed.basis.size() == 4,
        "Constructive module extension");
  auto incomplete = invariantModule({one, x, y}, {f, e1, e2}, false);
  check(incomplete.status == "not_invariant" &&
            !incomplete.residuals[0][2].terms.empty(),
        "Residual witness");
  auto matrix = m.matrices[0];
  matrix[5] = "-99/100";
  auto defects = representationDefect(m.basis, f, matrix);
  check(equal(defects[1], p({{{1, 0}, "-1/100"}})),
        "Exact perturbation defect");
  auto bound = residualBoxBound(defects, {"2", "3"});
  check(bound.vanishesAtOrigin && bound.uniformL1 == "1/50" &&
            bound.lipschitzL1 == "1/100",
        "Certified proportional error sector");
  check(!residualBoxBound({one}, {"1", "1"}).vanishesAtOrigin,
        "Nonzero bias cannot have vanishing sector");
  rejects([&] { residualBoxBound({x}, {"0", "1"}); });
  VectorField duffing{y,
                      p({{{1, 0}, "-1"}, {{3, 0}, "-1"}, {{0, 1}, "-1/10"}})};
  VectorField vdp{y, p({{{1, 0}, "-1"}, {{0, 1}, "1"}, {{2, 1}, "-1"}})};
  check(generatedAlgebra({duffing, e2}).status == "proven_infinite",
        "Duffing force cascade");
  check(generatedAlgebra({vdp, e1, e2}).status == "proven_infinite",
        "Van der Pol gradient cascade");
  VectorField radial{p({{{0,1,0},"1"}},3),
                     p({{{0,0,5},"10"},{{0,0,1},"-9"},
                        {{0,1,2},"-1/10"},{{0,2,1},"-3/2"}},3),
                     p({{{0,1,2},"-1"}},3)};
  VectorField pressure{p({},3),p({{{0,0,1},"-9"}},3),p({},3)};
  check(generatedAlgebra({radial,pressure}).status == "proven_infinite",
        "Pressure-controlled inverse-radius cascade");
  radial[1]=p({{{0,2,1},"-1"}},3);
  check(generatedAlgebra({radial,pressure}).status == "closed",
        "The exceptional linearizable inertia coefficient must not be rejected");
  auto limited = generatedAlgebra({duffing, e2}, Limits{2, 100, 8});
  check(limited.status == "limit_reached",
        "Budget exhaustion is not nonexistence");
  VectorField oscillator{eta,
                         p({{{1, 0}, "-1"}, {{1, 1}, "2"}, {{3, 0}, "-2"}})};
  auto osc = invariantModule({one, x, eta}, {oscillator, e1, e2});
  check(osc.status == "closed" && osc.basis.size() == 3,
        "Nonlinear oscillator module");
  auto osc6 = invariantModule(
      {one, x, eta, multiply(x, x), multiply(x, eta), multiply(eta, eta)},
      {oscillator, e1, e2});
  check(osc6.status == "closed" && osc6.basis.size() == 6 &&
            osc6.linearStateReconstruction,
        "Quadratic oscillator lift");
  check(generatedAlgebra({oscillator, e1, e2}).status == "closed",
        "Oscillator algebra");
  check(equal(poissonBracket(x, eta, {"0", "1", "-1", "0"}), one),
        "Canonical shear preserves Poisson structure");
  check(equal(flowBracket(f, x, x2), p({{{3, 0}, "-1"}})),
        "Flow Jacobi bracket eigenvalue sum");
  auto rationalPoly = p({{{2, 0}, "2/3"}, {{2, 0}, "1/3"}});
  check(equal(differentiate(rationalPoly, 0), p({{{1, 0}, "2"}})),
        "Rational duplicate coefficients");
  check(equal(multiply(x, p({{{0, 0}, "0"}})), z), "Zero product");
  rejects([&] { differentiate(x, 2); });
  rejects([&] { multiply(x, p({}, 1)); });
  rejects([&] { multiply(p({{{0, 0}, "1/0"}}), x); });
  rejects([&] { poissonBracket(x, y, {"0", "1", "1", "0"}); });
  rejects([&] { invariantModule({x, x}, {f}); });
  rejects([&] { generatedAlgebra({f}, Limits{0, 1, 1}); });
  auto j = dynsyslab::core::json::parse(
      R"({"schema":"dynsyslab.lie-algebra-query/v1","dimension":2,"action":"derivative","axis":0,"observable":[{"powers":[2,0],"coefficient":"1/3"}]})");
  check(analyze(j).at("value").array()[0].at("coefficient").string() == "2/3",
        "JSON boundary");
  j.at("dimension") = -1.;
  rejects([&] { analyze(j); });
}
} // namespace
int main() {
  try {
    run();
    std::cout << "Lie algebra exact computation tests passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
