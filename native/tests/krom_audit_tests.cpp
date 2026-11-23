#include <cmath>
#include <dynsyslab/core/KromDifferentialAnalysis.hpp>
#include <dynsyslab/core/KromLocalControlAnalysis.hpp>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void check(bool v, const char *reason) {
  if (!v)
    throw std::runtime_error(reason);
}
void near(double a, double b) {
  check(std::abs(a - b) < 1e-10 * (1 + std::abs(b)),
        "Numerical identity failed");
}
} // namespace
int main() {
  try {
    SystemConfiguration c;
    c.identifier = "quadratic_triangular";
    c.parameters = {{"a", -1},
                    {"b", -2},
                    {"c", 0},
                    {"input_gain_1", 1},
                    {"input_gain_2", 1}};
    auto model = KoopmanModel::fromPolynomialCoefficients(
        c, 2, {1, 0, 0, 1}, {{1, 0}, {0, 0}, {0, 0}, {1, 0}},
        {{-1, 0}, {-2, 0}});
    model.setControlMap(
        "constant", ConstantControlMap{{1., 1.}, 2, 1, "exact linear test"});
    model.setReconstruction("x_1", {{1., 0.}, 1, 2});
    const auto d = analyzeKromDifferentials(model, {1, 2, -.4, .7}, "constant");
    check(d.samples == 2 && d.jacobians.size() == 8, "Differential shapes");
    for (double v : d.generatorDefect)
      near(v, 0);
    for (double v : d.inputDefect)
      near(v, 0);
    for (double v : d.smallestSingularValue)
      near(v, 1);
    for (double v : d.orientedDeterminant)
      near(v, 1);
    KromLocalControlOptions o;
    o.point = {0, 0};
    o.measurement = "x_1";
    o.statePenalty = 1;
    o.inputPenalty = 1;
    model.setReconstruction("state",
                            {{1., 0., 0., 1.}, 2, 2, 0, 0, 0, {.2, .3}});
    const auto local = analyzeKromLocalControl(model, o);
    near(local.decodedStateAtPoint[0], .2);
    near(local.decodedStateAtPoint[1], .3);
    near(local.generatorAtPoint[0], 0);
    near(local.measuredAtPoint[0], 0);
    check(local.physicalLqr.closedLoopStable &&
              local.physicalLqg.closedLoopStable,
          "Exact model lost stability");
    near(local.physicalLqg.minimumReturnDifference,
         local.surrogateLqg.minimumReturnDifference);
    near(local.physicalLqr.minimumReturnDifference,
         local.surrogateLqr.minimumReturnDifference);
    for (std::size_t i = 0; i < 2; ++i)
      near(local.gain[i], local.physicalGain[i]);
    const auto regression =
        regressInvariantAction({0, 1, 2}, {1, 3, 5}, {3, 4}, {7, 9}, 1, 1);
    near(regression.coefficients[0], 1);
    near(regression.coefficients[1], 2);
    near(regression.validationRelativeError, 0);
    const auto nonclosed =
        regressInvariantAction({-1, 0, 1}, {1, 0, 1}, {2, 3}, {4, 9}, 1, 1);
    check(nonclosed.validationRelativeError > .5,
          "Nonlinear action misreported closed");
    model.setControlMap(
        "constant", ConstantControlMap{{2., 3.}, 2, 1, "deliberately wrong"});
    const auto bad = analyzeKromDifferentials(model, {.2, .3}, "constant");
    near(bad.inputDefect[0], -1);
    near(bad.inputDefect[1], -2);
    bool rejected = false;
    try {
      (void)analyzeKromDifferentials(model, {NAN, 0});
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    check(rejected, "NaN state accepted");
    std::cout << "krom_audit_tests: PASS\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
