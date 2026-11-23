#pragma once
#include <array>
#include <vector>

namespace dynsyslab::core::special {
// All elliptic APIs use the parameter m=k^2, never the modulus k.
struct JacobiValues {
    double sn, cn, dn;
};
JacobiValues jacobi(double u, double m);
double ellipticK(double m);
double ellipticF(double phi, double m);
double jacobiAmplitude(double u, double m);
double theta(unsigned index, double z, double nome);
// Real Euler-integral branch: c>a>0, x_j<=1, and integrable combined endpoints.
// Negative x_j are supported; complex branches and analytic continuation across x_j>1 are not.
double lauricellaFD(double a, const std::vector<double>& b, double c,
                     const std::vector<double>& x);
struct WeierstrassValues {
    double value, derivative, realPeriod;
};
WeierstrassValues weierstrass(double z, double g2, double g3);
// Exact quartic-flow representation, evaluated with elliptic functions.
// Coefficients are ascending: f(x)=sum coefficients[j]*x^j; x'(0)=velocity.
class QuarticOrbit {
  public:
    QuarticOrbit(std::array<double, 5> coefficients, double initial, double velocity);
    std::array<double, 2> evaluate(double elapsed) const;
    double realPeriod() const;

  private:
    std::array<double, 5> coefficients_;
    double x_, v_, f1_, f2_, f3_, f4_, g2_, g3_, period_;
    // Stable Jacobi chart for even confining quartics crossing zero. Avoids
    // cancellation of the general Weierstrass quotient near opposite turns.
    bool evenJacobi_{};
    double amplitude_{},frequency_{},parameter_{},phase_{};
};
// Evaluation of the matrix exponential, not a numerical ODE solve.
std::vector<double> affineExponential(const std::vector<double> &matrix,
                                      const std::vector<double> &constant,
                                      const std::vector<double> &initial, double elapsed);
} // namespace dynsyslab::core::special
