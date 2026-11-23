#pragma once
#include <complex>
#include <map>
#include <vector>
namespace dynsyslab::math {
using FourierIndex=std::vector<int>;
struct TorusFourierSeries {
    std::vector<double> angularFrequencies;
    std::map<FourierIndex,std::complex<double>> coefficients;
};
// Row-major product mesh, last angle fastest; phase_j=2*pi*i_j/shape_j.
// Finite Fourier approximation, not an assertion about unrepresented tails.
TorusFourierSeries torusFourierSeries(const std::vector<std::complex<double>>& samples,
    const std::vector<std::size_t>& shape,const std::vector<double>& frequencies,
    const std::vector<int>& maximumHarmonics);
std::complex<double> evaluate(const TorusFourierSeries&,const std::vector<double>& phases);
TorusFourierSeries torusTimeDerivative(const TorusFourierSeries&);
TorusFourierSeries torusPartialDerivative(const TorusFourierSeries&,std::size_t axis);
double harmonicFrequency(const FourierIndex&,const std::vector<double>&);
double meanSquare(const TorusFourierSeries&);
struct CohomologicalSolution {
    TorusFourierSeries solution,retained;
    double maximumInverseDivisor{};
};
// omega . grad(solution) = forcing - retained. Small/zero divisors are retained,
// never silently divided out. This is the finite Fourier homological step only.
CohomologicalSolution solveCohomological(const TorusFourierSeries&,double divisorTolerance);
}
