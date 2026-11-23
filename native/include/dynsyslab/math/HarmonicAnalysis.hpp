#pragma once
#include <complex>
#include <map>
#include <vector>
namespace dynsyslab::math {
struct FourierSeries {
    double angularFrequency{1};
    std::map<int,std::complex<double>> coefficients;
};
// phase is the dimensionless angular coordinate, not elapsed time.
std::complex<double> evaluate(const FourierSeries& series,double phase);
std::vector<std::complex<double>> evolveFourierSeries(const FourierSeries& series,
    const std::vector<double>& times,double initialPhase=0,double initialTime=0);
FourierSeries timeDerivative(const FourierSeries& series,unsigned order=1);
FourierSeries periodicPrimitive(const FourierSeries& series);
FourierSeries multiply(const FourierSeries& a,const FourierSeries& b);
FourierSeries projectCyclicSector(const FourierSeries& series,unsigned degree,int residue);
double meanSquare(const FourierSeries& series);
double omittedUniformBound(const FourierSeries& series,int maximumHarmonic);
// Samples cover one physical period without repeating its endpoint. Requires
// 2*maximumHarmonic < samples.size(), avoiding the Nyquist sign ambiguity.
// Reuses the native arbitrary-length FFT, normalized to exp(i*n*phase).
FourierSeries periodicFourierSeries(const std::vector<std::complex<double>>& samples,
                                   double angularFrequency,int maximumHarmonic);
// Inverse FFT at phases 2*pi*j/count. Rejects unresolved/Nyquist harmonics.
std::vector<std::complex<double>> uniformFourierSamples(const FourierSeries&,std::size_t count);
// Evaluate the supplied finite series at phi0+2*pi*j/count. Unlike the
// reconstruction API above, harmonics may exceed Nyquist: they are folded
// exactly for point evaluation, not discarded or interpreted as new coefficients.
std::vector<std::complex<double>> periodicGridValues(const FourierSeries&,std::size_t count,double phaseOrigin=0);
}
