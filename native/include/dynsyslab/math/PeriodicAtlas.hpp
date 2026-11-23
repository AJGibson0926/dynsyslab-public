#pragma once
#include <dynsyslab/math/HarmonicAnalysis.hpp>
namespace dynsyslab::math {
struct PeriodicAtlasKnot {
    double parameter{};
    std::vector<FourierSeries> channels;
};
// Piecewise-linear interpolation of frequency and Fourier coefficients across
// a declared orbit parameter. All knots must use a consistent phase origin.
// This is an approximate action/angle-style model, not a finite LTI KROM.
class PeriodicAtlas {
  public:
    explicit PeriodicAtlas(std::vector<PeriodicAtlasKnot>);
    std::vector<FourierSeries> interpolate(double parameter) const;
    std::vector<std::complex<double>> evaluate(double parameter,double phase) const;
    std::vector<std::complex<double>> predict(double parameter,double phase,double elapsed) const;
    const std::vector<PeriodicAtlasKnot>& knots() const { return knots_; }
  private:
    std::vector<PeriodicAtlasKnot> knots_;
};
// Cached periodic tables for repeated readout/estimation. Linear interpolation
// in parameter and phase; frequency uses the same atlas parameter interpolation.
class PeriodicAtlasTable {
  public:
    PeriodicAtlasTable(const PeriodicAtlas&,std::size_t samples=4096);
    double frequency(double parameter) const;
    std::vector<std::complex<double>> evaluate(double parameter,double phase) const;
  private:
    std::vector<double> parameters_,frequencies_;
    std::vector<std::complex<double>> values_;
    std::size_t samples_,channels_;
    std::pair<std::size_t,double> bracket(double parameter) const;
};
}
