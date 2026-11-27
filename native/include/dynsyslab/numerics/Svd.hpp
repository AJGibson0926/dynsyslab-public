#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::numerics {

// The dense SVD algorithm (per thread): divide and conquer (Gu & Eisenstat, SIAM J. Matrix Anal.
// Appl. 16:79, 1995; LAPACK gesdd's method, the default) or one-sided Jacobi (Demmel & Veselic,
// SIAM J. Matrix Anal. Appl. 13:1204, 1992; slower, with high relative accuracy for small singular
// values, the choice the gesvd driver stands for). Names: "gesdd" / "divide_and_conquer",
// "gesvd" / "jacobi".
enum class SvdDriver { divideAndConquer, jacobi };
void setSvdDriver(SvdDriver driver);
[[nodiscard]] SvdDriver svdDriver();
[[nodiscard]] SvdDriver svdDriverFromName(const std::string& name);
[[nodiscard]] const char* svdDriverName(SvdDriver driver);
class SvdDriverScope {
public:
    explicit SvdDriverScope(SvdDriver driver) : previous_(svdDriver()) { setSvdDriver(driver); }
    ~SvdDriverScope() { setSvdDriver(previous_); }
    SvdDriverScope(const SvdDriverScope&) = delete;
    SvdDriverScope& operator=(const SvdDriverScope&) = delete;
private:
    SvdDriver previous_;
};

struct RealSvd {
    std::vector<double> left;             // rows x rank
    std::vector<double> singularValues;   // rank
    std::vector<double> rightTranspose;   // rank x columns
    std::size_t rows{};
    std::size_t columns{};
    std::size_t rank{};
};

[[nodiscard]] RealSvd singularValueDecomposition(
    const std::vector<double>& matrix,std::size_t rows,std::size_t columns,
    std::size_t requestedRank=0,double relativeTolerance=0.0);

// Complex SVD M = U diag(s) V^H (thin): singular values descending; right vectors
// columns x min(rows, columns), row-major, singular vector k in column k.
struct ComplexSvd {
    std::vector<double> singularValues;
    std::vector<std::complex<double>> right;
    std::size_t rows{}, columns{}, rank{};
};
[[nodiscard]] ComplexSvd complexSingularValueDecomposition(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,std::size_t columns);

[[nodiscard]] std::vector<double> singularValues(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,
    std::size_t columns);

}  // namespace dynsyslab::numerics
