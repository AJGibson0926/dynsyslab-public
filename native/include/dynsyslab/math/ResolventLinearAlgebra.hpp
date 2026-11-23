#pragma once
#include <complex>
#include <cstddef>
#include <vector>
// Small Eigen-free complex dense kernels for the resolvent tools. dynsyslab_math must
// build without the optional Eigen numerics library, so these stay local to it.
// Algorithms: Golub & Van Loan, Matrix Computations (4th ed., 2013): LU with partial
// pivoting (sec. 3.4), Householder-Hessenberg and shifted QR (sec. 7.4-7.5), one-sided
// Jacobi SVD (Hestenes, J. SIAM 6 (1958); sec. 8.6).
namespace dynsyslab::math {
using ComplexVector=std::vector<std::complex<double>>;
struct ComplexMatrix {
    std::size_t rows{},cols{};
    ComplexVector data; // row-major
    ComplexMatrix()=default;
    ComplexMatrix(std::size_t r,std::size_t c):rows(r),cols(c),data(r*c){}
    std::complex<double>& operator()(std::size_t i,std::size_t j){return data[i*cols+j];}
    const std::complex<double>& operator()(std::size_t i,std::size_t j) const {return data[i*cols+j];}
    static ComplexMatrix identity(std::size_t n);
};
ComplexMatrix product(const ComplexMatrix& a,const ComplexMatrix& b);
ComplexVector product(const ComplexMatrix& a,const ComplexVector& x);
ComplexMatrix adjoint(const ComplexMatrix& a);
double frobeniusNorm(const ComplexMatrix& a);
double euclideanNorm(const ComplexVector& x);
// s I - A.
ComplexMatrix shiftedMatrix(const ComplexMatrix& a,std::complex<double> s);
// LU with partial pivoting. solve throws std::domain_error on an exactly singular pivot.
class ComplexLu {
  public:
    explicit ComplexLu(ComplexMatrix a);
    ComplexVector solve(ComplexVector b) const;
    ComplexMatrix solve(const ComplexMatrix& b) const;
    // Solves A^H x = b with the same factors.
    ComplexVector solveAdjoint(const ComplexVector& b) const;
    ComplexMatrix solveAdjoint(const ComplexMatrix& b) const;
    bool singular() const {return singular_;}
  private:
    ComplexMatrix lu_;
    std::vector<std::size_t> pivots_;
    bool singular_{false};
};
// Descending singular values.
std::vector<double> complexSingularValues(ComplexMatrix a);
// Eigenvalues by Householder-Hessenberg reduction and Wilkinson-shifted QR.
ComplexVector complexEigenvalues(ComplexMatrix a);
}
