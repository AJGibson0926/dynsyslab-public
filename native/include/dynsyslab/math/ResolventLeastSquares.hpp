#pragma once
#include <dynsyslab/math/ResolventLinearAlgebra.hpp>
// Thin complex SVD by one-sided Jacobi, with singular vectors, and minimum-norm least squares.
// Hestenes, J. SIAM 6 (1958) 51-90; Golub & Van Loan, Matrix Computations (2013), sec. 8.6.
namespace dynsyslab::math {
struct ThinSvd {ComplexMatrix U,V;std::vector<double> singularValues;}; // A=U diag(s) V^H, descending.
ThinSvd thinSvd(const ComplexMatrix& a);
struct LeastSquaresSolution {ComplexMatrix x;std::size_t rank{};double condition{};};
// Minimizes ||A x - B||_F with singular values below rankTolerance*s_max dropped.
LeastSquaresSolution leastSquares(const ComplexMatrix& a,const ComplexMatrix& b,double rankTolerance=1e-12);
}
