#include <dynsyslab/linalg_c_api.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        require(dsl_linalg_abi_version()==1,"linear-algebra ABI version");
        const double matrix[]{1,0,0,2};
        double left[4]{},singular[2]{},right[4]{}; size_t rank=0;
        require(dsl_linalg_svd_real(matrix,2,2,0,left,4,singular,2,right,4,&rank)
            ==DSL_STATUS_OK,"SVD status");
        require(rank==2&&std::abs(singular[0]-2.0)<1e-12,"SVD values");
        const double design[]{1,0,1,1,1,2}; const double observations[]{1,3,5};
        double coefficients[2]{};
        require(dsl_linalg_least_squares_real(
            design,3,2,observations,1,0.0,coefficients,2)==DSL_STATUS_OK,
            "least squares status");
        require(std::abs(coefficients[0]-1.0)<1e-12
            &&std::abs(coefficients[1]-2.0)<1e-12,"least squares values");
        const dsl_complex64 complexDesign[]{{1,1},{1,-1}};
        const dsl_complex64 complexObservations[]{{1,3},{3,-1}};
        dsl_complex64 complexCoefficient{};
        require(dsl_linalg_least_squares_complex(complexDesign,2,1,
            complexObservations,1,0.0,&complexCoefficient,1)==DSL_STATUS_OK,
            "complex least squares status");
        require(std::abs(complexCoefficient.real-2.0)<1e-12
            &&std::abs(complexCoefficient.imag-1.0)<1e-12,
            "complex least squares values");
        dsl_complex64 eigenvalues[2]{},vectors[4]{};
        const double rotation[]{0,-1,1,0};
        require(dsl_linalg_eigen_real(rotation,2,eigenvalues,2,vectors,4)
            ==DSL_STATUS_OK,"eigensystem status");
        require(std::abs(std::abs(eigenvalues[0].imag)-1.0)<1e-12,
                "eigensystem values");
        std::cout<<"dynsyslab_linalg_c_api_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        std::cerr<<"dynsyslab_linalg_c_api_tests: FAIL: "<<error.what()<<'\n'; return 1;
    }
}
