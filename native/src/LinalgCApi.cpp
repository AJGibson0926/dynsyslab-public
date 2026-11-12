#include <dynsyslab/linalg_c_api.h>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string lastError;

dsl_status fail(dsl_status status,const std::string& message) {
    lastError=message; return status;
}

dsl_status copyText(const std::string& value,char* buffer,size_t capacity,size_t* required) {
    const size_t needed=value.size()+1;
    if (required) *required=needed;
    if (!buffer||capacity<needed) return DSL_STATUS_BUFFER_TOO_SMALL;
    std::memcpy(buffer,value.c_str(),needed); return DSL_STATUS_OK;
}

template<class Function> dsl_status guarded(Function&& operation) {
    try { operation(); lastError.clear(); return DSL_STATUS_OK; }
    catch (const std::invalid_argument& error) {
        return fail(DSL_STATUS_INVALID_ARGUMENT,error.what());
    } catch (const std::exception& error) {
        return fail(DSL_STATUS_RUNTIME_ERROR,error.what());
    } catch (...) { return fail(DSL_STATUS_RUNTIME_ERROR,"unknown linear-algebra exception"); }
}

void validateShape(const void* matrix,size_t rows,size_t columns) {
    if (!matrix||!rows||!columns||rows>(std::numeric_limits<size_t>::max)()/columns)
        throw std::invalid_argument("invalid linear-algebra matrix shape");
}

}  // namespace

extern "C" {

uint32_t dsl_linalg_abi_version(void) { return DSL_LINALG_ABI_VERSION; }
dsl_status dsl_linalg_last_error(char* buffer,size_t capacity,size_t* required) {
    return copyText(lastError,buffer,capacity,required);
}

dsl_status dsl_linalg_svd_real(
    const double* matrix,size_t rows,size_t columns,size_t requestedRank,
    double* left,size_t leftCapacity,double* singular,size_t singularCapacity,
    double* right,size_t rightCapacity,size_t* outputRank) {
    return guarded([&] {
        validateShape(matrix,rows,columns);
        if (!left||!singular||!right||!outputRank)
            throw std::invalid_argument("SVD output pointer is null");
        const auto result=dynsyslab::numerics::singularValueDecomposition(
            {matrix,matrix+rows*columns},rows,columns,requestedRank);
        if (leftCapacity<rows*result.rank||singularCapacity<result.rank
            ||rightCapacity<result.rank*columns)
            throw std::invalid_argument("SVD output capacity is too small");
        std::copy(result.left.begin(),result.left.end(),left);
        std::copy(result.singularValues.begin(),result.singularValues.end(),singular);
        std::copy(result.rightTranspose.begin(),result.rightTranspose.end(),right);
        *outputRank=result.rank;
    });
}

dsl_status dsl_linalg_least_squares_real(
    const double* matrix,size_t rows,size_t columns,const double* rhs,size_t outputs,
    double ridge,double* solution,size_t capacity) {
    return guarded([&] {
        if (!solution||capacity<columns*outputs)
            throw std::invalid_argument("real least-squares output capacity is too small");
        validateShape(matrix,rows,columns); validateShape(rhs,rows,outputs);
        dynsyslab::numerics::LeastSquaresOptions options;options.ridge=ridge;
        const std::vector<double> activeMatrix(matrix,matrix+rows*columns);
        const std::vector<double> activeRhs(rhs,rhs+rows*outputs);
        const auto value=dynsyslab::numerics::solveLeastSquares(
            activeMatrix,rows,columns,activeRhs,outputs,options);
        std::copy(value.begin(),value.end(),solution);
    });
}

dsl_status dsl_linalg_least_squares_complex(
    const dsl_complex64* matrix,size_t rows,size_t columns,
    const dsl_complex64* rhs,size_t outputs,double ridge,
    dsl_complex64* solution,size_t capacity) {
    return guarded([&] {
        if (!solution||capacity<columns*outputs)
            throw std::invalid_argument("complex least-squares output capacity is too small");
        validateShape(matrix,rows,columns); validateShape(rhs,rows,outputs);
        std::vector<std::complex<double>> activeMatrix(rows*columns);
        std::vector<std::complex<double>> activeRhs(rows*outputs);
        for (std::size_t index=0;index<activeMatrix.size();++index)
            activeMatrix[index]={matrix[index].real,matrix[index].imag};
        for (std::size_t index=0;index<activeRhs.size();++index)
            activeRhs[index]={rhs[index].real,rhs[index].imag};
        dynsyslab::numerics::LeastSquaresOptions options;options.ridge=ridge;
        const auto value=dynsyslab::numerics::solveLeastSquares(
            activeMatrix,rows,columns,activeRhs,outputs,options);
        for (std::size_t index=0;index<value.size();++index)
            solution[index]={value[index].real(),value[index].imag()};
    });
}

dsl_status dsl_linalg_eigen_real(
    const double* matrix,size_t dimension,dsl_complex64* eigenvalues,
    size_t eigenvalueCapacity,dsl_complex64* vectors,size_t vectorCapacity) {
    return guarded([&] {
        validateShape(matrix,dimension,dimension);
        if (!eigenvalues||eigenvalueCapacity<dimension||!vectors
            ||vectorCapacity<dimension*dimension)
            throw std::invalid_argument("eigensystem output capacity is too small");
        const auto result=dynsyslab::numerics::eigensystem(
            {matrix,matrix+dimension*dimension},dimension,true);
        for (size_t index=0;index<dimension;++index) {
            const auto value=result.eigenvalues[index];
            eigenvalues[index]={value.real(),value.imag()};
        }
        for (size_t row=0;row<dimension;++row)
            for (size_t column=0;column<dimension;++column) {
                const auto value=result.rightEigenvectors[row*dimension+column];
                vectors[row*dimension+column]={value.real(),value.imag()};
            }
    });
}

}  // extern "C"
