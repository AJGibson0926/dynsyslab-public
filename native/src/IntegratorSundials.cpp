#include "IntegratorBackend.hpp"
#include <dynsyslab/core/Execution.hpp>

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_context.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

// SUNDIALS CVODE steppers (Hindmarsh et al., "SUNDIALS: Suite of nonlinear and
// differential/algebraic equation solvers", ACM TOMS 31 (2005) 363-396):
//   bdf   - variable-order BDF (orders 1-5) with Newton iteration and a dense direct solver;
//   adams - variable-order Adams-Moulton (orders 1-12), Newton iteration;
//   lsoda - automatic nonstiff/stiff switching in the spirit of Petzold's LSODA
//           (SIAM J. Sci. Stat. Comput. 4 (1983) 136-148). At every interval end the stepper
//           compares the last step h with the explicit stability limit 1/||J||_inf: Adams
//           switches to BDF after three intervals whose step is stability limited
//           (h ||J|| > 0.5), and BDF returns to Adams once h ||J|| < 0.5. The test differs
//           from ODEPACK's step-size comparison; switches are counted in methodSwitches.
// The stepper stops exactly at each interval end (CVodeSetStopTime), so held inputs never
// leak across a sampled-control tick.

namespace dynsyslab::core::detail {
namespace {

void require(int status, const char* operation) {
    if (status < 0) throw std::runtime_error(std::string("SUNDIALS ") + operation
                                             + " failed with code " + std::to_string(status));
}

struct Callback {
    RhsContext* context{};
    std::exception_ptr error;
};

int rightHandSide(sunrealtype time, N_Vector state, N_Vector derivative, void* raw) {
    auto& data = *static_cast<Callback*>(raw);
    try {
        data.context->rhs(time, N_VGetArrayPointer(state), N_VGetArrayPointer(derivative));
        return 0;
    } catch (...) { data.error = std::current_exception(); return -1; }
}

class CvodeStepper final : public Stepper {
public:
    CvodeStepper(const IntegrationOptions& options, std::size_t n) : options_(options), n_(n) {}
    ~CvodeStepper() override { release(); }

    void restart(RhsContext& context, double time, const std::vector<double>& state) override {
        if (time == options_.initialTime || !cvode_) {
            stiff_ = options_.method == IntegrationMethod::bdf;
            limitedIntervals_ = 0;
            create(context, time, state);
            return;
        }
        absorbStatistics(context);
        std::copy(state.begin(), state.end(), N_VGetArrayPointer(vector_));
        require(CVodeReInit(cvode_, time, vector_), "sampled-control reinitialization");
        lastSteps_ = lastFailures_ = 0;
    }

    void advance(RhsContext& context, double time, double end, std::vector<double>& state) override {
        callback_.context = &context;
        double* active = N_VGetArrayPointer(vector_);
        std::copy(state.begin(), state.end(), active);
        require(CVodeSetStopTime(cvode_, end), "stop-time configuration");
        sunrealtype reached{};
        const int status = CVode(cvode_, end, vector_, &reached, CV_NORMAL);
        if (callback_.error) {
            auto error = callback_.error; callback_.error = nullptr;
            std::rethrow_exception(error);
        }
        require(status, "integration");
        absorbStatistics(context);
        std::copy(active, active + n_, state.begin());
        context.accept(end, state.data());
        if (options_.method == IntegrationMethod::lsoda) switchIfNeeded(context, end, state, time);
    }

private:
    void create(RhsContext& context, double time, const std::vector<double>& state) {
        release();
        require(SUNContext_Create(SUN_COMM_NULL, &sunContext_), "context creation");
        vector_ = N_VNew_Serial(static_cast<sunindextype>(n_), sunContext_);
        if (!vector_) throw std::runtime_error("SUNDIALS state allocation failed");
        matrix_ = SUNDenseMatrix(static_cast<sunindextype>(n_), static_cast<sunindextype>(n_), sunContext_);
        if (!matrix_) throw std::runtime_error("SUNDIALS matrix allocation failed");
        solver_ = SUNLinSol_Dense(vector_, matrix_, sunContext_);
        if (!solver_) throw std::runtime_error("SUNDIALS linear-solver allocation failed");
        cvode_ = CVodeCreate(stiff_ ? CV_BDF : CV_ADAMS, sunContext_);
        if (!cvode_) throw std::runtime_error("SUNDIALS CVODE allocation failed");
        std::copy(state.begin(), state.end(), N_VGetArrayPointer(vector_));
        callback_.context = &context;
        require(CVodeInit(cvode_, rightHandSide, time, vector_), "initialization");
        require(CVodeSetUserData(cvode_, &callback_), "user-data binding");
        require(CVodeSStolerances(cvode_, options_.relativeTolerance, options_.absoluteTolerance),
                "tolerance configuration");
        require(CVodeSetLinearSolver(cvode_, solver_, matrix_), "dense linear solver binding");
        require(CVodeSetMaxNumSteps(cvode_, static_cast<long int>(
            std::min<std::size_t>(options_.maximumInternalSteps, 2000000000))), "step budget");
        if (options_.initialStep > 0.0)
            require(CVodeSetInitStep(cvode_, options_.initialStep), "initial-step configuration");
        if (options_.maximumStep > 0.0)
            require(CVodeSetMaxStep(cvode_, options_.maximumStep), "maximum-step configuration");
        lastSteps_ = lastFailures_ = 0;
    }

    void absorbStatistics(RhsContext& context) {
        long int steps = 0, failures = 0;
        require(CVodeGetNumSteps(cvode_, &steps), "step statistics");
        require(CVodeGetNumErrTestFails(cvode_, &failures), "error statistics");
        const auto newSteps = static_cast<std::size_t>(std::max(0L, steps - lastSteps_));
        const auto newFailures = static_cast<std::size_t>(std::max(0L, failures - lastFailures_));
        lastSteps_ = steps; lastFailures_ = failures;
        context.statistics.acceptedSteps += newSteps;
        context.statistics.rejectedSteps += newFailures;
        context.trajectorySteps += newSteps + newFailures;
        if (context.trajectorySteps > options_.maximumInternalSteps)
            throw std::runtime_error("native integrator exceeded its step budget (max_internal_steps)");
    }

    void switchIfNeeded(RhsContext& context, double end, const std::vector<double>& state, double start) {
        sunrealtype last{};
        if (CVodeGetLastStep(cvode_, &last) < 0 || !(last > 0.0)) return;
        double norm = 0.0;
        try {
            const auto jacobian = context.system.jacobian(end, state.data());
            if (jacobian.size() != n_ * n_) return;
            for (std::size_t r = 0; r < n_; ++r) {
                double row = 0.0;
                for (std::size_t c = 0; c < n_; ++c) row += std::abs(jacobian[r * n_ + c]);
                norm = std::max(norm, row);
            }
        } catch (const std::exception&) { return; }
        const double product = last * norm;
        // A step clipped by the interval end says nothing about stability.
        const bool clipped = last >= 0.999 * (end - start);
        if (!stiff_) {
            limitedIntervals_ = (!clipped && product > 0.5) ? limitedIntervals_ + 1 : 0;
            if (limitedIntervals_ < 3) return;
        } else if (product >= 0.5) {
            return;
        }
        absorbStatistics(context);
        stiff_ = !stiff_;
        limitedIntervals_ = 0;
        ++context.statistics.methodSwitches;
        create(context, end, state);
    }

    void release() {
        if (cvode_) CVodeFree(&cvode_);
        if (solver_) SUNLinSolFree(solver_);
        if (matrix_) SUNMatDestroy(matrix_);
        if (vector_) N_VDestroy(vector_);
        if (sunContext_) SUNContext_Free(&sunContext_);
        cvode_ = nullptr; solver_ = nullptr; matrix_ = nullptr; vector_ = nullptr; sunContext_ = nullptr;
    }

    const IntegrationOptions& options_;
    std::size_t n_;
    SUNContext sunContext_{};
    N_Vector vector_{};
    SUNMatrix matrix_{};
    SUNLinearSolver solver_{};
    void* cvode_{};
    Callback callback_;
    long int lastSteps_{}, lastFailures_{};
    bool stiff_{true};
    int limitedIntervals_{};
};

}  // namespace

std::unique_ptr<Stepper> makeCvodeStepper(const IntegrationOptions& options, std::size_t n) {
    return std::make_unique<CvodeStepper>(options, n);
}

}  // namespace dynsyslab::core::detail
