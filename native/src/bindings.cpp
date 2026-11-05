#include <dynsyslab/core/KROM.hpp>
#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/ODE.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Version.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace pybind11::literals;
using dynsyslab::core::FixedMethod;

namespace {

py::tuple integrateFixed(
    const std::string& model,
    py::array_t<double, py::array::c_style | py::array::forcecast> parametersArray,
    py::array_t<double, py::array::c_style | py::array::forcecast> initialArray,
    double initialTime, double finalTime, double stepSize, const std::string& method) {
    const auto parameters = parametersArray.request();
    const auto initial = initialArray.request();
    if (parameters.ndim != 1 || (initial.ndim != 1 && initial.ndim != 2))
        throw std::invalid_argument("native integration arrays have invalid dimensions");
    const auto batch = static_cast<std::size_t>(initial.ndim == 1 ? 1 : initial.shape[0]);
    const auto dimension = static_cast<std::size_t>(initial.ndim == 1 ? initial.shape[0] : initial.shape[1]);
    dynsyslab::core::FixedIntegrationRequest request;
    request.model = model;
    request.parameters.assign(static_cast<const double*>(parameters.ptr),
                              static_cast<const double*>(parameters.ptr) + parameters.shape[0]);
    request.initialStates.assign(static_cast<const double*>(initial.ptr),
                                 static_cast<const double*>(initial.ptr) + batch * dimension);
    request.batch = batch;
    request.initialTime = initialTime;
    request.finalTime = finalTime;
    request.stepSize = stepSize;
    if (method == "Euler") request.method = FixedMethod::euler;
    else if (method == "RK4") request.method = FixedMethod::rk4;
    else throw std::invalid_argument("native fixed integration supports Euler or RK4");
    dynsyslab::core::TrajectoryBatch result;
    {
        py::gil_scoped_release release;
        result = dynsyslab::core::integrateFixed(request);
    }
    py::array_t<double> times(result.time.size());
    std::copy(result.time.begin(), result.time.end(), static_cast<double*>(times.request().ptr));
    if (initial.ndim == 1) {
        py::array_t<double> values({static_cast<py::ssize_t>(result.samples),
                                    static_cast<py::ssize_t>(result.stateDimension)});
        std::copy(result.values.begin(), result.values.end(), static_cast<double*>(values.request().ptr));
        return py::make_tuple(times, values);
    }
    py::array_t<double> values({static_cast<py::ssize_t>(result.batch),
                                static_cast<py::ssize_t>(result.samples),
                                static_cast<py::ssize_t>(result.stateDimension)});
    std::copy(result.values.begin(), result.values.end(), static_cast<double*>(values.request().ptr));
    return py::make_tuple(times, values);
}

py::tuple integrateSystem(
    const std::string& identifier, const std::string& preset,
    const std::map<std::string, double>& parameters,
    py::array_t<double, py::array::c_style | py::array::forcecast> initialArray,
    double initialTime, double finalTime, double outputStep, const std::string& method,
    double relativeTolerance, double absoluteTolerance) {
    const auto initial = initialArray.request();
    if (initial.ndim != 1 && initial.ndim != 2)
        throw std::invalid_argument("native system initial state must be one- or two-dimensional");
    const auto batch = static_cast<std::size_t>(initial.ndim == 1 ? 1 : initial.shape[0]);
    const auto dimension = static_cast<std::size_t>(initial.ndim == 1 ? initial.shape[0] : initial.shape[1]);
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = identifier;
    configuration.preset = preset;
    configuration.parameters = parameters;
    auto system = dynsyslab::core::createSystem(configuration);
    if (dimension != system->descriptor().stateDimension)
        throw std::invalid_argument("native system initial-state dimension mismatch");
    dynsyslab::core::IntegrationOptions options;
    options.initialTime = initialTime;
    options.finalTime = finalTime;
    options.outputStep = outputStep;
    options.relativeTolerance = relativeTolerance;
    options.absoluteTolerance = absoluteTolerance;
    if (method == "Euler") options.method = dynsyslab::core::IntegrationMethod::euler;
    else if (method == "RK4") options.method = dynsyslab::core::IntegrationMethod::rk4;
    else if (method == "RK45") options.method = dynsyslab::core::IntegrationMethod::dormandPrince45;
    else throw std::invalid_argument("native system integration supports Euler, RK4, or RK45");
    std::vector<double> initialStates(static_cast<const double*>(initial.ptr),
        static_cast<const double*>(initial.ptr) + batch * dimension);
    dynsyslab::core::NativeSolution result;
    {
        py::gil_scoped_release release;
        result = dynsyslab::core::integrate(*system, initialStates, batch, options);
    }
    py::array_t<double> times(result.time.size());
    std::copy(result.time.begin(), result.time.end(), static_cast<double*>(times.request().ptr));
    py::array_t<double> values = initial.ndim == 1
        ? py::array_t<double>({static_cast<py::ssize_t>(result.sampleCount),
                              static_cast<py::ssize_t>(result.stateDimension)})
        : py::array_t<double>({static_cast<py::ssize_t>(result.trajectoryCount),
                              static_cast<py::ssize_t>(result.sampleCount),
                              static_cast<py::ssize_t>(result.stateDimension)});
    std::copy(result.values.begin(), result.values.end(), static_cast<double*>(values.request().ptr));
    py::dict metadata(
        "coordinate_system"_a = result.coordinateSystem,
        "accepted_steps"_a = result.statistics.acceptedSteps,
        "rejected_steps"_a = result.statistics.rejectedSteps,
        "rhs_calls"_a = result.statistics.rightHandSideCalls);
    return py::make_tuple(times, values, metadata);
}

py::tuple kromRollout(
    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> eigenvalueArray,
    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> modeArray,
    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> initialArray,
    py::array_t<double, py::array::c_style | py::array::forcecast> timeArray) {
    const auto eigenvalues = eigenvalueArray.request();
    const auto modes = modeArray.request();
    const auto initial = initialArray.request();
    const auto times = timeArray.request();
    if (eigenvalues.ndim != 1 || modes.ndim != 2 || initial.ndim != 1 || times.ndim != 1)
        throw std::invalid_argument("KROM arrays have invalid dimensions");
    dynsyslab::core::KROMRolloutRequest request;
    const auto* lambda = static_cast<const std::complex<double>*>(eigenvalues.ptr);
    const auto* mode = static_cast<const std::complex<double>*>(modes.ptr);
    const auto* coordinate = static_cast<const std::complex<double>*>(initial.ptr);
    const auto* time = static_cast<const double*>(times.ptr);
    request.eigenvalues.assign(lambda, lambda + eigenvalues.shape[0]);
    request.modes.assign(mode, mode + modes.shape[0] * modes.shape[1]);
    request.initialCoordinates.assign(coordinate, coordinate + initial.shape[0]);
    request.time.assign(time, time + times.shape[0]);
    request.stateDimension = static_cast<std::size_t>(modes.shape[0]);
    dynsyslab::core::KROMRollout result;
    {
        py::gil_scoped_release release;
        result = dynsyslab::core::rolloutKROM(request);
    }
    py::array_t<std::complex<double>> lifted({static_cast<py::ssize_t>(result.samples),
                                              static_cast<py::ssize_t>(result.modes)});
    py::array_t<std::complex<double>> reconstructed({static_cast<py::ssize_t>(result.samples),
                                                     static_cast<py::ssize_t>(result.stateDimension)});
    std::copy(result.lifted.begin(), result.lifted.end(),
              static_cast<std::complex<double>*>(lifted.request().ptr));
    std::copy(result.reconstructed.begin(), result.reconstructed.end(),
              static_cast<std::complex<double>*>(reconstructed.request().ptr));
    return py::make_tuple(lifted, reconstructed);
}

py::tuple ising2DRun(
    py::array_t<std::int64_t, py::array::c_style | py::array::forcecast> spinArray,
    double coupling, double field, double temperature, std::uint64_t seed, py::ssize_t sweeps) {
    const auto input = spinArray.request();
    if (input.ndim != 2 || input.shape[0] != input.shape[1] || sweeps < 0)
        throw std::invalid_argument("native Ising array or sweep count is invalid");
    dynsyslab::core::IsingRequest request;
    request.size = static_cast<std::size_t>(input.shape[0]);
    const auto* spins = static_cast<const std::int64_t*>(input.ptr);
    request.spins.assign(spins, spins + request.size * request.size);
    request.coupling = coupling;
    request.field = field;
    request.temperature = temperature;
    request.seed = seed;
    request.sweeps = static_cast<std::size_t>(sweeps);
    dynsyslab::core::IsingResult result;
    {
        py::gil_scoped_release release;
        result = dynsyslab::core::runIsing2D(request);
    }
    py::array_t<double> observations({sweeps, py::ssize_t{2}});
    py::array_t<std::int64_t> finalSpins({input.shape[0], input.shape[1]});
    std::copy(result.observations.begin(), result.observations.end(),
              static_cast<double*>(observations.request().ptr));
    std::copy(result.spins.begin(), result.spins.end(),
              static_cast<std::int64_t*>(finalSpins.request().ptr));
    return py::make_tuple(observations, finalSpins);
}

}  // namespace

PYBIND11_MODULE(_dynsyslab_native, module) {
    module.doc() = "Optional Python bindings for the standalone DynSysLab native core";
    module.def("capabilities", []() {
        return py::dict(
            "version"_a = std::string(dynsyslab::core::version),
            "fixed_step_ode"_a = true,
            "adaptive_ode"_a = true,
            "batched_ode"_a = true,
            "fixed_step_models"_a = dynsyslab::core::fixedStepModels(),
            "system_registry"_a = true,
            "krom_rollout"_a = true,
            "krom_rollout_auto"_a = false,
            "pde"_a = false,
            "lattice"_a = true,
            "ising2d"_a = true,
            "ising2d_auto"_a = false,
            "python_required"_a = false);
    });
    module.def("integrate_fixed", &integrateFixed, py::arg("model"), py::arg("parameters"),
               py::arg("initial_states"), py::arg("t0"), py::arg("tf"), py::arg("step_size"),
               py::arg("method"));
    module.def("integrate_system", &integrateSystem, py::arg("identifier"), py::arg("preset"),
               py::arg("parameters"), py::arg("initial_states"),
               py::arg("t0"), py::arg("tf"), py::arg("output_step"), py::arg("method") = "RK45",
               py::arg("rtol") = 1.0e-6, py::arg("atol") = 1.0e-9);
    module.def("systems", []() {
        static_cast<void>(dynsyslab::core::mackeyGlassDescriptor());
        static_cast<void>(dynsyslab::core::ising2DDescriptor());
        static_cast<void>(dynsyslab::core::heatEquationDescriptor());
        py::list result;
        for (const auto& descriptor : dynsyslab::core::availableSystems()) {
            py::dict item("id"_a = descriptor.identifier, "label"_a = descriptor.label,
                          "state_dimension"_a = descriptor.stateDimension,
                          "family"_a = descriptor.family,
                          "coordinate_system"_a = descriptor.preferredCoordinates,
                          "state_names"_a = descriptor.stateNames);
            result.append(std::move(item));
        }
        return result;
    });
    module.def("krom_rollout", &kromRollout, py::arg("eigenvalues"), py::arg("modes"),
               py::arg("initial_coordinates"), py::arg("time"));
    module.def("ising2d_run", &ising2DRun, py::arg("spins"), py::arg("coupling"), py::arg("field"),
               py::arg("temperature"), py::arg("seed"), py::arg("sweeps"));
}
