#include <dynsyslab/extension_api.h>

#include <cmath>
#include <new>

namespace {

struct Oscillator {
    double frequency{2.0};
};

void* create(const char*, char*, size_t) {
    return new (std::nothrow) Oscillator;
}

void destroy(void* instance) {
    delete static_cast<Oscillator*>(instance);
}

const char* resolvedConfiguration(const void*) {
    return R"({"equations_version":"example.oscillator/1","parameters":{"frequency":2.0},"parameter_units":{"frequency":"rad/time"},"initial_state":[1.0,0.0]})";
}

int rightHandSide(const void* raw, double, const double* state, size_t stateCount,
                  double input, double* derivative, size_t derivativeCount,
                  char*, size_t) {
    if (!raw || !state || !derivative || stateCount != 2 || derivativeCount != 2)
        return 1;
    const double frequency = static_cast<const Oscillator*>(raw)->frequency;
    derivative[0] = state[1];
    derivative[1] = -frequency * frequency * state[0] + input;
    return 0;
}

int inputMatrix(const void*, const double*, size_t stateCount, double* matrix,
                size_t matrixCount, char*, size_t) {
    if (!matrix || stateCount != 2 || matrixCount != 2) return 1;
    matrix[0] = 0.0;
    matrix[1] = 1.0;
    return 0;
}

int invalid(const void*, double, const double* state, size_t stateCount) {
    return !state || stateCount != 2 || !std::isfinite(state[0])
        || !std::isfinite(state[1]);
}

size_t observableCount(const void*) {
    return 1;
}

const char* observableDescriptor(const void*, size_t index) {
    return index == 0
        ? R"({"identifier":"energy","label":"Energy","units":"dimensionless","coordinate_system":"dimensional","output_dimension":1})"
        : nullptr;
}

int evaluateObservable(const void* raw, size_t index, const double* state,
                       size_t stateCount, double* values, size_t valueCount,
                       char*, size_t) {
    if (!raw || index != 0 || !state || stateCount != 2 || !values || valueCount != 1)
        return 1;
    const double frequency = static_cast<const Oscillator*>(raw)->frequency;
    values[0] = 0.5 * (state[1] * state[1]
                       + frequency * frequency * state[0] * state[0]);
    return 0;
}

const dsl_system_extension_v1 systems[]{
    {sizeof(dsl_system_extension_v1),
     R"({"identifier":"example_oscillator","label":"Example oscillator","state_dimension":2,"preferred_coordinates":"dimensional","state_names":["position","velocity"],"state_labels":["x","dx/dt"],"presets":[]})",
     create, destroy, resolvedConfiguration, nullptr, rightHandSide, nullptr,
     inputMatrix, invalid, observableCount, observableDescriptor,
     evaluateObservable, nullptr}
};

const dsl_extension_v1 extension{
    sizeof(dsl_extension_v1), DSL_EXTENSION_ABI_VERSION,
    "dynsyslab.example.oscillator", "1.0.0",
    R"({"systems":["example_oscillator"],"observables":["energy"]})",
    1, systems, nullptr};

}  // namespace

extern "C" DSL_EXTENSION_EXPORT const dsl_extension_v1* dynsyslab_extension_v1(void) {
    return &extension;
}
