#include <dynsyslab/extension_api.h>

#include <cmath>
#include <cstring>
#include <new>

namespace {
struct Oscillator { double omega{2.0}; };

void* create(const char*, char*, size_t) { return new (std::nothrow) Oscillator; }
void destroy(void* value) { delete static_cast<Oscillator*>(value); }
const char* resolved(const void*) {
    return R"({"equations_version":"test.oscillator/1","parameters":{"omega":2},"parameter_units":{"omega":"rad/time"},"initial_state":[1,0]})";
}
int rhs(const void* raw,double,const double* x,size_t n,double input,
        double* out,size_t m,char*,size_t) {
    if (!raw||n!=2||m!=2) return 1;
    const double omega=static_cast<const Oscillator*>(raw)->omega;
    out[0]=x[1]; out[1]=-omega*omega*x[0]+input; return 0;
}
int inputMatrix(const void*,const double*,size_t n,double* out,size_t m,char*,size_t) {
    if (n!=2||m!=2) return 1; out[0]=0.0; out[1]=1.0; return 0;
}
int invalid(const void*,double,const double* x,size_t n) {
    return n!=2||!std::isfinite(x[0])||!std::isfinite(x[1]);
}
size_t observableCount(const void*) { return 1; }
const char* observableDescriptor(const void*,size_t index) {
    return index==0 ? R"({"identifier":"energy","label":"Energy","units":"dimensionless","coordinate_system":"dimensional","output_dimension":1})" : nullptr;
}
int observable(const void* raw,size_t index,const double* x,size_t n,
               double* out,size_t m,char*,size_t) {
    if (!raw||index||n!=2||m!=1) return 1;
    const double omega=static_cast<const Oscillator*>(raw)->omega;
    out[0]=0.5*(x[1]*x[1]+omega*omega*x[0]*x[0]); return 0;
}
int execute(const char*,const char*,dsl_extension_output_writer writer,
            void* user,char*,size_t) {
    const char output[]=R"({"extension_result":true})";
    writer(output,sizeof(output)-1,user); return 0;
}

const dsl_system_extension_v1 systems[]{
    {sizeof(dsl_system_extension_v1),
     R"({"identifier":"test_oscillator","label":"Test oscillator","state_dimension":2,"preferred_coordinates":"dimensional","state_names":["position","velocity"],"state_labels":["x","dx/dt"],"presets":[]})",
     create,destroy,resolved,nullptr,rhs,nullptr,inputMatrix,invalid,
     observableCount,observableDescriptor,observable,nullptr}
};
const dsl_extension_v1 extension{
    sizeof(dsl_extension_v1),DSL_EXTENSION_ABI_VERSION,"dynsyslab.test","1.0.0",
    R"({"systems":["test_oscillator"],"observables":["energy"],"postprocessors":["test"],"controllers":["test"]})",
    1,systems,execute};
}

extern "C" DSL_EXTENSION_EXPORT const dsl_extension_v1* dynsyslab_extension_v1(void) {
    return &extension;
}
