// rayleigh-plesset-dimensional-baseline: the frozen pre-nondimensionalisation RPE of the
// Python dynsyslib/lib/rpe_dimensional_baseline.py, without a second copy of the physics.
//
// The baseline integrates, in SI units,
//     R'' = [ p_g0 (R0/R)^(3 kappa) + p_v - (p_inf - p_d(t)) - 2 sigma/R - 4 mu R'/R ] / (rho R)
//           - 3 R'^2 / (2 R),   p_g0 = p_inf - p_v + 2 sigma / R0,
// with p_inf = 101325 Pa, R0 = 10 um, the air-water 20 C medium and the driver entering as
// +p_d(t) (the external pressure is p_inf - p_d). That is the native `rpe` model with
// ambient_pressure_pa = 101325 and pressure_input_sign = +1, integrated through the exact
// dimensional change of variables (DimensionalSystem.hpp). Rayleigh, Phil. Mag. 34 (1917)
// 94-98; Plesset, J. Appl. Mech. 16 (1949) 277-282; scales after Gibson, Calvisi & Yee,
// SIAM J. Appl. Dyn. Syst. 21 (2022), doi:10.1137/21M1456078.
#include "BubbleSupport.hpp"

#include <dynsyslab/core/DimensionalSystem.hpp>

namespace dynsyslab::core::detail {

std::unique_ptr<ODESystem> makeRPEDimensionalBaseline(const SystemConfiguration& configuration) {
    auto model = configuration;
    model.identifier = "rpe";
    if (model.preset.empty()) model.preset = "air_water_20c";
    model.parameters.emplace("ambient_pressure_pa", 101325.0);
    model.parameters.emplace("pressure_input_sign", 1.0);
    // An initial state is given in metres and m/s; the model takes R/R0 and dR/dtau.
    const auto physicalInitial = model.initialState;
    model.initialState.clear();
    auto rpe = makeRPE(model);
    if (!physicalInitial.empty()) {
        const auto scales = rpe->coordinates().physicalStatePerModelState;
        if (physicalInitial.size() != scales.size())
            throw std::invalid_argument("rayleigh-plesset-dimensional-baseline initial state needs radius and velocity");
        for (std::size_t index = 0; index < scales.size(); ++index)
            model.initialState.push_back(physicalInitial[index] / scales[index]);
        rpe = makeRPE(model);
    }
    return dimensionalSystem(std::move(rpe), descriptorFor("rayleigh-plesset-dimensional-baseline"));
}

}  // namespace dynsyslab::core::detail
