#include "PhaseSurfaceInternal.hpp"
#include <algorithm>
namespace dynsyslab::core::phase_detail {
void traceSurface(const ODESystem &system, PhaseSpaceResult &result, std::size_t index, bool stable,
                  const ExecutionContext *execution) {
    const auto &o = result.options;
    PhaseSurface mesh;
    mesh.status = "unresolved";
    try {
        const auto basis = surface::tangentBasis(result.fixedPoints[index], stable, o);
        mesh = surface::grow(system, result.fixedPoints[index].state, basis, stable ? -1. : 1., o,
                             execution);
        auto fine = o;
        fine.surface.seedRadius *= .5;
        fine.surface.step *= .5;
        fine.surface.meshSize *= .5;
        fine.surface.chordTolerance *= .25;
        fine.surface.initialPoints *= 2;
        fine.surface.maximumFronts = std::min(std::size_t(512), 2 * o.surface.maximumFronts + 16);
        fine.surface.maximumPoints = std::min(std::size_t(65536), 4 * o.surface.maximumPoints);
        fine.maximumStep *= .5;
        fine.relativeTolerance *= .125;
        fine.absoluteTolerance *= .125;
        const auto refined = surface::grow(system, result.fixedPoints[index].state, basis,
                                           stable ? -1. : 1., fine, execution);
        surface::validate(system, mesh, refined, o, execution);
    } catch (const ExecutionCancelled &) {
        throw;
    } catch (const std::exception &error) {
        mesh.status = "unresolved";
        mesh.reason = error.what();
    }
    mesh.kind = stable ? "stable_manifold_surface" : "unstable_manifold_surface";
    mesh.baseObject = index;
    result.manifoldSurfaces.push_back(std::move(mesh));
}
} // namespace dynsyslab::core::phase_detail
