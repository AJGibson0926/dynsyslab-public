from __future__ import annotations

import unittest

import numpy as np

from dynsyslib import Observable
from postprocessing.koopman import (
    KoopmanPipelineSpecification,
    ScatteredEigenfunctionInterpolator,
    fit_boundary_eigenfunctions,
    optimize_trajectory_eigenvalues,
    run_composable_koopman_pipeline,
)
from time_series import TimeSeries


class ComposableKoopmanPipelineTests(unittest.TestCase):
    class _Line:
        name="line"; state_dimension=1; initial_state=np.zeros(1)
        def input_matrix(self,state): return np.ones((1,1))

    def _linear_ensemble(self):
        time=.1*np.arange(25); initial=np.asarray((.4,.8,1.2))
        values=np.empty((time.size,initial.size,1))
        values[:,:,0]=np.exp(-time[:,None])*initial[None,:]
        return TimeSeries(time,values,pooled_members=initial.size)

    def test_scattered_representation_uses_automatic_finite_difference(self):
        states=np.column_stack((np.linspace(-1.,1.,41),np.linspace(-.5,.5,41)))
        values=states[:,0]+2.*states[:,1]
        interpolator=ScatteredEigenfunctionInterpolator(states,values,neighbors=4)
        from postprocessing.koopman.spectrum import DynamicalSystemReference, KoopmanEigenfunction, KoopmanEigenvalue
        reference=DynamicalSystemReference(None,"test",{})
        function=KoopmanEigenfunction(KoopmanEigenvalue.from_continuous(-1.,.1,reference),
            "scattered",values=values,state_callable=interpolator,representation_kind="scattered-interpolation")
        self.assertEqual(function.gradient_backend,"finite_difference")
        self.assertEqual(function.summary()["representation_kind"],"scattered-interpolation")
        self.assertTrue(np.all(np.isfinite(function.nabla_psi(np.asarray((.1,.05))))))

    def test_optimization_is_eigenvalue_only_and_boundary_stage_reuses_values(self):
        series=self._linear_ensemble(); system=self._Line()
        observable=Observable.component(1,0,name="x")
        result=optimize_trajectory_eigenvalues(series,observable=observable,
            initial_eigenvalue_groups=((np.exp(-.09),),),system=system,
            metadata={"system":"line","system_parameters":{}},max_iterations=60)
        self.assertEqual(result.spectrum.eigenfunctions,[])
        self.assertAlmostEqual(result.spectrum.eigenvalues[0].discrete_time.real,np.exp(-.1),places=5)
        boundary=fit_boundary_eigenfunctions(result.spectrum,series,observable=observable)
        self.assertEqual(boundary.spectrum.eigenvalues[0],result.spectrum.eigenvalues[0])
        self.assertEqual(boundary.spectrum.eigenfunctions[0].representation_kind,"scattered-interpolation")
        self.assertEqual(boundary.spectrum.eigenfunctions[0].gradient_backend,"finite_difference")

    def test_pipeline_spec_selects_optimization_then_boundary(self):
        series=self._linear_ensemble(); system=self._Line()
        observable=Observable.component(1,0,name="x")
        result=run_composable_koopman_pipeline(
            KoopmanPipelineSpecification("optimization","boundary-functions"),
            eigenvalue_options={"series":series,"observable":observable,
                "initial_eigenvalue_groups":((np.exp(-.09),),),"system":system,
                "metadata":{"system":"line","system_parameters":{}},"max_iterations":60},
            eigenfunction_options={"series":series,"observable":observable})
        self.assertEqual(result.eigenvalues.method,"optimization")
        self.assertEqual(result.eigenfunctions.method,"boundary-functions")


if __name__=="__main__": unittest.main()
