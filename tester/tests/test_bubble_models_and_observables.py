from __future__ import annotations

import unittest

import numpy as np

from archives.store import serializable_system
from dynsyslib import (MARMOTTANT_PRESETS, REDDY_SZERI_PRESETS, Observable,
                       load_system, observable_from_summary, resolve_archived_system,
                       state_variables)


class BubbleModelAndObservableTests(unittest.TestCase):
    def test_every_ode_has_state_components_and_full_state(self):
        system=load_system("van_der_pol")
        catalog=system.observable_catalog("native")
        self.assertEqual(tuple(catalog)[:3],("x0","x1","state"))
        np.testing.assert_allclose(catalog["state"].evaluate((2.,3.)),(2.,3.))
        np.testing.assert_allclose(catalog["x[1]"].evaluate((2.,3.)),(3.,))

    def test_rpe_catalog_contains_matlab_balance_observables(self):
        system=load_system("rpe")
        catalog=system.observable_catalog("dimensionless")
        self.assertTrue({"radius","radial_velocity","state","momentum_balance",
                         "far_field_acoustic_echo"}.issubset(catalog))
        self.assertAlmostEqual(float(catalog["momentum_balance"].evaluate((1.,0.))[0]),0.,places=11)

    def test_rpe_echo_is_distance_scaled_monopole_in_either_chart(self):
        system=load_system("rpe")
        for state in ((.55,1.2),(1.4,-.6)):
            native=system.to_dimensional_state(state)
            R,V=native
            acceleration=system.rhs_with_input(0.,native,0.)[1]
            pressure_distance=system.preset.rho_kg_m3*(2*R*V**2+R**2*acceleration)
            for chart,point in (("dimensionless",state),("dimensional",native)):
                value=system.observable_catalog(chart)["far_field_acoustic_echo"].evaluate(point)[0]
                np.testing.assert_allclose(value,pressure_distance,rtol=1e-12)

    def test_duffing_catalog_preserves_complex_matlab_observation(self):
        system=load_system("duffing",preset="gibson_damped_stabilization_2022")
        catalog=system.observable_catalog("dimensionless")
        value=catalog["random_observation"].evaluate((.25,-.5))[0]
        self.assertTrue(np.iscomplexobj(value))
        self.assertAlmostEqual(catalog["potential_energy"].evaluate((.5,0.))[0],-.125+.25*4*.5**4)

    def test_marmottant_is_two_state_with_resettable_discrete_hysteresis(self):
        system=load_system("marmottant",preset="matlab_cleanup_linear")
        self.assertEqual(system.state_dimension,2)
        np.testing.assert_allclose(system.dimensionless_rhs(0.,(1.,0.)),(0.,0.),atol=1e-12)
        groups=system.dimensionless_parameters
        system.surface_tension(groups.break_radius+.01)
        self.assertTrue(system.hysteresis_ruptured)
        ruptured=system.surface_tension((groups.rupture_radius+groups.break_radius)/2,update_hysteresis=False)
        self.assertAlmostEqual(ruptured,system.preset.sigma_n_m)
        system.reset_hysteresis()
        intact=system.surface_tension((groups.rupture_radius+groups.break_radius)/2,update_hysteresis=False)
        self.assertNotAlmostEqual(intact,ruptured)
        system.surface_tension(groups.rupture_radius-.01)
        self.assertFalse(system.hysteresis_ruptured)

    def test_reddy_szeri_retains_twelve_states_and_full_control_map(self):
        system=load_system("reddy_szeri",preset="matlab_cleanup_reference")
        self.assertEqual(system.state_dimension,12)
        matrix=system.dimensionless_input_matrix(system.initial_state)
        self.assertEqual(matrix.shape,(12,1))
        self.assertNotEqual(matrix[5,0],0.)
        catalog=system.observable_catalog("dimensionless")
        self.assertTrue({"shape_5","compressibility_loss","radial_acceleration","k2","k5","g1","state"}.issubset(catalog))

    def test_analytic_observable_stashes_exact_partial_observables(self):
        x0,x1=state_variables(2)
        observable=Observable.from_expressions("custom",(x0*x1+x0**2,),
            validation_state=(1.,2.),derivative_axis_names=("position","velocity"))
        probe=np.asarray((3.,4.))
        np.testing.assert_allclose(observable.evaluate(probe),(21.,))
        np.testing.assert_allclose(observable.jacobian(probe),((10.,3.),))
        self.assertEqual(observable.derivative_observables[0][0].name,"custom_position_derivative")
        restored=observable_from_summary(observable.summary())
        np.testing.assert_allclose(restored.evaluate(probe),observable.evaluate(probe))
        np.testing.assert_allclose(restored.jacobian(probe),observable.jacobian(probe))

    def test_reddy_observable_coordinates_have_exact_jacobians_and_presets(self):
        system=load_system("reddy_szeri",preset="matlab_cleanup_reference")
        state=system.initial_state.copy(); state[1]=.07
        for name in ("c","k2","k5"):
            observable=system.observable(name,coordinate_system="dimensionless")
            analytic=observable.jacobian(state)
            numerical=observable.jacobian(state,backend="finite_difference",relative_step=1e-5)
            np.testing.assert_allclose(analytic,numerical,rtol=2e-4,atol=2e-5)
            self.assertIsNotNone(observable.derivative_observables)
        radial=system.observable_coordinates("radial_state",coordinate_system="dimensionless")
        np.testing.assert_allclose(radial.jacobian(state)[:,:2],np.eye(2))
        np.testing.assert_allclose(radial.jacobian(state)[:,2:],0.)
        transformed=system.observable_coordinates("c_k5",coordinate_system="dimensionless")
        self.assertEqual(transformed.output_dimension,2)
        self.assertEqual(transformed.jacobian(state).shape,(2,12))
        restored=observable_from_summary(transformed.summary())
        np.testing.assert_allclose(restored.evaluate(state),transformed.evaluate(state))
        np.testing.assert_allclose(restored.jacobian(state),transformed.jacobian(state))
        legacy=system._legacy_observable_terms(state)
        self.assertAlmostEqual(system.observable("c_time_derivative").evaluate(state)[0],legacy[3])
        for name in ("k2","k5"):
            observable=system.observable(name)
            expected=observable.jacobian(state)[0,0]*state[1]+observable.jacobian(state)[0,1]*legacy[1]
            self.assertAlmostEqual(system.observable(f"{name}_time_derivative").evaluate(state)[0],expected)

    def test_matlab_training_provenance_is_encoded_in_presets(self):
        marm=MARMOTTANT_PRESETS["matlab_cleanup_linear"].processing_metadata["training"]
        reddy=REDDY_SZERI_PRESETS["matlab_cleanup_reference"].processing_metadata["training"]
        self.assertEqual(marm["seed"],2141442)
        self.assertEqual(marm["initial_condition_ranges"]["radial_velocity"],[-10.,10.])
        self.assertEqual(reddy["delta_tau"],.001)

    def test_new_bubble_models_round_trip_from_archive_metadata(self):
        for key,preset in (("marmottant","matlab_cleanup_linear"),("reddy_szeri","matlab_cleanup_reference")):
            original=load_system(key,preset=preset)
            restored=resolve_archived_system(original.name,serializable_system(original))
            self.assertEqual(restored.state_dimension,original.state_dimension)
            self.assertEqual(restored.preset.key,original.preset.key)


if __name__=="__main__":
    unittest.main()
