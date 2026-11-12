from __future__ import annotations

import tempfile
import unittest

import numpy as np

from archives.store import serializable_system
from dynsyslib import load_system, resolve_archived_system
from dynsyslib.drivers import SinusoidalDriver
from solver import SolveRequest, SolverOptions, solve
from postprocessing.core import Trajectory, integration_coordinates
from time_series import TimeSeries
from viewer.sources import axis_context


class RPENondimensionalizationTests(unittest.TestCase):
    def test_default_physical_parameters_and_natural_period(self):
        system=load_system("rpe",preset="air_water_20c")
        physical=system.dimensional_parameters; groups=system.dimensionless_parameters
        self.assertEqual(physical.equilibrium_radius_m,1e-5)
        self.assertEqual(physical.kappa,4/3)
        self.assertEqual(physical.ambient_pressure_pa,1e5)
        self.assertEqual(physical.vapor_pressure_pa,2330)
        self.assertEqual(physical.rho_kg_m3,998.)
        self.assertEqual(physical.sigma_n_m,.0725)
        self.assertEqual(physical.mu_pa_s,.001)
        expected=np.sqrt(3*physical.kappa*physical.initial_gas_pressure_pa-2*physical.sigma_n_m/physical.equilibrium_radius_m)/(2*np.pi*physical.equilibrium_radius_m*np.sqrt(physical.rho_kg_m3))
        self.assertAlmostEqual(groups.natural_frequency_hz,expected)
        self.assertAlmostEqual(groups.time_scale_s,1/expected)

    def test_dimensionless_rhs_is_equivalent_to_dimensional_rhs(self):
        system=load_system("rpe",preset="air_water_20c",driver=SinusoidalDriver(210000.,1700.,.3))
        state=np.array((1.37e-5,-2.4)); transformed=system.to_dimensionless_state(state)
        dimensional=system.rhs(2e-6,state); dimensionless=system.dimensionless_rhs(float(system.to_dimensionless_time(2e-6)),transformed)
        scale=system.dimensionless_parameters.time_scale_s; radius=system.equilibrium_radius_m
        expected=np.array((scale/radius*dimensional[0],scale*scale/radius*dimensional[1]))
        np.testing.assert_allclose(dimensionless,expected,rtol=2e-13,atol=2e-13)

    def test_dimensionless_input_matrix_uses_control_pressure_ratio(self):
        system=load_system("rpe",preset="air_water_20c"); state=np.asarray((1.37,0.2)); drive=.01234
        delta=system.dimensionless_rhs_with_input(.3,state,drive)-system.dimensionless_rhs_with_input(.3,state,0.)
        np.testing.assert_allclose(delta,system.dimensionless_input_matrix(state)@np.asarray((drive,)),rtol=2e-13,atol=2e-13)
        physical=system.to_dimensional_state(state)
        # The second dimensionless state derivative is with respect to tau,
        # hence one additional T; u*=p_control/p0 supplies the p0 factor.
        expected=np.asarray((0.,(system.dimensionless_parameters.time_scale_s**2/system.equilibrium_radius_m)*system.input_matrix(physical)[1,0]*system.ambient_pressure_pa)).reshape(2,1)
        np.testing.assert_allclose(system.dimensionless_input_matrix(state),expected,rtol=2e-13,atol=2e-13)
        self.assertAlmostEqual(system.dimensionless_input_matrix(state)[1,0],-system.dimensionless_parameters.ambient_pressure/state[0])
        self.assertAlmostEqual(float(system.to_dimensionless_input(system.to_dimensional_input(drive))),drive)

    def test_positive_external_pressure_compresses_and_legacy_mode_is_reproducible(self):
        physical=load_system("rpe",preset="air_water_20c")
        equilibrium=physical.initial_state
        pressure=2500.
        base=physical.rhs_with_input(0.,equilibrium,0.)
        driven=physical.rhs_with_input(0.,equilibrium,pressure)
        self.assertEqual(physical.pressure_input_convention,"external-acoustic-pressure/v2")
        self.assertLess(physical.input_matrix(equilibrium)[1,0],0.)
        self.assertLess(driven[1],base[1])
        dimensionless=physical.dimensionless_initial_state
        dimensionless_pressure=physical.to_dimensionless_input(pressure)
        self.assertLess(physical.dimensionless_rhs_with_input(0.,dimensionless,dimensionless_pressure)[1],
                        physical.dimensionless_rhs_with_input(0.,dimensionless,0.)[1])

        legacy=physical.with_dimensional_parameters(pressure_input_convention="legacy-positive-acceleration/v1")
        self.assertGreater(legacy.input_matrix(equilibrium)[1,0],0.)
        self.assertGreater(legacy.rhs_with_input(0.,equilibrium,pressure)[1],
                           legacy.rhs_with_input(0.,equilibrium,0.)[1])

    def test_archives_without_pressure_convention_reconstruct_legacy_semantics(self):
        parameters=serializable_system(load_system("rpe",preset="air_water_20c"))
        parameters.pop("pressure_input_convention")
        restored=resolve_archived_system("rayleigh-plesset",parameters)
        self.assertEqual(restored.pressure_input_convention,"legacy-positive-acceleration/v1")
        explicit=resolve_archived_system("rayleigh-plesset",serializable_system(load_system("rpe",preset="air_water_20c")))
        self.assertEqual(explicit.pressure_input_convention,"external-acoustic-pressure/v2")

    def test_dimensionless_jacobian_matches_finite_difference(self):
        system=load_system("rpe",preset="air_water_20c",driver=SinusoidalDriver(220000.,1800.,.4))
        tau=0.37; state=np.array((1.08,-0.015)); step=1e-6
        numerical=np.column_stack([
            (system.dimensionless_rhs(tau,state+step*np.eye(2)[column])-
             system.dimensionless_rhs(tau,state-step*np.eye(2)[column]))/(2*step)
            for column in range(2)
        ])
        np.testing.assert_allclose(system.dimensionless_jacobian(tau,state),numerical,rtol=2e-7,atol=2e-8)

    def test_solver_archives_scaling_and_viewer_prefers_it(self):
        system=load_system("rpe",preset="air_water_20c"); scale=system.dimensionless_parameters.time_scale_s
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=system,kind="ode",time_span=(0.,2*scale),archive_root=directory,
                options=SolverOptions(method="DOP853",use_dimensionless_form=True)))
        self.assertEqual(result.status,"complete")
        self.assertEqual(result.metadata["nondimensionalization"]["convention"],"rpe-natural-period/v1")
        context=axis_context({"system":"rayleigh-plesset","integration_form":"dimensionless",
            "system_parameters":{"equilibrium_radius_m":system.equilibrium_radius_m,"ambient_pressure_pa":system.ambient_pressure_pa,
                                 "preset":{"rho_kg_m3":system.preset.rho_kg_m3}},
            "nondimensionalization":result.metadata["nondimensionalization"]})
        self.assertAlmostEqual(context["time_scale"],1/scale)

    def test_postprocessing_recovers_dimensionless_integration_coordinates_without_mutation(self):
        system=load_system("rpe",preset="air_water_20c"); scale=system.dimensionless_parameters.time_scale_s
        physical=TimeSeries(np.array((0.,scale,2*scale)),np.array(((2*system.equilibrium_radius_m,0.),
            (1.5*system.equilibrium_radius_m,system.equilibrium_radius_m/scale),
            (system.equilibrium_radius_m,0.))))
        metadata={"system":system.name,"system_parameters":serializable_system(system),
                  "integration_form":"dimensionless","uniform_time_steps":True,"delta_t":scale,
                  "nondimensionalization":{"time_scale_s":scale}}
        source=Trajectory(physical,metadata,None); transformed=integration_coordinates(source)
        np.testing.assert_allclose(transformed.time,np.array((0.,1.,2.)))
        np.testing.assert_allclose(transformed.values,np.array(((2.,0.),(1.5,1.),(1.,0.))))
        np.testing.assert_array_equal(source.values,physical.values)
        self.assertEqual(transformed.metadata["postprocessing_preprocessing"]["coordinate_system"],"dimensionless")


if __name__=="__main__":
    unittest.main()
