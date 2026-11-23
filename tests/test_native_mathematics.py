"""Public contract tests: the Python API delegates numerical work to C++."""
import math
from fractions import Fraction
import os
from pathlib import Path
import unittest
from dynsyslab_client.mathematics import query

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("DYNSYSLAB_NATIVE_BIN", ROOT / "build/native-suite-fast/bin"))
CLI = BIN / ("dynsyslab-cli.exe" if os.name == "nt" else "dynsyslab-cli")


@unittest.skipUnless(CLI.is_file(), "Build the native CLI first")
class MathematicsContractTest(unittest.TestCase):
    def call(self, **request):
        return query(request, executable=CLI)["result"]

    def test_cover_forcing_and_resonant_branch_resolution(self):
        exact=self.call(operation="phase_cover_resonance",prime=2,base_to_drive_ratio="4/3")
        self.assertEqual((exact['character_ratio'],exact['drive_harmonic'],exact['minimum_depth']),
                         ('3/4','-1',2))
        data=self.call(operation="phase_cover_forcing",prime=2,depth=3,
            base_angular_frequency=1,forcing_frequencies=[1],divisor_tolerance=1e-10,
            coefficients=[{'index':[8,-1],'coefficient':3},
                          {'index':[1,0],'coefficient':[0,2]}])
        self.assertFalse(data['kam_certified'])
        self.assertEqual(data['retained_mean_square'],9)
        self.assertEqual(data['solution']['coefficients'][0]['coefficient'],[16,0])
        resonant=next(m for m in data['modes'] if m['index']==[8,-1])
        self.assertIsNone(resonant['correction_magnitude'])
        self.assertTrue(resonant['retained'])

    def test_exact_arithmetic_and_topology(self):
        self.assertEqual(self.call(operation="continued_fraction", value="355/113")["convergents"],
                         ["3", "22/7", "355/113"])
        self.assertEqual(self.call(operation="padic", a="4/9", b="1/9",
                                  prime=3, precision=4)["distance"], "3")
        result = self.call(operation="graph", vertices=4, edges=[[0, 1], [1, 2]],
                           source=0, vertex_values=[0, 2, 1, 3])
        self.assertEqual(result["components"], 2)
        self.assertEqual(result["distances"], [0, 1, 2, None])

    def test_torus_clock_eigenfunctions_and_time_inverse(self):
        clock=dict(angular_frequencies=[2],coefficients=[
            dict(index=[0],coefficient=1.3),
            dict(index=[1],coefficient=.2),dict(index=[-1],coefficient=.2)])
        result=self.call(operation="torus_clock",clock_rate=clock,initial_phases=[.7],
                         time=[0,3,-2],characters=[[1]])
        self.assertFalse(result["infinite_limit_certified"])
        self.assertFalse(result["interval_certified"])
        self.assertAlmostEqual(result["clock_lower_bound"],.9)
        initial=result["samples"][0]["straightened_phases"][0]
        for row in result["samples"]:
            self.assertAlmostEqual(row["straightened_phases"][0],initial+2*row["time"]/1.3,places=10)
            self.assertLess(abs(row["time_residual"]),1e-10)
        with self.assertRaises(RuntimeError):
            self.call(operation="torus_clock",clock_rate=dict(
                angular_frequencies=[1,1],coefficients=[
                    dict(index=[0,0],coefficient=1),
                    dict(index=[1,-1],coefficient=.1),
                    dict(index=[-1,1],coefficient=.1)]))

    def test_fourier_and_distributions(self):
        s = {"angular_frequency": 2, "coefficients": [[-1, .5, 0], [1, .5, 0]]}
        ds = self.call(operation="fourier", action="differentiate", series=s)
        value = self.call(operation="fourier", action="evaluate", series=ds, phase=.4)["value"]
        self.assertAlmostEqual(value[0], -2*math.sin(.4), places=13)
        d = {"location": 0, "coefficients": [0, 0, 1]}
        product = self.call(operation="distribution", action="multiply_smooth",
                            distribution=d, smooth_derivatives=[0, 1, 0])
        self.assertEqual(product["coefficients"], [[0, 0], [-2, 0], [0, 0]])

    def test_exact_phase_cover_prefix(self):
        levels=self.call(operation="phase_tower",prime=2,digits=[1,0,1])["levels"]
        self.assertEqual([r['residue'] for r in levels],['0','1','1','5'])
        self.assertEqual(levels[-1]['phase_offset_turns'],'5/8')
        for j in range(3):
            self.assertTrue(self.call(operation='phase_compatibility',prime=2,lower_depth=j,
                lower=levels[j]['residue'],upper=levels[j+1]['residue'])['compatible'])
        deep=self.call(operation='phase_tower',prime=3,digits=[2]*80)['levels'][-1]
        self.assertEqual(int(deep['modulus']),3**80)
        self.assertEqual(int(deep['residue']),3**80-1)
        for p,d in [(4,[1]),(2,[2]),(2,[.5])]:
            with self.assertRaises(RuntimeError):self.call(operation='phase_tower',prime=p,digits=d)

    def test_phase_cover_flow_preserves_coverings_and_reports_hidden_carry(self):
        initial=dict(operation='phase_cover_flow',prime=2,digits=[1,1],
                     base_phase_turns='5/6')
        result=self.call(**initial,elapsed_base_periods='1/3')
        self.assertEqual(result['base_phase_turns'],'1/6')
        self.assertEqual(result['digits'],[0,0])
        self.assertEqual(result['prefix_carry'],'1')
        self.assertEqual(result['levels'][-1]['phase_turns'],'1/24')
        # The displayed prefix returns after four base periods; the infinite
        # address does not. Preserve that distinction at the public boundary.
        wrap=self.call(**initial,elapsed_base_periods='4')
        self.assertEqual(wrap['digits'],initial['digits'])
        self.assertEqual(wrap['base_phase_turns'],initial['base_phase_turns'])
        self.assertEqual(wrap['prefix_carry'],'1')
        self.assertTrue(wrap['finite_prefix'])
        backwards=self.call(operation='phase_cover_flow',prime=3,digits=[0,0],
                            base_phase_turns='1/6',elapsed_base_periods='-1/3')
        self.assertEqual(backwards['return_count'],'-1')
        self.assertEqual(backwards['digits'],[2,2])
        self.assertEqual(backwards['prefix_carry'],'-1')
        phases=[Fraction(row['phase_turns']) for row in backwards['levels']]
        for lower,upper in zip(phases,phases[1:]):
            self.assertEqual((3*upper-lower).denominator,1)

    def test_periodic_clock_coefficients(self):
        w=dict(angular_frequency=1,coefficients=[[-1,-.3,0],[0,1,0],[1,-.3,0]])
        g=dict(angular_frequency=1,coefficients=[[-1,.5,0],[1,.5,0]])
        r=self.call(operation="periodic_clock",action="transform",clock_rate=w,
                    channels=[g],samples=256,maximum_harmonic=12)
        self.assertAlmostEqual(r["physical_period"],2*math.pi)
        coefficients={n:complex(a,b) for n,a,b in r["channels"][0]["coefficients"]}
        self.assertAlmostEqual(coefficients[0].real,-.3,places=13)
        t=self.call(operation="periodic_clock",action="coefficient_expansion",clock_rate=w,
                    series=g,harmonic=1,order=24)
        self.assertAlmostEqual(t["coefficient"][0],coefficients[1].real,places=12)
        self.assertFalse(r["certified"])

    def test_character_requires_observed_factor_depth(self):
        request=dict(operation='phase_cover_character',prime=2,base_phase_turns='1/6')
        unresolved=self.call(**request,digits=[],frequency_ratio='1/2')
        self.assertFalse(unresolved['resolved'])
        self.assertIsNone(unresolved['phase_turns'])
        self.assertEqual(unresolved['minimum_depth'],1)
        self.assertEqual(unresolved['required_additional_digits'],1)
        phases=[self.call(**request,digits=[b],frequency_ratio='1/2')['phase_turns'] for b in (0,1)]
        self.assertEqual(phases,['1/12','7/12'])
        for b in (0,1):
            self.assertEqual(self.call(**request,digits=[b],frequency_ratio='1')['phase_turns'],'1/6')
        with self.assertRaises(RuntimeError):
            self.call(**request,digits=[0],frequency_ratio='1/3')

    def test_rejects_wrong_types_and_critical_density(self):
        for request in [
            dict(operation="padic", a="1", b="0", prime=4, precision=2),
            dict(operation="graph", vertices=1.5, edges=[]),
            dict(operation="spectral_density", roots=[dict(source_density=1,
                 frequency_derivative=0, harmonic=1, coefficient=1)]),
            dict(schema="unknown/v1", operation="continued_fraction", value="1/2"),
        ]:
            with self.subTest(request=request), self.assertRaises(RuntimeError):
                self.call(**request)

    def test_hamiltonian_first_order_and_twist(self):
        torus = dict(angular_frequencies=[2], coefficients=[
            dict(index=[-1], coefficient=.5), dict(index=[1], coefficient=.5)])
        result = self.call(operation="hamiltonian_first_order", perturbation=torus,
                           action_gradient=[dict(angular_frequencies=[2], coefficients=[])],
                           hessian=[1], divisor_tolerance=1e-12)
        self.assertFalse(result["kam_certified"])
        self.assertEqual(result["mean_frequency_shift"], [0])
        value = self.call(operation="torus_fourier", action="evaluate",
                          series=result["angle_correction"][0], phases=[.7])["value"]
        self.assertAlmostEqual(value[0], -math.sin(.7)/4, places=13)
        twist = self.call(operation="hamiltonian_twist", angular_frequencies=[2], hessian=[0])
        self.assertFalse(twist["kolmogorov_numerically_nondegenerate"])
        self.assertTrue(twist["isoenergetic_numerically_nondegenerate"])
        lax = self.call(operation="inverse_scattering", action="toda_lax",
                        positions=[0, 0], momenta=[0, 0], periodic=True)
        self.assertEqual(lax["matrix"], [0, 1, 1, 0])
        self.assertEqual(lax["generator"], [0, 0, 0, 0])

    def test_invariant_readout_floor(self):
        request=dict(operation="invariant_readout_floor",
                     target_component_means=[1,3],component_weights=[1,3])
        result=self.call(**request)
        self.assertAlmostEqual(result["mean_squared_floor"],.75)
        self.assertEqual(result["best_constant"],[2.5,0])
        result=self.call(**request,allow_constant=False)
        self.assertAlmostEqual(result["mean_squared_floor"],7)

    def test_kepler_mean_anomaly_contract(self):
        # Verify the observable at apoapsis, not only individual table entries.
        result=self.call(operation="analytic_fourier",kind="kepler_cos_eccentric",
                         parameter=.6,maximum_harmonic=127,angular_frequency=2)
        self.assertIsNone(result["uniform_tail_bound"])
        self.assertEqual(result["series"]["angular_frequency"],2)
        value=self.call(operation="fourier",action="evaluate",series=result["series"],phase=math.pi)
        self.assertAlmostEqual(value["value"][0],-1,places=12)
        with self.assertRaises(RuntimeError):
            self.call(operation="analytic_fourier",kind="kepler_cos_eccentric",parameter=1,maximum_harmonic=10)


if __name__ == "__main__":
    unittest.main()
