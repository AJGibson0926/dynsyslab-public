from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from dynsyslib import Observable
from postprocessing.koopman.spectrum import (
    DynamicalSystemReference, KoopmanEigenfunction, KoopmanEigenvalue, KoopmanSpectrum,
)
from rom import (
    BilinearControlMatrix, ConstantControlMatrix, ControlledDataset,
    ControlledTrajectory, KROM, ObservableReconstruction,
    estimate_symmetric_pulse_control_map,
    fit_bilinear_control_map_samples, fit_constant_control_map_samples,
    fit_constant_control_matrix, load_krom_archive,
    load_controlled_dataset, save_controlled_dataset,
    fit_conjugate_bilinear_control_map,
)
from postprocessing.koopman.workflows import (
    DuffingControlledTrainingSOP, generate_duffing_controlled_training,
    RPEControlledTrainingSOP, RPESymmetricPulseSOP,
    generate_rpe_controlled_training, generate_rpe_symmetric_pulse_training,
)


class _Plane:
    name="observable-plane"; initial_state=np.zeros(2); state_dimension=2
    def input_matrix(self,state): return np.asarray(((1.,),(2.,)))


def _observable_krom():
    system=_Plane(); reference=DynamicalSystemReference(None,system.name,{})
    values=[KoopmanEigenvalue.from_continuous(-.1+1j,.01,reference),
            KoopmanEigenvalue.from_continuous(-.1-1j,.01,reference)]
    functions=[KoopmanEigenfunction(values[0],"observable coordinate",
            state_callable=lambda y:complex(y[0],y[1]),gradient_callable=lambda y:np.asarray((1.,1j))),
        KoopmanEigenfunction(values[1],"observable coordinate",
            state_callable=lambda y:complex(y[0],-y[1]),gradient_callable=lambda y:np.asarray((1.,-1j)))]
    observable=Observable("quadratic-first",2,2,lambda x:np.asarray((x[0]**2,x[1])),
        lambda x:np.asarray(((2*x[0],0.),(0.,1.))),{"kind":"test"})
    return KROM(KoopmanSpectrum(reference,values,functions,[],{"method":"test"},
        dynamical_system=system),system,coordinate_observable=observable)


class _Line:
    name="line"; initial_state=np.zeros(1); state_dimension=1
    def input_matrix(self,state): return np.asarray(((2.,),))


def _line_krom():
    system=_Line(); reference=DynamicalSystemReference(None,system.name,{})
    value=KoopmanEigenvalue.from_continuous(-1.,.1,reference)
    function=KoopmanEigenfunction(value,"identity",state_callable=lambda x:complex(x[0]),
        gradient_callable=lambda x:np.asarray((1.,)))
    return KROM(KoopmanSpectrum(reference,[value],[function],[],{"method":"test"},
        dynamical_system=system),system)


class KoopmanControlRepresentationTests(unittest.TestCase):
    def test_exact_input_map_includes_observable_jacobian(self):
        model=_observable_krom()
        np.testing.assert_allclose(model.psi(np.asarray((2.,3.))),np.asarray((4.+3j,)))
        np.testing.assert_allclose(model.B_psi(np.asarray((2.,3.)),output="complex"),
            np.asarray(((4.+2j,),)))
        np.testing.assert_allclose(model.B_psi(np.asarray((2.,3.)),output="split"),
            np.asarray(((4.,),(2.,))))

    def test_all_input_representations_coexist_and_selection_is_explicit(self):
        source=_observable_krom()
        constant=ConstantControlMatrix(np.asarray(((7.+8j,),)))
        bilinear=BilinearControlMatrix(np.asarray(([[2.+0j]],)),np.asarray(((3.+4j,),)))
        model=source.with_control_representation("constant",constant)
        model=model.with_control_representation("bilinear",bilinear)
        self.assertEqual(model.default_control_representation,"exact")
        self.assertEqual(model.available_control_representations(),["exact","bilinear","constant"])
        np.testing.assert_allclose(model.B_psi(np.asarray((2.,3.)),output="complex",representation="constant"),constant.matrix)
        expected=bilinear.evaluate(model.psi(np.asarray((2.,3.))))
        np.testing.assert_allclose(model.B_psi(np.asarray((2.,3.)),output="complex",representation="bilinear"),expected)

    def test_conjugate_complete_bilinear_representation_is_general_real_linear_map(self):
        representation=BilinearControlMatrix(
            np.asarray(([[2.+1j]],)),np.asarray(((.1j,),)),
            conjugate_operators=np.asarray(([[.5-.25j]],)),
        )
        psi=np.asarray((.7-.3j,))
        expected=(2.+1j)*psi[0]+(.5-.25j)*np.conjugate(psi[0])+.1j
        np.testing.assert_allclose(representation.evaluate(psi),np.asarray(((expected,),)))
        self.assertTrue(representation.summary()["has_conjugate_operators"])

    def test_conjugate_complete_pointwise_fit_recovers_widely_linear_truth(self):
        class Function:
            representation_kind="analytic"
        class Synthetic:
            system=type("System",(),{"initial_state":np.zeros(2)})()
            input_dimension=1; lifted_dimension=1; delta_t=.01
            eigenfunctions=(Function(),)
            def psi(self,state,output="complex",**kwargs):
                value=complex(state[0],state[1])
                return np.asarray((value,))
            def B_psi(self,state,output="complex",representation="exact",**kwargs):
                value=complex(state[0],state[1])
                result=(2.+1j)*value+(.5-.25j)*np.conjugate(value)+.1j
                return np.asarray(((result,),))
        grid=np.linspace(-1.,1.,11)
        states=np.asarray([(x,y) for x in grid for y in grid])
        result=fit_conjugate_bilinear_control_map(Synthetic(),states)
        np.testing.assert_allclose(result.representation.operators,np.asarray(([[2.+1j]],)),atol=2e-14)
        np.testing.assert_allclose(result.representation.conjugate_operators,
                                   np.asarray(([[.5-.25j]],)),atol=2e-14)
        np.testing.assert_allclose(result.representation.offsets,np.asarray(((.1j,),)),atol=2e-14)
        self.assertLess(result.diagnostics["relative_residual"],1e-13)

    def test_multistep_output_regression_recovers_continuous_constant_B(self):
        model=_line_krom(); observable=Observable.component(1,0,name="x")
        reconstruction=ObservableReconstruction(observable,np.asarray(((1.+0j,),)),np.zeros(1),
            {"passed":True},{"method":"exact-test"})
        dt=.1; Ad=np.exp(-dt); G=1.-Ad; Bc=2.; controls=np.asarray((1.,-.5,.25,1.5,-1.))[:,None]
        states=[.3]
        for control in controls[:,0]: states.append(Ad*states[-1]+G*Bc*control)
        trajectory=ControlledTrajectory(np.arange(controls.shape[0]+1)*dt,
            np.asarray(states)[:,None],controls)
        result=fit_constant_control_matrix(model,ControlledDataset((trajectory,)),
            reconstruction=reconstruction,method="multistep-output")
        np.testing.assert_allclose(result.representation.matrix,np.asarray(((2.+0j,),)),rtol=1e-10,atol=1e-10)
        self.assertEqual(result.representation.time_domain,"continuous")
        self.assertIsNotNone(result.representation.discrete_matrix)

    def test_duffing_controlled_training_reconciles_source_input_scaling(self):
        policy=DuffingControlledTrainingSOP(trajectory_count=2,transition_count=3)
        source=np.asarray((((1.,),(-1.,),(.5,)),((-.25,),(.75,),(0.,))))
        initial=np.asarray(((.5,0.),(-.5,0.)))
        dataset,system=generate_duffing_controlled_training(source,
            initial_conditions=initial,sop=policy)
        self.assertEqual(dataset.transition_count,6)
        np.testing.assert_allclose(dataset.trajectories[0].inputs[:,0],(.5,-.5,.25))
        np.testing.assert_allclose(system.control_matrix,np.asarray(((0.,),(1.,))))
        self.assertEqual(dataset.metadata["source_to_dynsyslab_input_scale"],.5)
        self.assertEqual(dataset.metadata["accepted_trajectory_count"],2)

    def test_rpe_controlled_training_retains_exact_aligned_commands(self):
        policy=RPEControlledTrainingSOP(trajectory_count=4,transition_count=20,delta_tau=.001)
        dataset,system=generate_rpe_controlled_training(sop=policy)
        self.assertEqual(len(dataset.trajectories),4)
        self.assertEqual(dataset.transition_count,80)
        self.assertTrue(dataset.metadata["command_tape_exact"])
        self.assertEqual(dataset.trajectories[0].inputs.shape,(20,1))
        self.assertEqual(dataset.trajectories[0].states.shape,(21,2))
        self.assertEqual(system.control_input_units("dimensionless"),"p_control/p0")

    def test_rpe_symmetric_pulses_are_complete_opposite_pairs(self):
        policy=RPESymmetricPulseSOP(target_anchor_count=2,rectangle_anchor_count=2,
            transition_count=2,delta_tau=.001)
        dataset,_=generate_rpe_symmetric_pulse_training(sop=policy)
        self.assertEqual(len(dataset.trajectories),8)
        pairs={}
        for trajectory in dataset.trajectories:
            pairs.setdefault(trajectory.metadata["pair_id"],set()).add(trajectory.metadata["polarity"])
            self.assertTrue(trajectory.metadata["command_tape_exact"])
        self.assertEqual(set(pairs),set(range(4)))
        self.assertTrue(all(value=={-1,1} for value in pairs.values()))

    def test_symmetric_pulse_estimator_and_sample_regressions_recover_constant_map(self):
        model=_line_krom(); dt=.01; trajectories=[]
        for pair_id,anchor in enumerate((-.5,.25,.8)):
            for polarity in (-1,1):
                command=float(polarity); states=np.asarray((anchor,anchor+2.*command*dt,
                    anchor+4.*command*dt))[:,None]
                trajectories.append(ControlledTrajectory(np.asarray((0.,dt,2.*dt)),states,
                    np.full((2,1),command),{"pair_id":pair_id,"polarity":polarity}))
        samples=estimate_symmetric_pulse_control_map(model,ControlledDataset(trajectories))
        np.testing.assert_allclose(samples.matrices,np.full((3,1,1),2.+0j),atol=1e-13)
        constant=fit_constant_control_map_samples(model,samples)
        bilinear=fit_bilinear_control_map_samples(model,samples)
        np.testing.assert_allclose(constant.representation.matrix,np.asarray(((2.+0j,),)),atol=1e-13)
        np.testing.assert_allclose(bilinear.representation.offsets,np.asarray(((2.+0j,),)),atol=1e-13)
        np.testing.assert_allclose(bilinear.representation.operators,0.,atol=1e-13)

    def test_controlled_dataset_canonical_round_trip(self):
        trajectory=ControlledTrajectory(np.asarray((0.,.1,.2)),np.asarray(((1.,0.),(1.1,.2),(1.2,.1))),
            np.asarray(((.3,),(-.2,))),{"member":7})
        dataset=ControlledDataset((trajectory,),{"purpose":"round-trip"})
        with tempfile.TemporaryDirectory() as directory:
            archive,manifest=save_controlled_dataset(dataset,Path(directory)/"controlled.npz")
            restored=load_controlled_dataset(archive)
            self.assertTrue(manifest.is_file())
            self.assertEqual(restored.metadata,dataset.metadata)
            np.testing.assert_array_equal(restored.trajectories[0].states,trajectory.states)
            np.testing.assert_array_equal(restored.trajectories[0].inputs,trajectory.inputs)
            self.assertEqual(restored.trajectories[0].metadata,{"member":7})

    def test_archived_independent_bilinear_representation_remains_executable(self):
        from archives.store import serializable_system
        from dynsyslib import load_system
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); primary=root/"solution.npz"; archive=root/"inherited-krom.npz"
            system=load_system("rpe",preset="air_water_20c")
            solution_metadata={"format":"dynsyslab.solution/v1","system":"rayleigh-plesset",
                "system_parameters":serializable_system(system),"delta_t":.01,"integration_form":"dimensionless"}
            np.savez(primary,metadata_json=np.asarray(json.dumps(solution_metadata)))
            bilinear_summary={"encoding":"bilinear-koopman-input/v1","delta_t":.01,
                "provenance":{"inheritance_mode":"independent-closed-principal-restriction"}}
            analysis={"method":"koopman-with-sindy","primary_source":str(primary),
                "dictionary":{"encoding":"polynomial-dictionary/v1","degree":1,"include_constant":True},
                "state_preprocessing":{"center":[0.,0.],"scale":[1.,1.],"multiplier":[1.,1.]},
                "preprocessing":{"coordinate_system":"dimensionless"},
                "krom":{"delta_t":.01,"default_control_representation":"bilinear",
                        "control_representations":{"bilinear":bilinear_summary}}}
            retained_discrete=np.exp(np.asarray((-.1+1j,-.1-1j))*.01)
            retained_coefficients=np.asarray(((0.,1.,1j),(0.,1.,-1j)),dtype=complex)
            np.savez(archive,eigenvalues=retained_discrete,coefficients=retained_coefficients,
                training_states=np.asarray(((-1.,-1.),(0.,0.),(1.,1.))),
                bilinear_control_operators=np.asarray(([[2.]],),dtype=complex),
                bilinear_control_offsets=np.asarray(((4.,),),dtype=complex),
                metadata_json=np.asarray(json.dumps({"format":"dynsyslab.derived/v1","analysis_metadata":analysis})))
            model=load_krom_archive(archive); probe=np.asarray((1.5,-.25))
            restricted=model.control_representations["bilinear"]
            self.assertEqual(restricted.operators.shape,(1,1,1))
            np.testing.assert_allclose(model.B_psi(probe,output="complex",representation="bilinear"),
                                       np.asarray(((7.-.5j,),)))
            from rom.storage import _control_representation_arrays
            serialized=_control_representation_arrays(model)
            self.assertEqual(set(serialized),
                {"bilinear_control_operators","bilinear_control_offsets"})
            np.testing.assert_array_equal(serialized["bilinear_control_operators"],
                                          np.asarray(([[2.]],),dtype=complex))

    def test_archived_conjugate_complete_bilinear_representation_round_trips(self):
        from archives.store import serializable_system
        from dynsyslib import load_system
        from rom.storage import _control_representation_arrays
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); primary=root/"solution.npz"; archive=root/"widely-linear-krom.npz"
            system=load_system("rpe",preset="air_water_20c")
            solution_metadata={"format":"dynsyslab.solution/v1","system":"rayleigh-plesset",
                "system_parameters":serializable_system(system),"delta_t":.01,"integration_form":"dimensionless"}
            np.savez(primary,metadata_json=np.asarray(json.dumps(solution_metadata)))
            summary={"encoding":"widely-linear-bilinear-koopman-input/v1","delta_t":.01,
                "has_conjugate_operators":True}
            analysis={"method":"koopman-with-sindy","primary_source":str(primary),
                "dictionary":{"encoding":"polynomial-dictionary/v1","degree":1,"include_constant":True},
                "state_preprocessing":{"center":[0.,0.],"scale":[1.,1.],"multiplier":[1.,1.]},
                "preprocessing":{"coordinate_system":"dimensionless"},
                "krom":{"delta_t":.01,"default_control_representation":"bilinear",
                        "control_representations":{"bilinear":summary}}}
            np.savez(archive,
                eigenvalues=np.exp(np.asarray((-.1+1j,-.1-1j))*.01),
                coefficients=np.asarray(((0.,1.,1j),(0.,1.,-1j)),dtype=complex),
                training_states=np.asarray(((-1.,-1.),(0.,0.),(1.,1.))),
                bilinear_control_operators=np.asarray(([[2.]],),dtype=complex),
                bilinear_control_conjugate_operators=np.asarray(([[.5]],),dtype=complex),
                bilinear_control_offsets=np.asarray(((4.,),),dtype=complex),
                metadata_json=np.asarray(json.dumps({"format":"dynsyslab.derived/v1","analysis_metadata":analysis})))
            model=load_krom_archive(archive); probe=np.asarray((1.5,-.25))
            np.testing.assert_allclose(
                model.B_psi(probe,output="complex",representation="bilinear"),
                np.asarray(((7.75-.375j,),)),
            )
            serialized=_control_representation_arrays(model)
            self.assertEqual(set(serialized),{
                "bilinear_control_operators","bilinear_control_conjugate_operators",
                "bilinear_control_offsets",
            })


if __name__=="__main__": unittest.main()
