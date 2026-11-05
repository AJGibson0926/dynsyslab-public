from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from postprocessing.koopman import (DynamicalSystemReference, FreeTrajectoryDataset,
                                    KoopmanEigenfunction, KoopmanEigenvalue,
                                    KoopmanSpectrum)
from rom import (ConstantControlMatrix, ControlledDataset, ControlledTrajectory,
                 ControlRegressionResult,
                 KROM, load_krom_error_evaluations,
                 save_krom_error_evaluations)
from time_series import TimeSeries


class _Spiral:
    name="validation-spiral"; initial_state=np.asarray((0.,0.)); state_dimension=2
    rate=-.2; frequency=1.
    def rhs(self,time,state):
        x,y=np.asarray(state,dtype=float)
        return np.asarray((self.rate*x-self.frequency*y,self.frequency*x+self.rate*y))
    def input_matrix(self,state): return np.asarray(((1.,),(2.,)))


def _krom():
    system=_Spiral(); reference=DynamicalSystemReference(None,system.name,{})
    positive=KoopmanEigenvalue.from_continuous(system.rate+1j*system.frequency,.01,reference)
    negative=KoopmanEigenvalue.from_continuous(system.rate-1j*system.frequency,.01,reference)
    functions=[KoopmanEigenfunction(positive,"linear",state_callable=lambda x:complex(x[0],x[1]),
            gradient_callable=lambda x:np.asarray((1.,1j))),
        KoopmanEigenfunction(negative,"linear",state_callable=lambda x:complex(x[0],-x[1]),
            gradient_callable=lambda x:np.asarray((1.,-1j)))]
    spectrum=KoopmanSpectrum(reference,[positive,negative],functions,[],{
        "method":"exact-test","training_domain":{"encoding":"axis-aligned-training-domain/v1",
        "minimum":[-1.,-1.],"maximum":[1.,1.],"sample_count":1,
        "coordinate_system":"model"}},dynamical_system=system)
    return KROM(spectrum,system)


def _trajectory(initial=(.5,.2),duration=2.,step=.01):
    time=np.arange(0.,duration+step/2,step); z=complex(*initial)*np.exp(complex(-.2,1.)*time)
    return TimeSeries(time,np.column_stack((z.real,z.imag)))


def _controlled_dataset(label):
    time=np.asarray((0.,.01,.02)); offset=0. if label=="rich" else .1
    states=np.asarray(((.2+offset,.1),(.21+offset,.11),(.22+offset,.12)))
    return ControlledDataset((ControlledTrajectory(time,states,np.asarray(((1.,),(-1.,))),{"member":0}),),
        {"label":label})


class KoopmanErrorAnalysisTests(unittest.TestCase):
    def test_eigenfunction_records_training_semigroup_and_generator_errors(self):
        model=_krom(); dataset=FreeTrajectoryDataset((_trajectory(),),label="free training",
            role="training",coordinate_system="model",provenance={"source":"synthetic exact"})
        evaluation=model.evaluate_eigenfunction_errors(dataset,horizons=(.1,.5))[0]
        self.assertLess(evaluation.one_step.relative_l2,2e-13)
        self.assertLess(evaluation.horizons[0].statistics.relative_l2,2e-13)
        self.assertLess(evaluation.anchored_long_term.relative_l2,2e-13)
        self.assertLess(evaluation.generator.relative_l2,5e-5)
        self.assertEqual(evaluation.dataset["role"],"training")
        self.assertEqual(evaluation.coverage["outside_fraction"],0.)
        self.assertIn(evaluation.evaluation_id,model.eigenfunctions[0].error_evaluations)

    def test_test_dataset_reports_extrapolation_outside_training_domain(self):
        model=_krom(); dataset=FreeTrajectoryDataset((_trajectory((2.,0.),duration=.2),),
            label="outside test",role="test")
        evaluation=model.evaluate_eigenfunction_errors(dataset,horizons=(.1,))[0]
        self.assertEqual(evaluation.dataset["role"],"test")
        self.assertGreater(evaluation.coverage["outside_fraction"],0.)
        self.assertGreater(evaluation.coverage["maximum_relative_excursion"],0.)

    def test_multiple_control_datasets_remain_attached_and_comparable(self):
        model=_krom(); exact=ConstantControlMatrix(np.asarray(((1.+2j,),)),provenance={"fit":"good"})
        poor=ConstantControlMatrix(np.asarray(((0.+0j,),)),provenance={"fit":"poor"})
        good=model.evaluate_control_approximation(exact,_controlled_dataset("rich"),
            ((-1.,-1.),(1.,1.)),label="rich-data constant",grid_shape=(9,9))
        bad=model.evaluate_control_approximation(poor,_controlled_dataset("weak"),
            ((-1.,-1.),(1.,1.)),label="weak-data constant",grid_shape=(9,9))
        self.assertEqual(len(model.control_error_evaluations),2)
        self.assertNotEqual(good.controlled_dataset["fingerprint_sha256"],
            bad.controlled_dataset["fingerprint_sha256"])
        comparison=model.compare_control_error_evaluations(metric="relative.mean")
        self.assertEqual(comparison["ranking"][0]["evaluation_id"],good.evaluation_id)
        self.assertLess(comparison["ranking"][0]["value"],comparison["ranking"][1]["value"])

    def test_regression_result_has_direct_formal_evaluation_bridge(self):
        model=_krom(); dataset=_controlled_dataset("rich")
        result=ControlRegressionResult(ConstantControlMatrix(np.asarray(((1.+2j,),))),
            {"method":"test-regression"})
        evaluation=result.evaluate(model,dataset,((-1.,-1.),(1.,1.)),label="direct bridge",
            grid_shape=(5,5))
        self.assertEqual(evaluation.fit_diagnostics["method"],"test-regression")
        self.assertIn(evaluation.evaluation_id,model.control_error_evaluations)

    def test_krom_error_bundle_round_trips_free_metrics_fits_and_heatmaps(self):
        model=_krom(); dataset=FreeTrajectoryDataset((_trajectory(),),label="training",role="training")
        model.evaluate_eigenfunction_errors(dataset,horizons=(.1,))
        model.evaluate_control_approximation(ConstantControlMatrix(np.asarray(((1.+2j,),))),
            _controlled_dataset("controlled"),((-1.,-1.),(1.,1.)),label="constant",
            grid_shape=(7,8))
        with tempfile.TemporaryDirectory() as directory:
            archive=save_krom_error_evaluations(Path(directory)/"errors.npz",model)
            restored=_krom(); result=load_krom_error_evaluations(archive,restored)
        self.assertEqual(len(restored.eigenfunctions[0].error_evaluations),1)
        self.assertEqual(len(restored.control_error_evaluations),1)
        evaluation=next(iter(restored.control_error_evaluations.values()))
        self.assertEqual(evaluation.error_field.relative_error.shape,(7,8))
        self.assertEqual(result["encoding"],"dynsyslab.krom-error-evaluations/v1")


if __name__=="__main__": unittest.main()
