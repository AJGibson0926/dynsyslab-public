"""Native phase covers: the checked analytic query and the Viewer's phase-cover pages.

The fixtures are synthetic circles, not private research or simulated solutions.
"""
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get('DYNSYSLAB_NATIVE_BIN', ROOT / 'build/native-suite-fast/bin'))
CLI = BIN / ('dynsyslab-cli.exe' if os.name == 'nt' else 'dynsyslab-cli')
VIEWER = BIN / ('DynSysLabViewerCli.exe' if os.name == 'nt' else 'DynSysLabViewerCli')


def circle_geometry(prime=2, samples=32):
    values = [v for i in range(samples) for v in (10 + math.cos(i * math.tau / samples),
                                                  20 + math.sin(i * math.tau / samples))]
    return dict(schema='dynsyslab.phase-cover-geometry/v1', prime=prime, samples=samples,
                orbits=[dict(state_dimension=2, values=values) for _ in range(3)])


def run(*arguments):
    return subprocess.run([str(item) for item in arguments], capture_output=True, text=True,
                          encoding='utf-8', timeout=120)


@unittest.skipUnless(CLI.is_file(), 'Build the native CLI first')
class PhaseCoverQueryTest(unittest.TestCase):
    def query(self, states):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'query.json'
            path.write_text(json.dumps(dict(schema='dynsyslab.analytic-query/v1', operation='phase_cover',
                                            system=dict(name='isotonic_oscillator'), initial_states=states,
                                            prime=2, samples=64)), encoding='utf-8')
            return run(CLI, 'analytic', path)

    def test_single_level_cover_is_sampled_with_exact_labels(self):
        completed = self.query([[1.4, 0]])
        self.assertEqual(completed.returncode, 0, completed.stderr)
        geometry = json.loads(completed.stdout)
        self.assertEqual(geometry['schema'], 'dynsyslab.phase-cover-geometry/v1')
        self.assertEqual(geometry['orbits'][0]['modulus'], '1')
        self.assertEqual(len(geometry['orbits'][0]['values']), 2 * 64)

    def test_a_failed_clock_check_is_an_error(self):
        completed = self.query([[1.4, 0], [2.0, 0]])  # isochronous: no 2-cover
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn('clock', completed.stderr)


@unittest.skipUnless(VIEWER.is_file(), 'Build the native Viewer CLI first')
class PhaseCoverViewerTest(unittest.TestCase):
    def test_tower_levels_and_compatible_phase_points(self):
        for prime in (2, 3):
            with self.subTest(prime=prime), tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / 'geometry.json'
                source.write_text(json.dumps(circle_geometry(prime)), encoding='utf-8')
                completed = run(VIEWER, 'inspect', '--source', source, '--archives', directory)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                figures = json.loads(completed.stdout)['figures']
                self.assertEqual(len(figures), 3)
                self.assertEqual([len(f['traces']) for f in figures], [1 + prime ** j for j in range(3)])
                svg = Path(directory) / 'cover.svg'
                completed = run(VIEWER, 'render', '--source', source, '--archives', directory, '--output', svg)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertTrue(svg.read_text(encoding='utf-8').startswith('<svg'))


if __name__ == '__main__':
    unittest.main()
