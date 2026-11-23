"""Display contracts for supplied phase covers; optional headless browser gate.

Set DYNSYSLAB_TEST_BROWSER to a Chromium/Edge executable to check interactions
and the bounded display clock against the exact native rational API. The
fixtures are synthetic circles, not private research or simulated solutions.
"""
from fractions import Fraction
import importlib.util
import json
import math
import os
from pathlib import Path
import unittest

from dynsyslab_client.mathematics import query
from viewer.phase_covers import phase_cover_fragment

ROOT=Path(__file__).resolve().parents[1]
CLI=Path(os.environ.get('DYNSYSLAB_NATIVE_BIN',ROOT/'build/native-suite-fast/bin'))/(
    'dynsyslab-cli.exe' if os.name=='nt' else 'dynsyslab-cli')
BROWSER=os.environ.get('DYNSYSLAB_TEST_BROWSER','')
HAS_BROWSER=bool(BROWSER and Path(BROWSER).is_file() and importlib.util.find_spec('playwright'))


def circle_geometry():
    """Minimal public geometry schema: no pre-enumerated branch tables."""
    values=[v for i in range(32) for v in (10+math.cos(i*math.tau/32),
                                          20+math.sin(i*math.tau/32))]
    return {str(p):dict(schema='dynsyslab.phase-cover-geometry/v1',prime=p,
                       samples=32,orbits=[dict(state_dimension=2,values=values) for _ in range(3)])
            for p in (2,3)}


class PhaseCoverPresentationTest(unittest.TestCase):
    def test_consumes_sampling_schema_without_branch_tables(self):
        geometry=circle_geometry();before=json.dumps(geometry)
        fragment=phase_cover_fragment(geometry,axis_labels=('position','momentum'))
        self.assertNotIn('__',fragment)
        self.assertIn('"bounds": [8.88, 11.12, 18.88, 21.12]',fragment)
        self.assertIn('position',fragment)
        self.assertEqual(json.dumps(geometry),before)

    def test_bad_display_contracts(self):
        for kwargs in (dict(root_id='bad id'),dict(bounds=[1,0,0,1]),
                       dict(bounds=[0,1,0,float('inf')]),dict(axis_labels=['only one'])):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                phase_cover_fragment(circle_geometry(),**kwargs)


@unittest.skipUnless(HAS_BROWSER and CLI.is_file(),'Set DYNSYSLAB_TEST_BROWSER and build the native CLI')
class PhaseCoverBrowserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from playwright.sync_api import sync_playwright
        cls.runtime=sync_playwright().start()
        cls.browser=cls.runtime.chromium.launch(executable_path=BROWSER,headless=True)

    @classmethod
    def tearDownClass(cls):
        cls.browser.close();cls.runtime.stop()

    def test_display_clock_against_exact_native_forward_and_backward_flow(self):
        page=self.browser.new_page()
        try:
            page.add_script_tag(content=(ROOT/'viewer/phase_cover_clock.js').read_text())
            for p in (2,3):
                for degrees,ticks in ((0,0),(60,1),(359,1),(60,40),(60,160),
                                      (0,-1),(60,-40),(359,-161)):
                    with self.subTest(p=p,degrees=degrees,ticks=ticks):
                        digits=[p-1,p-1]
                        actual=page.evaluate('args=>phaseCoverDisplayClock(...args)',[p,digits,degrees,ticks])
                        exact=query(dict(operation='phase_cover_flow',prime=p,digits=digits,
                            base_phase_turns=str(Fraction(degrees,360)),
                            elapsed_base_periods=str(Fraction(ticks,40))),executable=CLI)['result']
                        self.assertEqual(actual['digits'],exact['digits'])
                        self.assertEqual(actual['turns'],int(exact['return_count']))
                        self.assertEqual(actual['prefixCarry'],int(exact['prefix_carry']))
                        for value,level in zip(actual['phases'],exact['levels']):
                            self.assertAlmostEqual(value,float(Fraction(level['phase_turns'])),places=14)
        finally:page.close()

    def test_time_branch_prime_zoom_drag_and_state_restore(self):
        page=self.browser.new_page(viewport=dict(width=1024,height=1400));errors=[]
        page.on('pageerror',lambda error:errors.append(str(error)))
        try:
            page.set_content(phase_cover_fragment(circle_geometry(),root_id='phase-test'))
            root=page.locator('#phase-test');root.wait_for()
            state=lambda:json.loads(root.get_attribute('data-clock'))
            initial=state()
            page.locator('[data-control=advance]').click()
            self.assertEqual(state()['base'],initial['base'])
            self.assertEqual(state()['residue'],1)
            page.locator('[data-control=time]').fill('160')
            self.assertEqual(state()['residue'],0)
            self.assertEqual(state()['prefixCarry'],1)
            before=state()['phases']
            page.locator('[data-control=digit2]').select_option('1')
            self.assertEqual(state()['phases'][:2],before[:2])
            self.assertAlmostEqual((state()['phases'][2]-before[2])%1,.5)
            page.locator('[data-control=prime]').select_option('3')
            self.assertEqual(state()['residue'],4)
            self.assertEqual(len(page.locator('[data-control=digit1] option').all()),3)
            page.evaluate("window.dispatchEvent(new CustomEvent('openai:set_globals',{detail:{globals:{widgetState:{modelContent:{prime:2,phaseDegrees:30,digits:[1,0],elapsedTicks:40}}}}}))")
            self.assertEqual(state()['residue'],2)
            self.assertAlmostEqual(state()['base'],1/12)
            canvas=page.locator('.cover-canvas').first;before=state()
            pixels=canvas.evaluate('el=>el.toDataURL()')
            canvas.hover(position=dict(x=150,y=150));page.mouse.wheel(0,-120)
            page.wait_for_function("previous=>document.querySelector('.cover-canvas').toDataURL()!==previous",arg=pixels)
            self.assertEqual(state(),before)
            pixels=canvas.evaluate('el=>el.toDataURL()');box=canvas.bounding_box()
            page.mouse.move(box['x']+150,box['y']+150);page.mouse.down()
            page.mouse.move(box['x']+180,box['y']+170);page.mouse.up()
            self.assertNotEqual(canvas.evaluate('el=>el.toDataURL()'),pixels)
            self.assertEqual(state(),before)
            self.assertEqual(errors,[])
        finally:page.close()


if __name__=='__main__':unittest.main()
