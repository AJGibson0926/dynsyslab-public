"""Exercise display bounds on a sharp turn, repeated points and backwards strides."""
import sys
import unittest
import os
import numpy as np
from resolution import compare_clocks, DisplayVertices
selector=DisplayVertices(os.environ['TRAJECTORY_RESOLUTION_LIBRARY']) if 'TRAJECTORY_RESOLUTION_LIBRARY' in os.environ else None


class ResolutionTests(unittest.TestCase):
    def test_bounds_and_extrema(self):
        if selector is None:
            self.skipTest('Supply the built shared library to test native display selection')
        t=np.linspace(-2,2,32769)
        for x,y in [(np.cos(t*12)[::-1],np.tanh(t*500)),
                    (np.sin(5*t),np.ones(len(t))),(np.zeros(1),np.zeros(1))]:
            keep=selector(x,y,1e-4)
            self.assertEqual((keep[0],keep[-1]),(0,len(x)-1))
            self.assertTrue(np.all(np.diff(keep)>0))
            for a in (x,y):
                self.assertIn(np.argmin(a),keep)
                self.assertIn(np.argmax(a),keep)
                error=abs(np.interp(np.arange(len(a)),keep,a[keep])-a)
                self.assertLessEqual(error.max(),np.ptp(a)*1e-4+1e-14)
        with self.assertRaises(ValueError):
            selector([1.,2.],[np.nan,3.])

    def test_compare_only_common_times(self):
        t=np.arange(9,dtype=float)
        y=np.column_stack([t,t**2])
        report=compare_clocks(t,y,t[1::2],y[1::2],(0.,6.))
        self.assertEqual(report['reference_count'],6)
        self.assertEqual(report['sample_count'],3)
        self.assertEqual(report['max_abs'],[0.,1.])
        self.assertEqual(report['compared_time'],[1.,6.])
        with self.assertRaises(ValueError):
            compare_clocks(t,y,t[::-1],y)


if __name__=='__main__':
    selector=DisplayVertices(sys.argv.pop(1))
    unittest.main()
