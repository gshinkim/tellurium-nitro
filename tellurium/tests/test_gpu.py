"""Batched control coefficients (te.controlCoefficients) must match roadrunner's own."""
import unittest

import numpy as np
import tellurium as te

MODEL = '''
J0: $X0 -> S1; k0*X0
J1: S1 -> S2; k1*S1 - k1r*S2
J2: S2 + S3 -> S4; k2*S2*S3
J3: S4 -> S3 + S2; k3*S4
J4: S2 -> $X1; k4*S2
X0=10; k0=1; k1=2; k1r=0.5; k2=0.3; k3=1.1; k4=0.7; S1=1; S2=1; S3=2; S4=1
'''


class TestControlCoefficients(unittest.TestCase):

    def setUp(self):
        self.r = te.loada(MODEL)
        self.k1 = np.array([0.5, 1.0, 2.0, 4.0])

    def expected(self, k1):
        r = te.loada(MODEL)
        r.k1 = k1
        r.steadyState()
        return (np.array(r.getScaledConcentrationControlCoefficientMatrix()),
                np.array(r.getScaledFluxControlCoefficientMatrix()))

    def check(self, res, rtol):
        self.assertTrue(res['ok'].all())
        for i, k in enumerate(self.k1):
            ccc, fcc = self.expected(k)
            np.testing.assert_allclose(res['ccc'][i], ccc, rtol=rtol, atol=rtol)
            np.testing.assert_allclose(res['fcc'][i], fcc, rtol=rtol, atol=rtol)

    def test_cpu_serial(self):
        te.gpu(False)
        self.check(te.controlCoefficients(self.r, {'k1': self.k1}), 1e-8)

    def test_cpu_processes(self):
        te.gpu(False)
        self.check(te.controlCoefficients(self.r, {'k1': self.k1}, processes=2), 1e-8)

    def test_gpu(self):
        dev = te.gpu(True)
        try:
            # MPS is float32 only
            self.check(te.controlCoefficients(self.r, {'k1': self.k1}), 1e-4 if dev == 'mps' else 1e-8)
        finally:
            te.gpu(False)

    def test_roadrunner_fallback(self):
        # without the Rust engine (roadrunner_rs) te.controlCoefficients uses roadrunner
        import sys
        from unittest import mock
        te.gpu(False)
        with mock.patch.dict(sys.modules, {'roadrunner_rs': None}):
            res = te.controlCoefficients(self.r, {'k1': self.k1})
        self.assertEqual(res['species'], list(self.r.getFloatingSpeciesIds()))
        self.check(res, 1e-8)

    def test_gpu_off_strings(self):
        self.assertEqual(te.gpu("false"), 'cpu')
        self.assertEqual(te.getGPU(), 'cpu')


if __name__ == '__main__':
    unittest.main()
