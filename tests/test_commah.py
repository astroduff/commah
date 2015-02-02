#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_commah
----------------------------------

Tests for `commah` module.
"""

import unittest
import numpy as np
import commah


class TestCommah(unittest.TestCase):

    def setUp(self):
        pass

    def test_cosmologies(self):
        cosmolist = ['WMAP1', 'WMAP3', 'WMAP5',
                     'WMAP7', 'WMAP9', 'Planck']
        conclist = [8.84952, 6.57093, 7.66308,
                    7.893508, 8.88391, 9.25026]
        ival = 0
        for cosmo in cosmolist:
            output = commah.run(cosmo, Mi=[1e12])
            assert(np.allclose(output['c'].flatten()[0],
                   conclist[ival], rtol=1e-3))
            ival += 1
        pass

    def test_evolution(self):
        zlist = np.array([0., 1., 2.])
        conclist = np.array([7.66308, 5.70009, 4.55295])
        output = commah.run('WMAP5', zi=[0.], Mi=[1e12], z=zlist)
        assert(np.allclose(output['c'].flatten(), conclist, rtol=1e-3))
        pass

    def test_startingz(self):
        zlist = np.array([0., 1., 2.])
        conclist = np.array([4.55295, 4.43175, 4.26342])
        output = commah.run('WMAP5', zi=zlist, Mi=[1e12], z=2.)
        assert(np.allclose(output['c'].flatten(), conclist, rtol=1e-3))
        pass

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
