#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_commah
----------------------------------

Tests for `commah` module.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import commah


class TestCommah(object):
    def test_cosmologies(self):
        cosmolist = ['WMAP1', 'WMAP3', 'WMAP5',
                     'WMAP7', 'WMAP9',
                     'Planck13', 'Planck15']
        conclist = [8.84952005507872,
                    6.570929526400634,
                    7.663081978810533,
                    7.906741464320479,
                    8.883912548303279,
                    9.250262507139139,
                    9.044999285760362]
        for ival, cosmo in enumerate(cosmolist):
            output = commah.run(cosmo, Mi=[1e12])
            assert(np.allclose(output['c'].flatten()[0],
                   conclist[ival], rtol=1e-3))

    def test_evolution(self):
        zlist = np.array([0., 1., 2.])
        conclist = np.array([7.66308, 5.70009, 4.55295])
        output = commah.run('WMAP5', zi=[0.], Mi=[1e12], z=zlist)
        assert(np.allclose(output['c'].flatten(), conclist, rtol=1e-3))

    def test_startingz(self):
        zlist = np.array([0., 1., 2.])
        conclist = np.array([4.55295, 4.43175, 4.26342])
        output = commah.run('WMAP5', zi=zlist, Mi=[1e12], z=2.)
        assert(np.allclose(output['c'].flatten(), conclist, rtol=1e-3))
