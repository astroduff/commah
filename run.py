#!/usr/bin/env python

from __future__ import absolute_import, division, print_function
import commah

# Here's an example run script, modify as preferred
commah.run('WMAP5',
           zi=0,
           Mi=[1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14],
           z=[0, 0.5, 1, 1.5, 2],
           filename='WMAP5_Test.txt')
