# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from .__version__ import version as __version__
from . import commah
from .commah import run
from .commah import getcosmo

__all__ = ['commah', 'run', 'getcosmo']

__author__ = 'Alan Duffy'
__email__ = 'mail@alanrduffy.com'
