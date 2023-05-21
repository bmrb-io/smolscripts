#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  __init__.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit


from __future__ import absolute_import

from . import boilerplate
from . import chemcomp
from . import chemshifts
from . import nmrfam
#from . import nmrstar
from . import sample
from . import topspin
from . import www
#from .nmrfam2json import FAMtoDAT
#from .nmrfam2incoming import Precheck
from .incoming2json import FAMtoJSN


#
#
#
__all__ = [
    "boilerplate",
    "chemcomp",
    "chemshifts",
    "nmrfam",
#    "nmrstar",
    "sample",
    "topspin",
    "www",
#    "Precheck",
    "FAMtoJSN"
]
#
#
