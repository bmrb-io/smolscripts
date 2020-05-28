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

# parameter sets in NMRBOT log to "human-readable" BMRB names
#
EXPERIMENTS = {
    "ZGCPPR"               : "1D 1H",
    "METABTOCSY"           : "2D 1H-1H TOCSY",
    "C13CPD"               : "1D 13C",
    "C13CPD_MMCD"          : "1D 13C",
    "C13DEPT90"            : "1D DEPT90",
    "C13DEPT90_MMCD"       : "1D DEPT90",
    "C13DEPT135_MMCD"      : "1D DEPT135",
    "C13DEPT135"           : "1D DEPT135",
    "HSQCETGP"             : "2D 1H-13C HSQC",
    "HSQCETGP_MMCD"        : "2D 1H-13C HSQC",
    "HSQCETGPASW"          : "2D 1H-13C HSQC SW small",
    "HMBCGPND"             : "2D 1H-13C HMBC",
    "HMBCETGPL3ND"         : "2D 1H-13C HMBC",
    "HMBCETGPL3ND_MMCD"    : "2D 1H-13C HMBC",
    "COSY45SW"             : "2D 1H-1H COSY",
    "HSQC_TOCSY_ADIA"      : "2D 1H-13C HSQC-TOCSY-ADIA",
    "HSQC_TOCSY_ADIA_MMCD" : "2D 1H-13C HSQC-TOCSY-ADIA",
}

from .botlog import NMRBotLog
from .exptfiles import ExperimentFiles
from .topspinpeaks import TopSpinPeakList

#
#
#
__all__ = [
        EXPERIMENTS,
        "NMRBotLog",
        "ExperimentFiles",
        "TopSpinPeakList"
    ]
#
#
