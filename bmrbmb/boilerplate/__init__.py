#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  __init__.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit

# this is simply a list of "boilerplate" files in here and a method to read one in
#
#

from __future__ import absolute_import

import os
import sys
import json

_HERE = os.path.realpath( os.path.split( __file__ )[0] )

# paths are relative to _HERE
# some saveframes are dicts of one, simply so I don't have to special-case them later
# ones where sf category is the same as sf name are non-replicable saveframes
#
FILES = {
    "entry_information" : { "entry_information" : "entry_information.json" },
    "citation" : {
        "entry_citation" : "citation_entry.json",
        "citation_pubmed" : "citation_pubmed.json",
        "citation_alatis" : "citation_alatis.json",
        "citation_nmrbot" : "citation_nmrbot.json"
    },
    "natural_source" : { "natural_source" : "natural_source.json" },
    "experimental_source" : { "experimental_source" : "experimental_source.json" },
    "software" : {
        "software_alatis" : "software_alatis.json",
        "software_nmrbot" : "software_nmrbot.json",
        "software_topspin" : "software_topspin.json"
    },
    "spectrometer" :  { "spectrometer_kerry" : "spectrometer_kerry.json" },
    "chem_shift_reference" : {
        "chem_shift_reference_dss" : "cs_reference_dss.json",
        "chem_shift_reference_tms" : "cs_reference_tms.json"
    }
}

#
#
def get_saveframe( category ) :
    global _HERE
    global FILES
    assert category in FILES.keys()

    rc = []
    for name in FILES[category].keys() :
        infile = os.path.join( _HERE, FILES[category][name] )
        if not os.path.exists( infile ) :
            raise IOError( "File not found: %s" % (infile,) )
        with open( infile, "rU" ) as f :
            rc.append( json.load( f ) )

    if len( rc ) < 1 : return None
    return rc

#
#
#
__all__ = [
        "get_saveframe"
    ]
#
#
