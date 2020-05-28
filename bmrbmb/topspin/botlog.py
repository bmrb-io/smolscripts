#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  botlog.py
#
#  Copyright 2018 Board of Regents University of Wisconsin - Madison
#  Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#

from __future__ import absolute_import

import sys
import os
import re
import json

_HERE = os.path.split( __file__ )[0]
sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
import bmrbmb.topspin

#
# Read NMRBot log file, extract
# - NMRBot version
# - experiment list
#   - expt name
#   - expt number
#   - spectral dimension list
#     - label
#     - nucleus
#     - sweep width
#
# Experiment "human-readable name" to NMRBot label map is in __init__.py
#
# the output is a dict
#  "version" : X,
#  "data" : { 1 : { "name" : "1D 1H",
#                            "dims" : { "F1" : { "nuc" : "1H", "sw" : 16.03 },
#  ...
#

class NMRBotLog( object ) :
    """NMRBot log reader"""

    # main
    #
    @classmethod
    def parse( cls, filename, verbose = False ) :
        """parse file, return instance of the class"""
        if verbose : sys.stdout.write( "%s.parse(%s)\n" % (cls.__name__, filename,) )

        infile = os.path.realpath( filename )
        dat = cls( verbose )

        with open( infile, "rU" ) as inf :
            expt_num = None
            for line in inf :
                if verbose :
                    sys.stdout.write( line )

                m = dat.version_pat.search( line )
                if m :
                    dat.version = m.group( 1 )
                    continue

                m = dat.expt_pat.search( line )
                if m :
                    expt_num = int( m.group( 1 ) )
                    par_set = m.group( 2 ).upper()

                    if not par_set in bmrbmb.topspin.EXPERIMENTS.keys() :
                        raise Exception( "Unknown experiment parameter set: %s" % (m.group( 2 ),) )

# adapted sweep width HSQC
#
                    if (par_set == "HSQCETGP") and (m.group( 3 ) is not None) :
                        expt_name = "2D 1H-13C HSQC SW small"
                    else :
                        expt_name = bmrbmb.topspin.EXPERIMENTS[par_set]

                    dat.data[expt_num] = { "name" : expt_name }

# next line should have experiment details
# 1 or 2D only
#

                m = dat.dim_pat.search( line )
                if m :
                    if expt_num is None :
                        raise Exception( "Experiment detail without parameter set" )

                    dims = { m.group( 1 ) : { "nuc" : m.group( 2 ), "sw" : m.group( 3 ) } }
                    if m.group( 4 ) is not None :
                        dims[m.group( 4 )] = { "nuc" : m.group( 5 ), "sw" : m.group( 6 ) }

                    dat.data[expt_num]["dims"] = dims

                    expt_num = None

        return dat

    #
    #
    def __init__( self, verbose = False ) :
        self._verbose = bool( verbose )

        self.version_pat = re.compile( r"^##\s+NMRbot\s+version\s+(\S+)\s+##" )

# e.g. "# Exp 7 - HSQCETGP (hsqcetgp) - SWadapted"
        self.expt_pat = re.compile( r"^#\s*[Ee]xp\s*(\d+)\s*-\s*(\w+)\s*\([^)]+\)(\s*-\s*SWadapted)?\s*$" )

# e.g. "-     dims %s- F1 13C, SW 221.8327; F2 1H, SW 12.0230'
        self.dim_pat = re.compile( r"^.*dims\s+.*-\s*(F\d)\s+(\d+.)(?:,\s*)SW\s+(\d*\.?\d+)(?:\s*;\s*(F\d)\s*(\d+.)\s*,\s*SW\s+(\d*\.?\d+))?\s*$" )

        self.data = {}
        self.version = None

    #
    #
    def __json__( self ) :
        return json.dumps( { "version" : self.version, "data" : self.data } )

    json = property( __json__ )

#
#
if __name__ == '__main__':

    p = NMRBotLog.parse( sys.argv[1] )
    sys.stdout.write( p.json )

#
#
