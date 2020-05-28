#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  make_conditions.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit


from __future__ import absolute_import

import os
import sys
import json
import pprint
import collections
import copy

#_HERE = os.path.split( __file__ )[0]
#sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
#import bmrbmb.sample
#import bmrbmb.boilerplate

# This one creates a single sveframe (dict), not a list.
#
class Conditions( object ) :

    # boilerplate
    #
    CONDITIONS = {
        "sf_category" : "sample_conditions",
        "sample_condition_variable" : [
        {
            "type" : "temperature",
            "val" : 298,
            "val_err" : 0.1,
            "val_units" : "K"
        },
        {
            "type" : "pressure",
            "val" : 1,
            "val_units" : "atm"
        },
        ]
    }

    PH = {
        "water" : {
            "type" : "pH",
            "val_err" : 0.37,
            "val_units" : "pD"
        },
        "d2o" : {
            "type" : "pD",
            "val_err" : 0.37,
            "val_units" : "pH"
        },
    }

    #
    #
    @classmethod
    def create( cls, num = 1, solvent = "d2o", pH = 7.4, verbose = False ) :
        if verbose : sys.stdout.write( "%s.create()\n" % (cls.__name__,) )

        rc = cls( num = num, solvent = solvent, pH = pH, verbose = verbose )
        rc.make_conditions()

        return rc

    #
    #
    def __init__( self, num, solvent, pH = 7.4, verbose = False ) :

        self._num = int( num )

# num = 0 actually breaks php summary page scripts... it has to be _1
#
        if self._num < 1 : self._num = 1

        self._solvent = solvent

        float( pH )
        self._pH = pH

        self._verbose = bool( verbose )

        self._dat = None

    #
    #
    @property
    def data( self ) :
        if len( self._dat ) < 1 : return None
        return self._dat

    #
    #
    def __json__( self ) :
        if len( self._dat ) < 1 : return None
        return json.dumps( self._dat )

    #
    #
    json = property( __json__ )

    #
    #
    def make_conditions( self ) :

        if self._verbose : sys.stdout.write( "%s._make_conditions()\n" % (self.__class__.__name__,) )

        self._dat = copy.deepcopy( self.CONDITIONS )
        self._dat["name"] = "sample_conditions_%s" % (self._num,)
        if self._solvent.lower() == "water" :
            self._dat["sample_condition_variable"].append( copy.deepcopy( self.PH["water"] ) )
        elif self._solvent.lower() == "d2o" :
            self._dat["sample_condition_variable"].append( copy.deepcopy( self.PH["d2o"] ) )
        if self._solvent.lower() in ("water","d2o") :
            for i in range( len( self._dat["sample_condition_variable"] ) ) :
                if self._dat["sample_condition_variable"][i]["type"] in ("pH","pD") :
                    self._dat["sample_condition_variable"][i]["val"] = self._pH

#
#
#
if __name__ == '__main__':

    d2o = Conditions.create( solvent = "d2o", verbose = True )
    pretty = json.loads( d2o.json )
    sys.stdout.write( json.dumps( pretty, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )
    ace = Conditions.create( num = 2, solvent = "acetone", verbose = True )
    pretty = json.loads( ace.json )
    sys.stdout.write( json.dumps( pretty, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )

#
#
