#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  make_sample.py
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

# this creates a single saveframne dict, not a list
#
#
class Sample( object ) :

# "known" sample types
# Notes: id 1 is the actual solute, to be added by make_sample()
# where type = reference, mol_common_name matches ID in boilerplate/cs_reference_*.json
#  i.e. if you add more here, also check/edit boilerplate/cs_reference_*.json
#
    SAMPLES = {
    "acetone" : {
        "solvent_system" : "acetone",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Acetone-d6",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "acetonitrile" : {
        "solvent_system" : "acetonitrile",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Acetonitrile-d3",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "benzene" : {
        "solvent_system" : "benzene",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Benzene-d6",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "chloroform" : {
        "solvent_system" : "chloroform",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Chloroform-d",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "chloro-meth" : {
        "solvent_system" : "65% chloroform/35% methanol",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Chloroform-d",
                "type" : "solvent",
                "concentration_val" : 65.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "Methanol-d4",
                "type" : "solvent",
                "concentration_val" : 35.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 4,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
   "dmso" : {
        "solvent_system" : "DMSO",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "DMSO",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "methanol" : {
        "solvent_system" : "methanol",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Methanol-d4",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "pyridine" : {
        "solvent_system" : "pyrridine",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "Pyridine-d5",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "TMS",
                "type" : "reference",
                "concentration_val" : 0.05,
                "concentration_val_units" : "%"
            }
        ]
    },
    "d2o" : {
        "solvent_system" : "100% D2O",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "D2O",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "DSS",
                "type" : "reference",
                "concentration_val" : 0.01,
                "concentration_val_units" : "mg/mL"
            },
            {
                "id" : 4,
                "mol_common_name" : "sodium phosphate",
                "type" : "buffer",
                "concentration_val" : 50,
                "concentration_val_units" : "mM"
            },
            {
                "id" : 5,
                "mol_common_name" : "sodium azide",
                "type" : "cytocide",
                "concentration_val" : 500,
                "concentration_val_units" : "uM"
            },
        ]
    },
    "water" : {
        "solvent_system" : "100% H2O",
        "sample_component" : [
            {
                "id" : 2,
                "mol_common_name" : "H2O",
                "type" : "solvent",
                "concentration_val" : 100.0,
                "concentration_val_units" : "%"
            },
            {
                "id" : 3,
                "mol_common_name" : "DSS",
                "type" : "reference",
                "concentration_val" : 0.01,
                "concentration_val_units" : "mg/mL"
            },
            {
                "id" : 4,
                "mol_common_name" : "sodium phosphate",
                "type" : "buffer",
                "concentration_val" : 50,
                "concentration_val_units" : "mM"
            },
            {
                "id" : 5,
                "mol_common_name" : "sodium azide",
                "type" : "cytocide",
                "concentration_val" : 500,
                "concentration_val_units" : "uM"
            },
        ]
    }
    }

    # boilerplate
    #
    SAMPLE = {
        "sf_category" : "sample",
        "type" : "solution"
    }

    # solvent: if it's H2O or D2O we add pH/pD row to sample conditions
    # pH: ignored if solvent is not one of the above
    # num: if it's 0, there is only one saveframe and it gets the same name as category
    #      otherwise the name is category_num
    # solute and catalog_name should be the same, usually
    #
    @classmethod
    def create( cls, solute, solvent, vendor, catalog_name, catalog_num, concentration = 100,
            units = "mM", num = 0, verbose = False ) :
        if verbose : sys.stdout.write( "%s.create()\n" % (cls.__name__,) )

        rc = cls( num = num, solvent = solvent, verbose = verbose )
        rc.make_sample( vendor = vendor, catalog_name = catalog_name, catalog_num = catalog_num,
                solute = solute, concentration = concentration, units = units )
        return rc

    #
    #
    def __init__( self, num, solvent, verbose = False ) :

        self._num = int( num )

# num = 0 actually breaks php summary page scripts... it has to be _1
#
        if self._num < 1 : self._num = 1

        self._solvent = str( solvent ).strip().lower()
        if self._solvent == "h2o" : self._solvent = "water"
        assert self._solvent in self.SAMPLES.keys()

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
    def make_sample( self, vendor, catalog_name, catalog_num, solute, concentration = 100, units = "mM" ) :

        if self._verbose : sys.stdout.write( "%s._make_sample()\n" % (self.__class__.__name__,) )

        self._dat = copy.deepcopy( self.SAMPLE )
        if self._num < 1 : self._dat["name"] = "sample"
        else : self._dat["name"] = "sample_%s" % (self._num,)
        self._dat["solvent_system"] = self.SAMPLES[self._solvent]["solvent_system"]
        self._dat["sample_component"] = copy.deepcopy( self.SAMPLES[self._solvent]["sample_component"] )

# component #1
#
        ein = {
            "id" : 1,
            "mol_common_name" : solute,
            "type" : "solute",
            "concentration_val" : concentration,
            "concentration_val_units" : units,
            "isotopic_labeling" : "natural abundance",
        }
        if vendor is not None : ein["vendor"] = vendor
        if catalog_num is not None : ein["vendor_product_code"] = catalog_num
        if catalog_name is not None : ein["vendor_product_name"] = catalog_name
        self._dat["sample_component"].append( ein )

#
#
#
if __name__ == '__main__':

    smpl = Sample.create( solute = sys.argv[1], solvent = sys.argv[2], vendor = "Sigma", 
        catalog_name = sys.argv[1], catalog_num = "123", concentration = 100, units = "mM", verbose = True )

    pretty = json.loads( smpl.json )
    sys.stdout.write( json.dumps( pretty, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )

#
#
