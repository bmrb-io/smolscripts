#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  bmrb.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#

from __future__ import absolute_import

import os
import sys
import json
import pprint
import collections

_HERE = os.path.split( __file__ )[0]
sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
import bmrbmb.www

# as of 2018-07-29 public chem comp database has Comp IDs as Entry IDs.
# so we can just use "get id" API
# (of course we could query the database locally too)
#
QRYURL = "http://webapi.bmrb.wisc.edu/current/search/get_id_by_tag_value/" \
        + "%(tag)s/%(value)s" \
        + "?database=chemcomps" # metabolomics, macromolecules
INCHITAG = "_Chem_comp.InChI_code"

# chem comps are supposed to be unique, there can be only one (or 0)
#
def id_from_inchi( inchi, verbose = False ) :
    global QRYURL
    global INCHITAG

    addr = QRYURL % { "tag" : INCHITAG, "value" : inchi }
    dat = bmrbmb.www.fetch_json( url = addr, verbose = verbose )
    if dat is None :
        return None
    if not isinstance( dat, collections.Iterable ) :
        sys.stderr.write( "BMRB ID from InChI: not an iterable:\n" )
        pprint.pprint( dat )
        return None
    if dat[0] is None:
        return None
    bmrbid = str( dat[0] ).strip()
    if bmrbid == "" :
        return None

    return { "db_code" : "BMRB Ligand Expo",
            "acc_code" : bmrbid,
            "acc_type" : "comp id" }

#
#
def as_json( inchi, verbose = False ) :
    dat = id_from_inchi( inchi, verbose )
    if dat is None : return None
    rc = { "chem_comp_db_link" : dat }
    return json.dumps( rc )

#
#
#
if __name__ == '__main__':
    sys.stdout.write( as_json( sys.argv[1], verbose = True ) )
    sys.stdout.write( "\n" )

#
#
