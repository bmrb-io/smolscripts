#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  cactvs.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#

#
# stuff we fetch from NIH
#

from __future__ import absolute_import

import os
import sys
import json
import pprint
from collections.abc import Iterable

_HERE = os.path.split( __file__ )[0]
sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
import bmrbmb.www

# cactus can output a lot of things but we only get the proper IUPAC name for now
# it tries to guess the input, better jjust use inchi to be on the safe side
#
NAMEURL = "http://cactus.nci.nih.gov/chemical/structure/%s/iupac_name"

#
#
def get_iupac_name( inchi, verbose = False ) :
    global NAMEURL
    addr = NAMEURL % (inchi,)
    dat = bmrbmb.www.fetch_text( url = addr, verbose = verbose )
    if dat is None :
        return None

    if not isinstance( dat, Iterable ):
        sys.stderr.write( "IUPAC name from InChI: not an iterable:\n" )
        pprint.pprint( dat )
        return None

    rc = []
    for i in dat :
        rc.append( { "identifier" : i,
            "type" : "IUPAC NAME",
            "program" : "http://cactus.nci.nih.gov/chemical/structure",
            "version" : "na" } )

    return rc

#
#
def as_json( inchi, verbose = False ) :
    dat = get_iupac_name( inchi, verbose )
    if dat is None : return None

    return json.dumps( { "chem_comp_identitier" : dat } )

#
#
#
if __name__ == '__main__':
    sys.stdout.write( as_json( sys.argv[1] ) )
    sys.stdout.write( "\n" )

#
#
