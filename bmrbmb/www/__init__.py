#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  __init__.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#

from __future__ import absolute_import

import os
import sys
import urllib.request
import urllib.error
import traceback
import json

# url fetch wrapper
#
def fetch_json( url, verbose = False ) :
    if url is None : return None
    rc = None
    try :
        f = urllib.request.urlopen( url )
        if (f.getcode() is not None) and (int( f.getcode() ) < 300) :
            try :
                rc = json.load( f )
            except Exception as e :
                for i in f : sys.stderr.write( "* " + i + "\n" )
                traceback.print_exc()
                rc = None
        f.close()
    except urllib.error.HTTPError as e :
        sys.stderr.write( "Error: %s at %s\n" %(e.code, url,) )
        if e.code not in (400, 404,) :
            raise
    if (rc is None) or (len( rc ) < 1) : return None
    return rc
#

# returns lines of text/plain response in an array
#
def fetch_text( url, verbose = False ) :
    if url is None : return None
    rc = []
    try :
        f = urllib.request.urlopen( url )
        if (f.getcode() is not None) and (int( f.getcode() ) < 300) :
            if f.info().gettype() != "text/plain" :
                sys.stderr.write( "Not text/plain from %s\n" % (f.geturl(),) )
            try :
                for i in f :
                    s = str( i ).strip()
                    if s != "" :
                        rc.append( s )
            except Exception as e :
                for i in f : sys.stderr.write( "* " + i + "\n" )
                traceback.print_exc()
                rc = []
        f.close()
    except urllib.error.HTTPError as e :
        if e.code != 404 :
            raise
    if len( rc ) < 1 : return None
    return rc

#

from .bmrb import as_json as bmrbid_json, id_from_inchi as bmrbid
from .cactvs import as_json as cactus_iupac_name_json, get_iupac_name as cactus_iupac_name
from .pubchem import as_json as pubchem_json, get_pubchem_data

#
#
#
__all__ = [ "fetch_json", "fetch_text",
    "bmrbid", "bmrbid_json",
    "cactus_iupac_name", "cactus_iupac_name_json",
    "get_pubchem_data", "pubchem_json",
    ]

#
#
