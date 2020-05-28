#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  pubchem.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#
# stuff we can fetch from PubChem
# see http://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access
#

from __future__ import absolute_import

import os
import sys
import json
import pprint

_HERE = os.path.split( __file__ )[0]
sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
import bmrbmb.www

# PUG View returns everything we want at once, but you have to get through layers and layers of sub-sections
#
PUGVURL = "http://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/?response_type=display"

# PUG REST needs one query per "stuff": its "full" record has names that are also in SDF but no
# synonyms, synonyms have no mesh synonyms, etc.
#
#PUGRURL = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/%s/JSON" #
#PUGSURL = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON" #


# this basically drills down through layers of headings and sub-headings. you'll have to look at
# the retrun from pug view and walk through it to make sense of this.
# it'll break if pubchem improves anything.
#
def pug_fetch( cid, url = None, verbose = False ) :

    global PUGVURL
    if url is None : addr = PUGVURL % (cid,)
    else : addr = str( url ).strip() % (cid,)

    dat = bmrbmb.www.fetch_json( url = addr, verbose = verbose )
    if dat is None : return None
    if not "Record" in dat.keys() : return None
    if not "Section" in dat["Record"].keys() : return None
    rc = {}

# each section is a dict with "TOCHeading" -- the name, and "Information" array that has values
#
    for s in dat["Record"]["Section"] :

        if not "TOCHeading" in s.keys() : continue
        if not "Section" in s.keys() : continue

        for subs in s["Section"] :

            if not "TOCHeading" in subs.keys() : continue

# at this point I don't care for "Record Title" (IUPAC name etc.) nor "Computed Descriptors" (InChI, SMILES)
#
# CAS, NSC numbers
#
            if subs["TOCHeading"] == "Other Identifiers" :
                if not "Section" in subs.keys() : continue
                for subsubs in subs["Section"] :

                    if not "TOCHeading" in subsubs.keys() : continue

                    if subsubs["TOCHeading"] == "CAS" :
                        if not "Information" in subsubs.keys() : continue
                        for val in subsubs["Information"] :
                            if ("Name" in val.keys()) and (val["Name"] == "CAS") \
                            and ("StringValue" in val.keys()) :
                                rc["cas_num"] = str( val["StringValue"] ).strip()

                    if subsubs["TOCHeading"] == "NSC Number" :
                        if not "Information" in subsubs.keys() : continue
                        for val in subsubs["Information"] :
                            if ("Name" in val.keys()) and (val["Name"] == "NSC Number") \
                            and ("StringValue" in val.keys()) :
                                rc["nsc_num"] = str( val["StringValue"] ).strip()

# synonyms: this has subsections for "Removed Synonyms" and "Depositor-Supplied Synonyms" (useless)
# and "MeSH Entry Terms" that are useful
#
            if subs["TOCHeading"] == "Synonyms" :
                if not "Section" in subs.keys() : continue
                for subsubs in subs["Section"] :
                    if not "TOCHeading" in subsubs.keys() : continue

                    if subsubs["TOCHeading"] in ("MeSH Entry Terms", "MeSH Synonyms",) :
                        if not "Information" in subsubs.keys() : continue
                        for val in subsubs["Information"] :
                            if ("Name" in val.keys()) and (val["Name"] in ("MeSH Entry Terms", "MeSH Synonyms",)) \
                            and ("StringValueList" in val.keys()) :
                                rc["mesh_terms"] = []
                                for syn in val["StringValueList"] :
                                    if all( ord( c ) < 128 for c in syn ) :
                                        rc["mesh_terms"].append( str( syn ).strip() )

    if verbose :
        sys.stdout.write( json.dumps( rc, indent = 4, sort_keys = True ) )
        sys.stdout.write( "\n" )

    return rc

# what we can get from PUG
#
def get_pubchem_data( cid, verbose = False ) :
    dat = pug_fetch( cid, verbose = verbose )
    if len( dat ) < 1 : return None

# CAS and NSC numbers, if any, go into identifiers
#
    ids = []
    if "cas_num" in dat.keys() :
        ids.append( { "identifier" : dat["cas_num"], "type" : "CAS REGISRY NUMBER", "program" : "na", "version" : "na" } )
    if "nsc_num" in dat.keys() :
        ids.append( { "identifier" : dat["nsc_num"], "type" : "NSC NUMBER", "program" : "na", "version" : "na" } )

# MeSH synonyms go into "common name" for now
#
    syns = []
    if "mesh_terms" in dat.keys() :
        for syn in dat["mesh_terms"] :
            syns.append( { "name" : syn, "type" : "MeSH term" } )

    rc = {}
    if len( ids ) > 0 : rc["identifiers"] = ids
    if len( syns ) > 0 : rc["common_name"] = syns
    if len( rc ) < 1 : return None
    return rc

# re-package into NMR-STAR-ish JSON struct, or None
#
def as_json( cid, verbose = False ) :
    rc = get_pubchem_data( cid, verbose )
    if rc is None : return None
    return json.dumps( rc )

#
#
#
if __name__ == '__main__':
    sys.stdout.write( as_json( sys.argv[1] ) )
    sys.stdout.write( "\n" )

#
#
