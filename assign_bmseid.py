#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  nmrfam_to_bmrb.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit

# convert files created by NMRFAM to NMR-STAR
#

from __future__ import absolute_import

import os
import sys
import json
import re
import sqlite3
import argparse

_HERE = os.path.split( os.path.realpath( __file__ ) )[0]
sys.path.append( _HERE )
import bmrbmb

DBFILE = "bmrbids.sqlt3"

#
#
#
if __name__ == '__main__':

    ap = argparse.ArgumentParser( description = "Assign a new BMSE ID or return the existing one" )
    ap.add_argument( "-v", "--verbose", default = False, action = "store_true",
        help = "print lots of messages to stdout", dest = "verbose" )
    ap.add_argument( "-i", "--index", help = "index file", dest = "index", required = True )
    ap.add_argument( "-d", "--directory", help = "incoming directory", dest = "directory", default = None )

    args = ap.parse_args()

# DB has 
#  id 
#  dirname - used to ID the entry directory on the NMRFAM side
#  inchi - the proper key, but not convenient for teh above
#
# alatis inchi string is stored in index file
# index file is normally in the entry directory coming from NRFAM
# so we get inchi and the last directory component out of the index
#

    indexfile = os.path.realpath( args.index )
    with open( indexfile, "rU" ) as f :
        dat = json.load( f )
    if dat is None : raise Exception( "no entry data: %s" % (indexfile,) )
    indir = os.path.split( indexfile )[0]

    inchi = None
    if "inchi" in dat.keys() :
        inchi = dat["inchi"]
    else :
        if not "sdf" in dat.keys() :
            raise Exception( "No sdf in index file" )
        sdf = os.path.realpath( os.path.join( indir, dat["sdf"] ) )
        if not os.path.exists( sdf ) :
            raise IOError( "File not found: %s" % (sdf,) )
        mol = bmrbmb.chemcomp.Molecule.from_file( filename = sdf, verbose = args.verbose )
        chem_comp = mol.chem_comp
 
# inchi string is needed for web queries
# 
        if "descriptors" in chem_comp.keys() :
            for i in chem_comp["descriptors"] :
                if (i["type"] == "InChI") and (i["program"] == "ALATIS") :
                    inchi = i["descriptor"]

    if inchi is None : raise Exception( "no inchi_code in %s or SDF" % (indexfile,) )

    if args.directory is None :
        path = os.path.split( os.path.split( indexfile )[0] )[1]
        if path == "" : raise Exception( "can't figure dirname from %s" % (indexfile,) )
    else :
        path = os.path.split( args.directory )[1]

    if not os.path.exists( DBFILE ) : raise Exception( "no DB fie %s" % (DBFILE,) )

    if args.verbose :
        sys.stdout.write( "Get ID for\n InChI %s\n dir %s\n" % (inchi,path,) )

    conn = sqlite3.connect( DBFILE )
    curs = conn.cursor()
    sql = "select id from ids where inchi=?"
#    sql = "select * from ids where inchi=?"
    bmse = None
    curs.execute( sql, (inchi,) )
    for row in curs :
        bmse = row[0]

#    if bmse is None :
#        sql = "select id from ids where dirname=?"
#        curs.execute( sql, (path,) )
#        for row in curs :
#            bmse = row[0]

    if bmse is not None :
        sys.stdout.write( "%s\n" % (bmse,) )
        curs.close()
        conn.close()
        sys.exit( 0 )

    bmid = 0
    pat = re.compile( r"^bmse(\d+)$" )
    sql = "select id from ids"
    curs.execute( sql )
    for row in curs :
        m = pat.search( row[0] )
        if m :
            if int( m.group( 1 ) ) > bmid :
                bmid = int( m.group( 1 ) )

    if bmid == 0 : raise Exception( "no bmse id (huh?)" )

    bmse = "bmse%06d" % ((bmid + 1),)
    sql = "insert into ids (id,dirname,inchi) values (?,?,?)"

    curs.execute( sql, (bmse,path,inchi,) )
    if curs.rowcount != 1 : raise Exception( "rowcount is %s" % (curs.rowcount,) )
    conn.commit()
    sys.stdout.write( "BMRB ID for %s is %s\n" % (path,bmse,) )
    curs.close()
    conn.close()


#
#
