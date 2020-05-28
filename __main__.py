#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  nmrfam_to_bmrb.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit

#        if self._verbose : sys.stdout.write( "%s._add_peak_lists(%s,%s)\n" % (self.__class__.__name__,exptname,filename,) )

# convert files created by NMRFAM to NMR-STAR
#

from __future__ import absolute_import

import os
import sys
import ConfigParser
import argparse
import shutil
import json

import bmrbmb

CONFIG = "bmrbmb.properties"

#
#
#
if __name__ == '__main__':

    ap = argparse.ArgumentParser( description = "Generate a BMRB metabolomics entry directory" )
    ap.add_argument( "-v", "--verbose", default = False, action = "store_true",
        help = "print lots of messages to stdout", dest = "verbose" )
    ap.add_argument( "-t", "--time", help = "time the operatons", dest = "time",
        action = "store_false", default = True )

    ap.add_argument( "-c", "--config", help = "config file", dest = "conffile", required = True )
    ap.add_argument( "-b", "--bmrbid", help = "BMRB ID", dest = "bmrbid", required = True )
    ap.add_argument( "-i", "--index", help = "index file", dest = "idx", required = True )
    ap.add_argument( "-o", "--outdir", help = "output directory", dest = "outdir" )

    ap.add_argument( "--keep-json", help = "make JSON data file", dest = "json",
        action = "store_true", default = False )
#    ap.add_argument( "--no-incoming", help = "don't copy original files to 'incoming'",
#        dest = "incoming", action = "store_false", default = True )

    args = ap.parse_args()

    cp = ConfigParser.SafeConfigParser()
    f = os.path.realpath( args.conffile )
    cp.read( f )

    indexfile = os.path.realpath( args.idx )
    if not os.path.exists( indexfile ) : raise IOError( "ERR: index file not found: %s" % (indexfile,) )
    idx = None
    with open( indexfile, "rU" ) as f :
        idx = json.load( f )
    assert( idx is not None )

    incoming = os.path.split( indexfile )[0]

    out = os.path.realpath( args.outdir )
    if not os.path.isdir( out ) : raise IOError( "ERR: not a directory: %s" % (out,) )

    tgtdir = os.path.join( out, args.bmrbid )
    if os.path.exists( tgtdir ) : shutil.rmtree( tgtdir )
    os.umask( 0o002 )
    os.makedirs( tgtdir )

    inc = bmrbmb.Precheck.make_incoming( indexfile = indexfile, outdir = tgtdir, verbose = args.verbose )
    dat = bmrbmb.FAMtoJSN.make_entry( indexfile = inc._indexfile, verbose = args.verbose )

    if args.json :
        outfile = os.path.join( tgtdir, ("%s.json" % (args.bmrbid,) ) )
        with open( outfile, "wb" ) as f :
            f.write( json.dumps( dat.data, indent = 4, sort_keys = True ) )

    star = bmrbmb.nmrstar.StarMaker.from_nmrfam( config = cp, data = dat.data, id = args.bmrbid,
            verbose = args.verbose )

    outfile = os.path.join( tgtdir, ("%s.str" % (args.bmrbid,) ) )
    with open( outfile, "wb" ) as f :
        star.to_star( out = f )

# other files
#
    d = bmrbmb.nmrstar.EntryDir( srcdir = incoming, outdir = tgtdir, bmrbid = args.bmrbid,
            index = idx, verbose = args.verbose )
    d.makedirs()
    d.copymeta()
    d.copypeaks()
    d.copyspectra()
    d.untartimedomain()

    bmrbmb.nmrstar.EntryExtras.makefiles( sdf = d._sdf, verbose = args.verbose )

#
#
