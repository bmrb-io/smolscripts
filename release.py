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
import argparse
import ConfigParser

_EGG = os.path.join( os.path.split( os.path.realpath( __file__ ) )[0], "validate-1.0-py2.7.egg" )
sys.path.append( _EGG )
import validator
import validator.functions

import pprint

#
#
#
def release_entry( filename, config, debug = 0 ) :
    errs = []
    db = validator.starobj.DbWrapper( config = cp, verbose = ( debug & 1 ) )
    db.connect()

    sd = validator.starobj.StarDictionary( db, verbose = ( debug & 2 ) )
    se = validator.starobj.NMRSTAREntry( db, verbose = ( debug & 4 ) )

    p = validator.starobj.StarParser.parse_file( db = db, dictionary = sd,
                    filename = filename, errlist = errs, types = False,
                    verbose = (debug & 8) )
    if len( errs ) > 0 :
        sys.stderr.write( "--------------- parse errors -------------------\n" )
        for e in errs :
            sys.stderr.write( str( e ) )
            sys.stderr.write( "\n" )

        del errs[:]
        return False

# check
#
    dirty = False

# saveframes
#
    fn = validator.functions.CheckSaveFrames.run( starentry = se, stardictionary = sd,
            errorlist = errs, verbose = (debug & 16) )
    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "check saveframes: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

        del errs[:]

# edits
# chicken-and-egg: need to add the missing stuff before checking what's missing.
#
# The function needs ETS to fetch dates for macromolecule entries, but here the dates are already set.
# Just fake it to keep it happy
#

    dsn = { "fake" : "invalid" }

    fn = validator.functions.FixEntryInfo( etsdsn = dsn, starentry = se,
            stardictionary = sd, errorlist = errs, verbose = (debug & 16) )

    fn.release_type = "original"
    fn._release_author = "BMRB"
    fn._fetch_ids()
    fn.fix_entry()

    relnum = fn.last_release_number()
#    sys.stdout.write( "** %s\n" % (relnum,) )
    if relnum == 0 :
        itr = fn._entry.iter_values( "Entry", columns = ("Original_release_date",), sfid = fn._entrysfid )
        i = itr.next()
#        pprint.pprint( i )
        if i is not None :
            fn._release_date = i[0]
        fn.insert_release_row()

    fn.count_data_sets()
    fn.count_datums()
    fn.insert_entryid()

    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "fix entry information: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

    del errs[:]

# tables
#
    fn = validator.functions.CheckTables.run( starentry = se, stardictionary = sd,
            errorlist = errs, verbose = (debug & 16) )
    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "check tables: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

        del errs[:]

# now that framecodes and local ids are OK,
#
    fn = validator.functions.FixChildTags.run( starentry = se,
            stardictionary = sd, errorlist = errs, verbose = (debug & 16) )

    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "fix parent-child: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

    del errs[:]

# DDL
#
    fn = validator.functions.CheckDDLTypes.run( starentry = se, stardictionary = sd,
                errorlist = errs, verbose = (debug & 16) )
    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "check DDL types: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

        del errs[:]

# DB
#
    fn = validator.functions.CheckDBTypes.run( starentry = se, stardictionary = sd,
            errorlist = errs, verbose = (debug & 16) )
    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "check DB types: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

        del errs[:]

# NOT NULLs
#
    fn = validator.functions.CheckNulls.run( starentry = se, stardictionary = sd,
            errorlist = errs, verbose = (debug & 16) )
    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "check NULLs: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )

    del errs[:]

# PKs
#
    fn = validator.functions.CheckPrimaryKeys.run( starentry = se, stardictionary = sd,
            errorlist = errs, verbose = (debug & 16) )
    dirty = (dirty or fn.dirty)

    if len( errs ) > 0 :
        sys.stdout.write( "check PKs: %d errors\n" % (fn.error_count,) )
        for e in errs :
            sys.stdout.write( str( e ) )
            sys.stdout.write( "\n" )


    if dirty :

# write
#
        sys.stdout.write( "File was modified\n" )
        out = open( os.path.realpath( filename ), "wb" )
        p = validator.starobj.StarWriter.pretty_print( entry = se, dictionary = sd, out = out,
                errlist = errs, verbose = (debug & 32) )
        out.close()
        if len( errs ) > 0 :
            sys.stdout.write( "-------- unparse errors -----------\n" )
            for e in errs :
                sys.stdout.write( str( e ) )
                sys.stdout.write( "\n" )


#
#
#
if __name__ == '__main__':

    ap = argparse.ArgumentParser( description = "Assign a new BMSE ID or return the existing one" )
    ap.add_argument( "-c", "--conffile", help = "config file", dest = "conffile", required = True )
    ap.add_argument( "-i", "--infile", help = "input file", dest = "infile", required = True )
    ap.add_argument( "-d", "--debug", help = "debug flag", dest = "debug", type = int, default = 0 )

    args = ap.parse_args()

    infile = os.path.realpath( args.infile )
    if not os.path.exists( infile ) :
        sys.stderr.write( "Not found: %s" % (infile,) )
        sys.exit( 1 )

    conffile = os.path.realpath( args.conffile )
    if not os.path.exists( conffile ) :
        sys.stderr.write( "Not found: %s" % (conffile,) )
        sys.exit( 2 )

    cp = ConfigParser.SafeConfigParser()
    cp.read( conffile )

    release_entry( filename = infile, config = cp, debug = args.debug )

#
#
