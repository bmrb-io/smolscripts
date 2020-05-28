#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  nmrfam2incoming.py
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
import shutil
import pprint

if __name__ != "__main__" :
    from . import chemcomp
    from . import topspin

"""Read the json index file and copy files to $BMRBID/incoming.
This combines the "precheck" for completeness and creation of "incoming" directory
that all subsequent operations will read from."""

#
#
class Precheck( object ) :

    INDEXNAME = "__index__.jsn"

    # self._idx is read from json file
    # self._index is sanitized version
    #
    def __init__( self, indexfile, outdir, verbose = False ) :

        self._indexfile = os.path.realpath( indexfile )
        self._srcdir = os.path.split( self._indexfile )[0]
        with open( self._indexfile, "rU" ) as f :
            self._idx = json.load( f )
        if outdir is None :
            tgt = os.path.realpath( "." )
        else :
            tgt = os.path.realpath( outdir )
        self._outdir = os.path.join( tgt, "incoming" )
        self._verbose = bool( verbose )
        self._isgood = False
        self._index = {}

    #
    #
    def __del__( self ) :
        return
        if not self._isgood :
            if os.path.isdir( self._outdir ) :
                shutil.rmtree( self._outdir )

    #
    #
    @property
    def __json__( self ) :
        return json.dumps( self._index, sort_keys = True, indent = 4 )

    #
    #
    @classmethod
    def make_incoming( cls, indexfile, outdir, verbose = False ) :
        if verbose : sys.stdout.write( "%s.make_incoming()\n" % (cls.__name__,) )

        rc = cls( indexfile, outdir, verbose )
        rc.precheck()
        if verbose : 
            sys.stdout.write( rc.__json__ )
            sys.stdout.write( "\n" )
        rc.copy_files()
        if verbose : 
            sys.stdout.write( "\nindex now\n" )
            sys.stdout.write( rc.__json__ )
            sys.stdout.write( "\n" )

        return rc

    #
    #
    def _saniflush( self, name ) :
        if self._verbose : sys.stdout.write( "%s._saniflush(%s)\n" % (self.__class__.__name__,name,) )

        (car, cdr) = os.path.split( name )
        if cdr == "" :
            raise Exception( "Invalid filename %s in index file" % (name,) )
        if car == "" :
            f = os.path.realpath( os.path.join( self._srcdir, cdr ) )
        else :
            f = os.path.realpath( name )
        if not os.path.exists( f ) :
            raise IOError( "File not found: %s" % (f,) )
        return f

    # pass in the sub-tree that has necessary keys
    #
    def _check_dataset( self, data = None ) :
        if self._verbose : sys.stdout.write( "%s._check_dataset()\n" % (self.__class__.__name__,) )

        if data is None :
            data = self._idx

        rc = {}

# sample
#
        if "solvent" in data.keys() :
            rc["solvent"] = data["solvent"]
        else :
            raise Exception( "No solvent in index file" )
        if "solute_amount" in data.keys() :
            rc["solute_amount"] = data["solute_amount"]
        else :
            sys.stderr.write( "WARN: no solute amount in index file, assuming 100\n" )
            rc["solute_amount"] = 100
        if "solute_units" in data.keys() :
            rc["solute_units"] = data["solute_units"]
        else :
            sys.stderr.write( "WARN: no solute_units in index file, assuming mM\n" )
            rc["solute_units"] = "mM"

# nmrbot log is needed for experiment lists etc.
#
        if "botlog" in data.keys() :
            rc["botlog"] = self._saniflush( name = data["botlog"] )
        else :
            raise Exception( "No botlog in index file" )

# topspin tarball is needed for experiment files
#
        if "timedomain" in data.keys() :
            rc["timedomain"] = self._saniflush( name = data["timedomain"] )
        else :
            raise Exception( "No timedomain in index file" )

# CSV w/ assignments
# could actually live without
#
        if "shifts" in data.keys() :
            rc["shifts"] = self._saniflush( name = data["shifts"] )
        else :
            raise Exception( "No shifts in index file" )

# peak files
# can live w/o those
#
        peaks = {}
        if "PEAKFILES" in data.keys() :
            for f in data["PEAKFILES"].keys() :
                if not data["PEAKFILES"][f] in topspin.EXPERIMENTS.values() :
                    raise Exception( "Unknown experiment: %s for file %s" % (data["PEAKFILES"][f],f,) )
                fname = self._saniflush( name = f )
                peaks[fname] = data["PEAKFILES"][f]
        else :
            sys.stderr.write( "WARN: no PEAKFILES in index file\n" )
        rc["PEAKFILES"] = peaks

# spectra images
# ditto
#
        pix = {}
        if "PNGFILES" in data.keys() :
            for f in data["PNGFILES"].keys() :
                if not data["PNGFILES"][f] in topspin.EXPERIMENTS.values() :
                    raise Exception( "Unknown experiment: %s for file %s" % (data["PNGFILES"][f],f,) )
                fname = self._saniflush( name = f )
                pix[fname] = data["PNGFILES"][f]
        else :
            sys.stderr.write( "WARN: no PNGFILES in index file\n" )
        rc["PNGFILES"] = pix

        if len( rc ) < 1 : return None
        return rc

    # sanitize paths and make sure stuff exists
    # populate self._index
    #
    def precheck( self ) :
        if self._verbose : sys.stdout.write( "%s.precheck()\n" % (self.__class__.__name__,) )

        for i in ("sdf", "molpic") :
            if not i in self._idx.keys() :
                raise Exception( "No %s in index file" % (i,) )
            self._index[i] = self._saniflush( name = self._idx[i] ) 

# inchi string
#  (early versions of alatis didn't put it in sdf, it had to be in the index)

        inchi = None
        mol = chemcomp.Molecule.from_file( filename = self._index["sdf"], verbose = self._verbose )
        if "descriptors" in mol.chem_comp.keys() :
            for i in mol.chem_comp["descriptors"] :
                if (i["type"] == "InChI") and (i["program"] == "ALATIS") :
                    inchi = i["descriptor"]
        if inchi is None :
            if not "inchi" in self._idx.keys() :
                raise Exception( "No InChI string in SDF or index file" )
            inchi = self._idx["inchi"]
        if not inchi.startswith( "InChI=" ) :
            raise Exception( "Invalid InChI string %s" % (inchi,) )
        self._index["inchi"] = inchi

# vendor details
#
        if "catalog_name" in self._idx.keys() :
            self._index["catalog_name"] = self._idx["catalog_name"]
        else :
            raise Exception( "No catalog_name in index file" )
        if "catalog_num" in self._idx.keys() :
            self._index["catalog_num"] = self._idx["catalog_num"]
        else :
            raise Exception( "No catalog_num in index file" )
        if "vendor" in self._idx.keys() :
            self._index["vendor"] = self._idx["vendor"]
        else :
            raise Exception( "No vendor in index file" )

# TODO: convert to boolean
#
        if "paramagnetic" in self._idx.keys() :
            self._index["paramagnetic"] = self._idx["paramagnetic"]
        else :
            self._index["paramagnetic"] = "no"

# this is required to run a PUG query
#
        if "cid" in self._idx.keys() :
            self._index["cid"] = self._idx["cid"]

# older entries have no "data" key
#
        self._index["data"] = {}
        if not "data" in self._idx.keys() :
            self._index["data"][1] = self._check_dataset()
        else :
            i = 0
            for k in sorted( self._idx["data"].keys() ) :
                i += 1
                self._index["data"][i] = self._check_dataset( data = self._idx["data"][k] )

    # copy files to incoming
    #
    def copy_files( self ) :
        if self._verbose : sys.stdout.write( "%s.copy_files()\n" % (self.__class__.__name__,) )

        if os.path.exists( self._outdir ) :
            shutil.rmtree( self._outdir )
        os.makedirs( self._outdir )

# make index file in incoming w/o full paths
#
        for i in ( "sdf", "molpic" ) :
            src = self._index[i]
            fname = os.path.split( src )[1]
            dst = os.path.join( self._outdir, fname )
            if self._verbose :
                sys.stdout.write( "copy %s to %s\n" % (src,dst,) )
            shutil.copy( src, dst )
            self._index[i] = fname

        for i in self._index["data"].keys() :
            for j in ("timedomain", "botlog", "shifts") :
                src = self._index["data"][i][j]
                fname = os.path.split( src )[1]
                dst = os.path.join( self._outdir, fname )

# in case same filename is used in different data sets
#
                if os.path.exists( dst ) :
                    raise IOError( "%s already exists in %s, won't overwrite" % (fname, self._outdir) )

                if self._verbose :
                    sys.stdout.write( "copy %s to %s\n" % (src,dst,) )
                shutil.copy( src, dst )
                self._index["data"][i][j] = fname

            for j in ("PEAKFILES", "PNGFILES") :
                if not j in self._index["data"][i].keys() :
                    continue

                tmplist = {}
# filename's the key here
#
                for src in self._index["data"][i][j] :
                    fname = os.path.split( src )[1]
                    dst = os.path.join( self._outdir, fname )
                    if os.path.exists( dst ) :
                        raise IOError( "%s already exists in %s, won't overwrite" % (fname, self._outdir) )
                    if self._verbose :
                        sys.stdout.write( "copy %s to %s\n" % (src,dst,) )
                    shutil.copy( src, dst )
                    tmplist[fname] = self._index["data"][i][j][src]

                del self._index["data"][i][j]
                self._index["data"][i][j] = tmplist

# finally
#
        self._indexfile = os.path.join( self._outdir, self.INDEXNAME )
        with open( self._indexfile, "wb" ) as out :
            out.write( self.__json__ )

        self._isgood = True

#
#
#
if __name__ == "__main__" :

    import chemcomp
    import topspin

    rc = Precheck.make_incoming( indexfile = sys.argv[1], outdir = sys.argv[2], verbose = True )

#
# eof
#
