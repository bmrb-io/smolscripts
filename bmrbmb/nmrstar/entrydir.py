#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  entrydir.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#
# 2018-09-06: bmrb entry directory is
#
# <BMRB ID>/
#          top-level README, .mol, .sdf, .png, .svg, .str
#          nmr/
#              setXX/
#                    NMRBot log file
#                    README to map experiment names to topspin directory names
#                    ... topspin files ...
#                    spectra/
#                    ... pngs ...
#                    transitions/
#                    ... peak lists ...
#

from __future__ import absolute_import

import os
import sys
import shutil
import json
import re
import tarfile

from . import README, DATASET_DIR, PEAKFILE_DIR, SPECTRA_DIR

#
#
#
class EntryDir( object ) :

    #
    #
    def __init__( self, srcdir, outdir, bmrbid, index, verbose = False ) :

        assert srcdir is not None
        assert outdir is not None
        assert bmrbid is not None

        self._srcdir = os.path.realpath( srcdir )
        if not os.path.isdir( self._srcdir ) :
            raise IOError( "Not a directory: %s" % (self._srcdir,) )
        self._outdir = os.path.realpath( outdir )
        self._bmrbid = bmrbid
        self._idx = index
        self._verbose = bool( verbose )

# we'll need sdf later
#
        self._sdf = os.path.join( self._outdir, "%s%s" % (self._bmrbid,os.path.splitext( self._idx["sdf"] )[1],) )
        os.umask( 0o002 )

    #
    #
    def makedirs( self ) :

        if not os.path.isdir( self._outdir ) :
            os.makedirs( self._outdir )


        nmrdir = os.path.split( os.path.join( self._outdir, DATASET_DIR ) )[0]
        if os.path.isdir( nmrdir ) :
            shutil.rmtree( nmrdir )

        os.makedirs( nmrdir )
        for i in self._idx["data"].keys() :
            num = int( i )
            os.makedirs( os.path.join( self._outdir, DATASET_DIR % (num,) ) )
            os.makedirs( os.path.join( self._outdir, PEAKFILE_DIR % (num,) ) )
            os.makedirs( os.path.join( self._outdir, SPECTRA_DIR  % (num,) ) )

    # copy files in top-level directory
    #
    def copymeta( self ) :

# sdf

        src = os.path.join( self._srcdir, self._idx["sdf"] )
        shutil.copy( src, self._sdf )

# svg - this one has atom labels, rename to "_nom"
# if it has a <rect> -- strip that off
#
        src = os.path.join( self._srcdir, self._idx["molpic"] )
        dst = os.path.join( self._outdir, "%s_nom%s" % (self._bmrbid,os.path.splitext( self._idx["molpic"] )[1],) )

# ... unless it's a png from marvin
#
        if os.path.splitext( self._idx["molpic"] )[1].lower() != ".svg" :
            shutil.copy( src, dst )
        else :
            rectpat = re.compile( r'^<svg:rect\s.+</svg:rect>$' )
#        dst2 = os.path.join( self._outdir, "%s_nom2%s" % (self._bmrbid,os.path.splitext( self._idx["molpic"] )[1],) )
            with open( src, "rU" ) as f :
                with open( dst, "wt" ) as out :
                    for line in f :
                        m = rectpat.search( line )
                        if m : continue
                        out.write( line )

# README
        with open( os.path.join( self._outdir, "README" ), "w" ) as out :
            out.write( README )
            out.write( "\n" )

    ################################################################################################
    # Chris had peaks and images in (spectra|transitions)/$expt_name/$expt_name.$ext
    # Keep the convention for now
    # copy peak lists
    #
    def copypeaks( self ) :
        pat = re.compile( r"\s+" )
        for i in self._idx["data"].keys() :
            set_num = int( i )
            peakdir = os.path.join( self._outdir, PEAKFILE_DIR % (set_num,) )

            for j in self._idx["data"][i]["PEAKFILES"].keys() :
                src = os.path.join( self._srcdir, j )
                dstname = pat.sub( "_", self._idx["data"][i]["PEAKFILES"][j] )
                dstdir = os.path.join( peakdir, dstname )
                dst = os.path.join( dstdir, "%s%s" % (dstname,os.path.splitext( j )[1],) )
                os.makedirs( dstdir )
                shutil.copy( src, dst )

    # copy pictures
    #
    def copyspectra( self ) :
        pat = re.compile( r"\s+" )
        for i in self._idx["data"].keys() :
            set_num = int( i )
            specdir = os.path.join( self._outdir, SPECTRA_DIR % (set_num,) )

            for j in self._idx["data"][i]["PNGFILES"].keys() :
                src = os.path.join( self._srcdir, j )
                dstname = pat.sub( "_", self._idx["data"][i]["PNGFILES"][j] )
                dstdir = os.path.join( specdir, dstname )
                dst = os.path.join( dstdir, "%s%s" % (dstname,os.path.splitext( j )[1],) )
                os.makedirs( dstdir )
                shutil.copy( src, dst )

    # untar timedomain data
    #
    def untartimedomain( self ) :

        for i in self._idx["data"].keys() :
            set_num = int( i )
            tarball = self._idx["data"][i]["timedomain"]
            src = os.path.join( self._srcdir, tarball )
            tgtdir = os.path.join( self._outdir, DATASET_DIR % (set_num,) )

# strip off top dir if present
#
            tar = tarfile.open( src, "r" )
            (car,cdr) = os.path.split( tar.getmembers()[0].name )
            while car != "" :
                (car,cdr) = os.path.split( car )
            root = "%s/" % (cdr,)
            extract = []
            for i in tar.getmembers() :
                if self._verbose :
                    sys.stdout.write( "root = %s, replace with '' in %s if startswith\n" % (root,i.name) )
                if i.name.startswith( root ) :
                    i.name = i.name.replace( root, "" )
                    extract.append( i )
            if self._verbose :
                sys.stdout.write( "Untar these to %s:\n" % (self._datadir,) )
                for i in extract : sys.stdout.write( " - %s\n" % (i.name,) )
            tar.extractall( members = extract, path = tgtdir )
            tar.close()

#
#
#
if __name__ == '__main__':

    indexfile = os.path.join( sys.argv[1], sys.argv[2] )
    with open( indexfile, "rU" ) as f :
        idx = json.load( f )
    d = EntryDir.makedir( srcdir = sys.argv[1], where = ".", bmrbid = "temp123", index = idx, verbose = True )

#
#
