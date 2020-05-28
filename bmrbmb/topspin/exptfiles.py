#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  exptfiles.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit

from __future__ import absolute_import

import sys
import os
import shutil
import tarfile
import tempfile
import pprint
import json
import re
import glob

# untar/zip topspin exeriment directory, make a list of files in there
# these are text/directory, type "timedomain data", and named 1, 2, ...
#
# output is a dict expt.num : { "name" : same as num, "path" : ditto,
#     "content" : "NMR experiment directory",
#     "type" : "text/directory" }
#
# this will be matched to NMRBotLog's list by expt.num to fill in experiment_file table
#

class ExperimentFiles( object ) :

    UMASK = 002

    # mainline wrapper: read the tarball and make the lists but don't extract
    #
    @classmethod
    def read( cls, tarfile, outdir = None, verbose = False ) :
        rc = cls( tarfile = tarfile, outdir = outdir, verbose = verbose )
        rc._listfiles()
        rc.list()
        return rc

    #
    #
    def __init__( self, tarfile, outdir = None, verbose = False ) :

        filename = os.path.realpath( tarfile )
        if not os.path.exists( filename ) :
            raise IOError( "Not found: %s" % (filename,) )
        self._tarfile = filename

        if outdir is not None :
            self._outdir = os.path.realpath( outdir )
        else :
            self._outdir = None

        self._verbose = bool( verbose )

        self._filelist = []
        self.data = {}

    #
    #
    def __json__( self ) :
        return json.dumps( self.data )

    json = property( __json__ )

    # if this works, self._filelist will be a list of tarfile.TarInfo objects
    #
    def _listfiles( self ) :
        if self._verbose : sys.stdout.write( "%s._listfiles()\n" % (self.__class__.__name__,) )

        with tarfile.open( self._tarfile, "r" ) as tar :

# grrr...
# try to strip off the top directory
#
            (car, cdr) = os.path.split( tar.getmembers()[0].name )
            while car != "" :
                (car, cdr) = os.path.split( car )
            root = "%s/" % (cdr,)

            for i in tar.getmembers() :
                if self._verbose :
                    sys.stdout.write( "root = %s, replace with '' in %s if startswith\n" % (root,i.name) )
                if i.name.startswith( root ) :
                    i.name = i.name.replace( root, "" )
                    self._filelist.append( i )
            if self._verbose :
                sys.stdout.write( "%d files in tarball\n" % (len( self._filelist ),) )
                pprint.pprint( self._filelist )

    # fill in self.data with { "name" : N, "experiment_id" : N, "content" : "NMR experiment directory",
    #   "type" : "text/directory" }
    #
    def list( self ) :
        if self._verbose : sys.stdout.write( "%s._list()\n" % (self.__class__.__name__,) )

        if len( self._filelist ) < 1 :
            self._listfiles()

# topspin experiment directories are named 1..10
#
        dirpat = re.compile( r"^(\d+)$" )
        botpat = re.compile( r"^NMRbot_.+$" )

        for i in self._filelist :

            m = dirpat.search( i.name )
            if m :
                if i.isdir() :
                    self.data[int( i.name )] = {  "name" : i.name,
                            "experiment_id" : i.name,
                            "content" : "NMR experiment directory",
                            "type" : "text/directory" }
                    continue

# we get NMRBot log separately ATM so ignore it here for now
#
            m = botpat.search( i.name )
            if m :
                continue

# nag about wrong files?
#

        return

    # extract to a directory
    # if none is specified here, or in c'tor, make a temporary one in the cwd
    #
    def extract( self, outdir = None ) :
        if self._verbose : sys.stdout.write( "%s.extract()\n" % (self.__class__.__name__,) )

        os.umask( self.UMASK )

        if outdir is not None :
            self._outdir = os.path.realpath( outdir )

        if self._outdir is None :
            self._outdir = tempfile.mkdtemp( dir = os.path.realpath( "." ) )
        else :
            if os.path.isdir( self._outdir ) :
                shutil.rmtree( self._outdir )
            os.makedirs( self._outdir )

        if len( self._filelist ) < 1 :
            self._listfiles()

        with tarfile.open( self._tarfile, "r" ) as tar :
            tar.extractall( members = self._filelist, path = self._outdir )

##############################################
#
#
if __name__ == '__main__':

#    outdir = os.path.join( os.path.realpath( "." ), "tmp" )
    efl = ExperimentFiles( tarfile = sys.argv[1], outdir = None, verbose = True )
    efl._listfiles()
    efl.list()
    pprint.pprint( efl.data )
    efl.extract()


#
