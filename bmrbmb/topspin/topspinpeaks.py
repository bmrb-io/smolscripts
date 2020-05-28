#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  topspinpeaks.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit


from __future__ import absolute_import

import sys
import os
import xml.parsers.expat
import re
import json

# wrapper for expat
# read peak list in topspin xml format.
# 1D list may have optional annotation: atom ID
#
# 1 or 2D only
# only one peak list per file!
#
# output is a dict (these are small molecules so it shouldn't be too big) of
#  experiment : NUM (should match numbers in NMRBot log)
#  peaks :
#   - either (1D)
#    num : { "shift1" : VAL, "intensity" : VAL, optional "atom" : ID }
#   - or (2D)
#    num : { "shift1" : VAL, "shift2" : VAL , "intensity" : VAL }
#
#
class TopSpinPeakList( object ) :
    """read topspin xml peak list"""

    # main
    #
    @classmethod
    def parse( cls, filename, verbose = False ) :
        """parse file, return instance of the class"""
        if verbose : sys.stdout.write( "%s.parse(%s)\n" % (cls.__name__, filename,) )

        infile = os.path.realpath( filename )
        dat = cls( verbose )
        p = xml.parsers.expat.ParserCreate()
        p.StartElementHandler = dat.start_element
        p.EndElementHandler = dat.end_element
        p.CharacterDataHandler = dat.char_data

        with open( infile, "rU" ) as inf :
            try :
                p.ParseFile( inf )
            except NotImplementedError :
                sys.stderr.write( "Filename: %s\n" % (infile,) )
                raise
        del p
        p = None

        return dat

    #
    #
    def __init__( self, verbose = False ) :
        self._verbose = bool( verbose )
        self._manual_pat = re.compile( r"^#\s+Manually\s+picked\s+peaks$" )
        self._data = {}
        self._reading_1D = False
        self._reading_2D = False
        self._manual_peaks = False
        self._check_manual = False
        self._peak_num = 0
        self._expt_num = -1
        self._have_peaks = False

    #
    #
    def __data__( self ) :
        if self._expt_num < 0 :
            return None
        if len( self._data ) < 1 :
            return None
        return { "experiment" : self._expt_num, "peaks" : self._data }

    data = property( __data__ )

    #
    #
    def __json__( self ) :
        return json.dumps( self.__data__() )

    #
    #
    json = property( __json__ )

#####################################

    #
    #
    def start_element( self, name, attrs ) :
        if self._verbose : sys.stdout.write( "%s.start_element( %s )\n" % (self.__class__.__name__, name,) )

        if name == "PeakList" :
            self._peak_num = 0
            return

        if name == "PeakList1D" :
#            if self._have_peaks :
#                sys.stderr.write( "***************************\n" )
#                sys.stderr.write( "More than one peak list in input file\n" )
#                sys.stderr.write( "This is not supported, please fix\n" )
#                sys.stderr.write( "***************************\n" )
#                raise NotImplementedError()
            self._have_peaks = True
            self._reading_1D = True
            self._reading_2D = False
            self._data.clear()
            return

        if name == "PeakList2D" :
#            if self._have_peaks :
#                sys.stderr.write( "***************************\n" )
#                sys.stderr.write( "More than one peak list in input file\n" )
#                sys.stderr.write( "This is not supported, please fix\n" )
#                sys.stderr.write( "***************************\n" )
#                raise NotImplementedError()
            self._have_peaks = True
            self._reading_1D = False
            self._reading_2D = True
            self._manual_peaks = False
            self._data.clear()
            return

# next is PeakListXDHeader that has expNo attribute we want
#
        if name in ( "PeakList1DHeader", "PeakList2DHeader" ) :
            if "expNo" in attrs.keys() :
                self._expt_num = attrs["expNo"]
            return

# for 2D only read manually picked peaks: auto-picked ones are useless
#   that is in char_data of PeakPickDetails
# for 1D PeakPickDetails does not have anything we don't already got from NMRBot log
#
        if self._reading_2D and (name == "PeakPickDetails") :
            self._check_manual = True
            return

# data is in element
#  1D peak may be annotated w/ atom ID
#
        if name == "Peak1D" :
            self._peak_num += 1
            self._data[self._peak_num] = { "shift1" : attrs["F1"], "intensity" : attrs["intensity" ] }
            if "annotation" in attrs.keys() :
                self._data[self._peak_num]["atom"] = attrs["annotation"]
            return

        if self._manual_peaks and (name == "Peak2D") :
            self._peak_num += 1
            self._data[self._peak_num] = { "shift1" : attrs["F1"], "shift2" : attrs["F2"], "intensity" : attrs["intensity" ] }

    #
    #
    def end_element( self, name ) :
        if self._verbose : sys.stdout.write( "%s.end_element( %s )\n" % (self.__class__.__name__, name,) )
        if self._reading_2D and (name == "PeakPickDetails") :
            self._check_manual = False

    # turn on the flag for reading 2D peaks
    #
    def char_data( self, buf ) :
        if self._verbose : sys.stdout.write( "%s.char_data( %s )\n" % (self.__class__.__name__, buf,) )
        if not self._check_manual :return
        m = self._manual_pat.search( buf.strip() )
        if m : self._manual_peaks = True

#
#
#
if __name__ == '__main__':
    peaks = TopSpinPeakList.parse( sys.argv[1] )
    sys.stdout.write( peaks.json )

#
#
