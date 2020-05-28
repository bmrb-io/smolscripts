#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  resonance.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#

from __future__ import absolute_import
import os
import sys
import csv
import re
import pprint
import json

# This is the "better" mapper, but it returns an unsorted table of assignments.
# The caller has to sort it by chem_comp_atom.ordinal, add row numbers, and then
# make the ambiguous assignments table from those and resonance_ids
#
# This is tuned for input format w/ one atom per line but possibly multiple
# assignments and "set IDs" for ambiguous ones. E.g.
#
# atom_ID,chem_shift,ambiguity,set_ID
# C1,24.9588,1,
# C2,"63.1833, 63.6985, 63.7526",4,1
# C3,"63.1833, 63.6985, 63.7526",4,1
# C4,"63.1833, 63.6985, 63.7526",4,1
# C5,62.7662,1,
#
# - chem_shit may be na or null in which case we drop the row
#
# (set ID is redundant with this mapper, superseded by resonance id, but that's
# the format originally agreed upon)
#
class Resonance( object ) :

    # wrapper
    #
    @classmethod
    def readfile( cls, filename, verbose = False ) :
        if verbose : sys.stdout.write( "%s.readfile(%s)\n" % (cls.__name__,filename,) )
        infile = os.path.realpath( filename )
        with open( infile, "rb" ) as f :
            return cls.read( f, verbose = verbose )

    #
    #
    @classmethod
    def read( cls, fp, verbose = False ) :
        if verbose : sys.stdout.write( "%s.read()\n" % (cls.__name__,) )
        rc = cls( verbose = verbose )
        rc._read( fp )
        return rc

    #
    #
    def __init__( self, verbose = False ) :
        self._verbose = bool( verbose )
        self._dat = {}
        self._pat = re.compile( r"[\s,]+" )

    #
    #
    def _read( self, fp ) :
        if self._verbose : sys.stdout.write( "%s._read()\n" % (self.__class__.__name__,) )
        rdr = csv.DictReader( fp )
        for row in rdr :
            if row["chem_shift"] in ("","na","n/a") :
                continue

            self.add( atom = row["atom_ID"], values = row["chem_shift"], ambiguity = row["ambiguity"] )

    # input is one atom/multiple values per row
    # set of values is the same for multiple atoms
    # need to group atoms based on values
    #
    def _find( self, ambiguity, values ) :
        if self._verbose : sys.stdout.write( "%s._find()\n" % (self.__class__.__name__,) )

        ambi = int( ambiguity )
        if ambi < 2 :
            raise NotImplementedError( "ambiguous assignments only please")

        if len( values ) < 1 :
            raise Exception( "No values" )

        if (ambi in (2,3,)) and (len( values ) > 2) :
            raise Exception( "more then 2 values for ambiguity %d" % (ambi,) )

        rc = []
        for i in self._dat.keys() :

            if self._dat[i]["ambiguity"] != ambi :
                continue

            if self._dat[i]["values"] != values :
                continue

            rc.append( i )

# sanity
#
        if len( rc ) > 1 :
            raise Exception( "Found %d matches (there can be only one)" % (len( rc ),) )

        if len( rc ) < 1 :
            return -1

        return rc[0]

    # values can be many
    #
    def add( self, atom, values, ambiguity = 1, setid = None ) :
        if self._verbose : sys.stdout.write( "%s.add(%s,%s)\n" % (self.__class__.__name__,atom,ambiguity) )

        if len( self._dat.keys() ) < 1 : reson = 1
        else :
            reson = int( max( self._dat.keys() ) ) + 1

        ambi = int( ambiguity )

        vals = set()
        for v in self._pat.split( values ) :
            vals.add( float( v ) )

# one atom, one value, no set id
#
        if ambi == 1 :
            if len( vals ) != 1 :
                raise Exception( "%d values for atom %s (need 1)" % (len( vals ),atom,) )
            self._dat[reson] = { "atoms" : set( [atom] ), "values" : vals, "ambiguity" : ambi }
            return

        res = self._find( ambiguity = ambi, values = vals )
        if res < 1 :
            self._dat[reson] = { "atoms" : set( [atom] ), "values" : vals, "ambiguity" : ambi }
        else :
            self._dat[res]["atoms"].add( atom )
        return

    # ambiguity 1 means 1 value, 1 atom
    # 2 or 3 can have 1 or 2 values for 2 atoms
    # others can't have > 1 value/atom
    # but many atoms per value is OK
    #
    def _is_good( self ) :
        if self._verbose : sys.stdout.write( "%s._is_good()\n" % (self.__class__.__name__,) )

        if len( self._dat ) < 1 : return False

        rc = True
        for i in self._dat.keys() :

            if (self._dat[i]["ambiguity"] is None) \
            or (int( self._dat[i]["ambiguity"] ) < 1 ) :
                sys.stderr.write( "Invalid ambiguity code in res. %d\n" % (i,) )
                rc = False
                continue

            if len( self._dat[i]["atoms"] ) < 1 :
                sys.stderr.write( "No atoms in res. %d\n" % (i,) )
                rc = False
                continue

            if len( self._dat[i]["values"] ) < 1 :
                sys.stderr.write( "No assignments in res. %d\n" % (i,) )
                rc = False
                continue

            if self._dat[i]["ambiguity"] == 1 :
                if len( self._dat[i]["atoms"] ) != 1 :
                    sys.stderr.write( "%d atoms in res. %d (there can be only one)\n" \
                                        % (len( self._dat[i]["atoms"] ), i,) )
                    rc = False
                    continue
                if len( self._dat[i]["values"] ) != 1 :
                    sys.stderr.write( "%d assignments in res. %d (there can be only one)\n" \
                                        % (len( self._dat[i]["values"] ), i,) )
                    rc = False
                    continue

            if self._dat[i]["ambiguity"] in (2,3,) :
                if len( self._dat[i]["atoms"] ) > 2 :
                    sys.stderr.write( "%d atoms in res. %d (can be 1 or 2)\n" \
                                        % (len( self._dat[i]["atoms"] ), i,) )
                    rc = False
                    continue
                if len( self._dat[i]["values"] ) > 2 :
                    sys.stderr.write( "%d assignments in res. %d (can be 1 or 2)\n" \
                                        % (len( self._dat[i]["values"] ), i,) )
                    rc = False
                    continue

            if len( self._dat[i]["atoms"] ) < len( self._dat[i]["values"] ) :
                sys.stderr.write( "%d assignments for %d atoms in res. %d\n" \
                        % (len( self._dat[i]["values"] ), len( self._dat[i]["atoms"] ),i,) )
                sys.stderr.write( "NMR-can't-STAR\n" )
                rc = False

        return rc

    # because json can't encode sets and python can't index them
    #
    def __data__( self ) :
        rc = {}
        for i in self._dat.keys() :
            rc[i] = { "ambiguity" : self._dat[i]["ambiguity"],
                    "atoms" : list( self._dat[i]["atoms"] ),
                    "values" : list( self._dat[i]["values"] ) }

        return rc
    data = property( __data__ )

    # return NMR-STAR-ish table of assignments, 1 atom per shift
    #
    def _to_table( self ) :

        if not self._is_good() : 
            sys.stderr.write( self.__str__() )
            raise Exception( "Bad table" )

        shifts = []
        dat = self.data
        for i in dat.keys() :

            if dat[i]["ambiguity"] == 1 :
                shifts.append( { "ambiguity_code" : 1, 
                                "atom_id" : dat[i]["atoms"][0],
                                "resonance_id" : i,
                                "value" : dat[i]["values"][0] } )
                continue

# if there's only 1 value, then ambiguity code is 1
#
            if dat[i]["ambiguity"] in (2,3,) :
                if len( dat[i]["values"] ) == 1 :
                    shifts.append( { "ambiguity_code" : 1, 
                                    "atom_id" : dat[i]["atoms"][0],
                                    "resonance_id" : i,
                                    "value" :  dat[i]["values"][0] } )
                    shifts.append( { "ambiguity_code" : 1,
                                    "atom_id" : dat[i]["atoms"][1],
                                    "resonance_id" : i,
                                    "value" : dat[i]["values"][0] } )

# otherwise ambiguity code may still be 1 if the assignments are unique/stereospcefic,
# but we assume that's already correct in the input (as is the above, but the above is easy)
#
                else :
                    shifts.append( { "ambiguity_code" : dat[i]["ambiguity"], 
                                    "atom_id" : dat[i]["atoms"][0],
                                    "resonance_id" : i,
                                    "value" : dat[i]["values"][0] } )
                    shifts.append( { "ambiguity_code" : dat[i]["ambiguity"], 
                                    "atom_id" : dat[i]["atoms"][1],
                                    "resonance_id" : i,
                                    "value" : dat[i]["values"][1] } )
                continue

# else spread n assignments over m (<= n) atoms
#
            s = 0
            for a in range( len( dat[i]["atoms"] ) ) :
                if s == len( dat[i]["values"] ) :
                    s = 0
                shifts.append( { "ambiguity_code" : dat[i]["ambiguity"], 
                                "atom_id" : dat[i]["atoms"][a],
                                "resonance_id" : i,
                                "value" :  dat[i]["values"][s] } )
                s += 1

        if len( shifts ) < 1 : return None
        return shifts

    #
    #
    shifts_table = property( _to_table )

#
#
if __name__ == '__main__':

    pretty = Resonance.readfile( filename = sys.argv[1], verbose = True )
    if not pretty._is_good() :
        sys.stderr.write( "Bad file, no cookie\n" )
        pprint.pprint( pretty.data )
        sys.stderr.write( "\n" )
    else :
        sys.stdout.write( "** data **\n" )
        pprint.pprint( pretty.data )
        sys.stdout.write( "** table **\n" )
        pprint.pprint( pretty.shifts_table )
        sys.stdout.write( "\n" )

#
#
