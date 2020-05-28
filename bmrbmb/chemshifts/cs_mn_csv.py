#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  cs_mn_csv.py
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
import sqlite3
import pprint
import json

# Assignments come in as CSV w/ following columns:
# atom_ID,chem_shift,ambiguity,set_ID
# - chem_shift may be a comma[+space]-separated list, in that case
# - ambiguity has to be > 1
# - set_id must be not null, there may or may not be as many atoms in "a set" as there are values
#
# - chem_shit may be na or null in which case we drop the row
#
# There's 2 parts to this: check/clean-up and make NMR-STAR-ish tables.
# Final NMR-STAR tables need a few boilerplate ids etc. that aren't known here.
#
# TODO: this needs to be rewritten as a proper nomenclature mapper
#
#
class AtomChemShift( object ) :

    #
    #
    @classmethod
    def readfile( cls, filename, num = "set01", verbose = False ) :
        if verbose : sys.stdout.write( "%s.readfile(%s)\n" % (cls.__name__,filename,) )
        infile = os.path.realpath( filename )
        with open( infile, "rb" ) as f :
            return cls.read( f, num = num, verbose = verbose )

    #
    #
    @classmethod
    def read( cls, fp, num = "set01", verbose = False ) :
        if verbose : sys.stdout.write( "%s.read()\n" % (cls.__name__,) )
#        verbose = True
        rc = cls( num = num, verbose = verbose )
        rc._read( fp )
        if verbose :
            rc._dump_src()
        rc._convert()
        if verbose :
            rc._dump_dst()
        return rc

    # num is added to saveframe name
    #
    def __init__( self, num = "set01", verbose = False ) :
        self._num = num
        self._verbose = bool( verbose )
        self._conn = sqlite3.connect( ":memory:" )
        self._curs = self._conn.cursor()

# value lists are delimited w/ space and/or comma
#
        self._pat = re.compile( r"[\s,]+" )

# for metabolomics, atom IDs are Type-Number and we sort by the number.
#
        self._atompat = re.compile( r"^([A-Za-z]+)(\d+)$" )

# new stuff
#
        self._dit = {}

    # the two tables: assignments and ambiguities, in nmr-star-ish struct
    # these shouldn't get too big for a regular small molecule
    #
    def __data__( self ) :
        rc = {}

# assignments
#
        rc["atom_chem_shift"] = []
        self._curs.execute( "select id,atomid,val,ambicode from shifts order by id" )
        for row in self._curs :
            rc["atom_chem_shift"].append( {"id" : row[0], "atom_id" : row[1], "value" : row[2], "ambiguity_code" : row[3]} )

        if len( rc["atom_chem_shift"] ) < 1 :
            return None

# ambiguous assignments
#
        rc["ambiguous_atom_chem_shift"] = []
        self._curs.execute( "select id,atomid from ambis order by id, atomid" )
        for row in self._curs :
            rc["ambiguous_atom_chem_shift"].append( {"set_id" : row[0], "shift_id" : row[1]} )

# star stuff
#
        rc["name"] = "chem_shift_%s" % (self._num,)
        rc["sf_category"] = "assigned_chemical_shifts"

        return rc

    data = property( __data__ )

    #json string
    # can get big-ish, don't use on large molecules
    #
    def __json__( self ) :

        rc = self.__data__()
        if rc is None : return None
        return json.dumps( rc )

    # output
    #
    json = property( __json__ )

    # read in w/ some basic checks
    #
    def _read( self, fp ) :

        if self._verbose : sys.stdout.write( "%s._read()\n" % (self.__class__.__name__,) )

        self._curs.execute( "create table src (atomid text,val float,ambicode integer,setid integer)" )
        sql = "insert into src (atomid,val,ambicode,setid) values (:atomid,:val,:ambi,:setid)"

        rdr = csv.DictReader( fp )
        for row in rdr :
            tmp = { "atomid" : None, "val" : None, "ambi" : None, "setid" : None }
            for col in row.keys() :
                if col.lower() == "atom_id" :
                    tmp["atomid"] = str( row[col] ).strip()

                elif col.lower() == "chem_shift" :
                    tmp["val"] = str( row[col] ).strip()
                    if tmp["val"] in ("","na","n/a") :
                        tmp["val"] = None

# may be "1.23, 4.56"
#
#                    else :
#                        try :
#                            float( tmp["val"] )
#                        except ValueError :
#                            raise Exception( "Shift value for %s is not a number: %s" % (tmp["atomid"],tmp["val"],) )

                elif col.lower() == "ambiguity" :
                    tmp["ambi"] = str( row[col] ).strip()
                    if tmp["ambi"] == "" :
                        tmp["ambi"] = None
                    else :
                        try :
                            int( tmp["ambi"] )
                        except ValueError :
                            raise Exception( "Ambiguity code for %s is not a number: %s" % (tmp["atomid"],tmp["ambi"],) )

                elif col.lower() == "set_id" :
                    tmp["setid"] = str( row[col] ).strip()
                    if tmp["setid"] == "" :
                        tmp["setid"] = None

# skip rows w/o assignments
#
            if tmp["val"] is None : continue

# check ambiguous assignments
#
            if tmp["ambi"] is None :
                vals = self._pat.split( tmp["val"] )
                if len( vals ) == 1 :
                    tmp["ambi"] = 1
                else :
                    raise Exception( "Missing ambiguity code for %s" % (tmp["atomid"],) )

            if int( tmp["ambi"] ) == 1 :

# unambiguous
#
                try :
                    float( tmp["val"] )
                except ValueError :
                    raise Exception( "Shift value for %s is not a number: %s" % (tmp["atomid"],tmp["val"],) )

# codes 2 & 3 are for methylene and ring protons that don't need set ids, it's implicit
# there should be 2 assignments, and atom symbol should be "H" or "C" (there could be ambiguous methyl groups)
#
            elif (int( tmp["ambi"] ) == 2) or (int( tmp["ambi"] ) == 3) :
                m = self._atompat.search( tmp["atomid"] )
                if not m :
                    raise Exception( "Invalid atom ID %s for ambiguity %s" % (tmp["atomid"],tmp["ambi"]) )
                if not m.group( 1 ) in ("H","C") :
                    raise Exception( "Invalid atom ID %s for ambiguity %s" % (tmp["atomid"],tmp["ambi"]) )
                vals = self._pat.split( tmp["val"] )
                if len( vals ) != 2 :
                    raise Exception( "%d assignments for atom ID %s, ambiguity %s (should be 2)" % (len( vals ),tmp["atomid"],tmp["ambi"],) )
                if tmp["setid"] is not None :
                    raise Exception( "You don't need set id for %s (ambiguity code %s)" % (tmp["atomid"],tmp["ambi"],) )

# else there must be set id
#
            elif int( tmp["ambi"] ) > 3 :
                if tmp["setid"] is None :
                    raise Exception( "Missing set id for %s (ambiguity code %s)" % (tmp["atomid"],tmp["ambi"],) )

                vals = self._pat.split( tmp["val"] )
                for val in vals :
                    if len( str( val ).strip() ) < 1 :
                        continue
                    try :
                        float( val )
                    except ValueError :
                        raise Exception( "Shift value for %s is not a number: %s" % (tmp["atomid"],val,) )

# if we're still here, basic checks are OK
#
            if self._verbose :
                sys.stdout.write( sql + "\n" )
                pprint.pprint( tmp )
            self._curs.execute( sql, tmp )

    # debug
    #
    def _dump_src( self ) :
        self._curs.execute( "select atomid,val,ambicode,setid from src" )
        for row in self._curs :
            sys.stdout.write( "%s, %s, %s, %s\n" % tuple( row ) )

    # convert to NMR-STAR-ish tables
    #
    def _convert( self ) :
        if self._verbose : sys.stdout.write( "%s._convert()\n" % (self.__class__.__name__,) )

        groups = {}
        sql = "select distinct setid, ambicode from src where setid is not null and setid>0"
        self._curs.execute( sql )
        for row in self._curs :
            if row[0] in groups.keys() :
                sys.stdout.write( "ERR:\n" )
                pprint.pprint( groups )
                sys.stdout.write( "row:\n" )
                pprint.pprint( row )
                raise Exception( "duplicate set ID %s" % (row[0],) )
            groups[row[0]] = { "ambicode" : int( row[1] ), "atoms" : [], "shifts" : [] }

        sql = "select distinct atomid from src where setid=?"
        for group in groups.keys() :
            self._curs.execute( sql, (group,) )
            for row in self._curs : groups[group]["atoms"].append( row[0] )

        sql = "select distinct val from src where setid=?"
        for group in groups.keys() :
            self._curs.execute( sql, (group,) )
            for row in self._curs :
                vals = self._pat.split( str( row[0] ) )
                for val in vals :
                    if len( str( val ).strip() ) < 1 :
                        continue
                    if not val in groups[group]["shifts"] :
                        groups[group]["shifts"].append( val )

# this is just to make ordering consistent
#
        for group in groups.keys() :
            groups[group]["atoms"].sort()
            groups[group]["shifts"].sort()
        if self._verbose :
            pprint.pprint( groups )

# we may have m assignments for n atoms
# however if there's more peaks than atoms, we must be looking at mixture of isomers that has to be
# converted differently...
#
        for group in groups.keys() :
            if len( groups[group]["atoms"] ) < len( groups[group]["shifts"] ) :
                raise Exception( "%d values for %d atoms in set %d" % (len( groups[group]["shifts"] ), len( groups[group]["atoms"] ), group) )

# there could be ambiguously assigned methyl groups, like in VAL or LEU, in which case there should
# be n values and n * 3 atoms. this is assuming they didn't give it ambiguity code 2 and no set id.
#
            if len( groups[group]["atoms"] ) > len( groups[group]["shifts"] ) :
                (gr, rem) = divmod( len( groups[group]["atoms"] ), len( groups[group]["shifts"] ) )
                if (len( groups[group]["atoms"] ) % 3 == 0) and (rem == 0) :

# they may be dragons
#
                    newvals = []
                    for i in range( len( groups[group]["shifts"] ) ) :
                        for j in range( gr ) :
                            newvals.append( groups[group]["shifts"][i] )
                    groups[group]["shifts"] = newvals

# else just pad w/ 1st value: there isn't a nice way to represent this in NMR-STAR anyway
#
                else :
                    for i in range( len( groups[group]["shifts"] ), len( groups[group]["atoms"] ) ) :
                        groups[group]["shifts"].append( groups[group]["shifts"][0] )

# at this point there should be shift value for every atom
#
        if self._verbose :
            pprint.pprint( groups )

        self._curs.execute( "create table shifts (id integer,atomid text,val float,ambicode integer)" )

        sql = "insert into shifts (id,atomid,val,ambicode) values (?,?,?,?)"
        qry = "select atomid,val,ambicode from src where ambicode in (2,3)"
        curs = self._conn.cursor()
        lastval = None
        curs.execute( qry )
        for row in curs :
            m = self._atompat.search( row[0] )
            if not m :
                raise Exception( "Invalid atom ID %s" % (row[0],) )
            vals = self._pat.split( str( row[1] ) )
            if lastval == vals[0] : lastval = vals[1]
            else : lastval = vals[0]
            try :
                float( lastval )
            except :
                raise Exception( "Invalid shift value %s for atom ID %s" % (lastval,row[0],) )
            self._curs.execute( sql, (m.group( 2 ),row[0],lastval,row[2]) )
        curs.close()

#        qry = "select atomid,val,ambicode from src where setid is null or setid=0"
        qry = "select atomid,val,ambicode from src where ambicode=1"
        curs = self._conn.cursor()
        curs.execute( qry )
        for row in curs :
            m = self._atompat.search( row[0] )
            if not m :
                raise Exception( "Invalid atom ID %s" % (row[0],) )
            try :
                float( row[1] )
            except :
                raise Exception( "Invalid shift value %s for atom ID %s" % (row[1],row[0],) )
            self._curs.execute( sql, (m.group( 2 ),row[0],row[1],row[2]) )
        curs.close()

# ambiguous assignments table
#
        self._curs.execute( "create table ambis (id integer,atomid text)" )
        sql2 = "insert into ambis (id,atomid) values (?,?)"
        setid = 0
        for group in sorted( groups.keys() ) :
            setid += 1
            for i in range( len( groups[group]["atoms"] ) ) :
                m = self._atompat.search( groups[group]["atoms"][i] )
                if not m :
                    raise Exception( "Invalid atom ID %s" % (groups[group]["atoms"][i],) )
                atomid = m.group( 2 )

# this should bomb if the lists are out of sync
#
                val = groups[group]["shifts"][i]
                self._curs.execute( sql, (atomid,groups[group]["atoms"][i],val,groups[group]["ambicode"],) )
                self._curs.execute( sql2, (setid,atomid,) )

    # debug
    #
    def _dump_dst( self ) :
        self._curs.execute( "select id,atomid,val,ambicode from shifts" )
        sys.stdout.write( "------------------------\n" )
        for row in self._curs :
            sys.stdout.write( "%s, %s, %s, %s\n" % tuple( row ) )
        sys.stdout.write( "------------------------\n" )
        self._curs.execute( "select id,atomid from ambis" )
        for row in self._curs :
            sys.stdout.write( "%s, %s\n" % tuple( row ) )
        sys.stdout.write( "------------------------\n" )


#
#
if __name__ == '__main__':

    cs = AtomChemShift.readfile( sys.argv[1], verbose = True )
    pretty = json.loads( cs.json )
    sys.stdout.write( json.dumps( pretty, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )

#
#
