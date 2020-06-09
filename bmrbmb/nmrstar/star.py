#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  star.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#
# make a bmrb entry out of a "from nmrfam" data structure and a bunch of extras
#

from __future__ import absolute_import

import os
import sys
import ConfigParser
import json
import datetime
import pprint
import re

STAROBJ_PATH = "/share/dmaziuk/projects/starobj"
sys.path.append( STAROBJ_PATH )
import starobj

# self
#
#_HERE = os.path.split( __file__ )[0]
#sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
#import bmrbmb.nmrstar
from . import DATASET_MASK, DATASET_DIR, PEAKFILE_DIR, SPECTRA_DIR
#
#
#
class StarMaker( object ) :

    # our DOI: 10.12018/BMSE12345
    #
    DOI = "10.12018/%s"

    #
    #
    @classmethod
    def from_nmrfam( cls, config, data, id = "temp123", verbose = False ) :

        wrp = starobj.DbWrapper( config = config, verbose = verbose )
        wrp.connect()

        sd = starobj.StarDictionary( wrp, verbose = verbose )
        sd.printable_tags_only = True
        sd.public_tags_only = True

        star = starobj.NMRSTAREntry( wrp, verbose = verbose )

        rc = cls( stardict = sd, entrydb = star, entrydata = data, entryid = id,
                verbose = verbose )
        rc._create_tables()

# hardcoded
# used by _make_entry_info() and _make_chem_comp()
# FIXME oughta run enrtydir() first abd make sure these actually exist...
#
        files = { "MOL" : "%s.mol" % (rc._id,), "PNG" : "%s.png" % (rc._id,),
            "SVG" : "%s.svg" % (rc._id,), "SVG_N" : "%s_nom.svg" % (rc._id,),
            "SDF" : "%s.sdf" % (rc._id,) }

# there's only one entry information
#
        tables = rc._get_tables( "entry_information" )
        if tables is None : raise Exception( "no entry information" )
        rc._make_entry_info( table = tables[0], files = files )

# we want entry citation to come 1st in the star file
#
        tables = rc._get_tables( "citations" )
        if tables is None : raise Exception( "no citations" )
        num = 2
        for table in tables :
            if table["class"] == "entry citation" :
                rc._make_citation( table, citid = 1 )
            else :
                rc._make_citation( table, citid = num )
                num += 1

# NOTE that there is only one natural/experimental source saveframes,
# in these entries there's only one row in their data tables and almost all values are 'na',
# and only one entity too.
# They shouldn't be mandatory in these entries in the 1st place but try telling Eldon...
#
        tables = rc._get_tables( "natural_source" )
        if tables is None : raise Exception( "no natural source" )
        rc._make_nat_src( tables[0] )

        tables = rc._get_tables( "experimental_source" )
        if tables is None : raise Exception( "no experimental source" )
        rc._make_expt_src( tables[0] )

# boilerplate saveframes
#
        tables = rc._get_tables( "software" )
        if tables is None : raise Exception( "no software" )
        num = 1
        for table in tables :
            rc._make_software( table, softid = num )
            num += 1

#FIXME: there is only one now, may need to expand
#
        tables = rc._get_tables( "NMR_spectrometer" )
        if tables is None : raise Exception( "no spectroemter" )
        rc._make_spectrometer( tables[0] )

# assembly and entity are built from chem comp, we need to make that first
# there is only one chem comp in these entries
# FIXME if we ever get to multiple ones
#   e.g. comp ids are "BMETxxxxxx" but for multiple ones it'll have to be BMETxxxxxx_01, *_02, ...
#   image files need to be renamed (in entrydir.py), etc.
#
        tables = rc._get_tables( "chem_comp" )
        if (tables is None) or (len( tables ) < 1) : raise Exception( "no chem comps" )
        if len( tables ) > 1 :
            raise NotImplementedError( "don't know how to make %d chem comps" % (len( tables ),) )

        rc._make_chem_comp( tables[0], files = files )
        rc._make_entity()
        rc._make_assembly()

# must come before experiments: experiments have sample & conditions names but not ids.
#  want to create those ids first
# sample
#
        tables = rc._get_tables( "sample" )
        if tables is None : raise Exception( "no samples" )
        i = 1
        for table in tables :
            rc._make_sample( table, sampleid = i, assemblyid = 1, entityid = 1 )
            i += 1

# conditions
#
        tables = rc._get_tables( "sample_conditions" )
        if tables is None : raise Exception( "no sample conditions" )
        i = 1
        for table in tables :
            rc._make_sample_conditions( table, condsid = i )
            i += 1

# there can be only one
# experiment_files are for many datasets, should contain "dataset_id"
#  they need sample and conditions ids
#
        tables = rc._get_tables( "experiment_list" )
        if tables is None : raise Exception( "no experiment list" )
        rc._make_experiment_list( tables[0] )

# there can be separate CSrefs for different CSs
# like w/ experiments, CS have ref names but not ids
#
        tables = rc._get_tables( "chem_shift_reference" )
        if tables is None : raise Exception( "no CS ref" )
        i = 1
        for table in tables :
            rc._make_cs_reference( table, refid = i )
            i += 1

# since CS is linked 1-1 to sample conditions & CSref, we may have to make two
#
        tables = rc._get_tables( "assigned_chemical_shifts" )
        if tables is None : raise Exception( "no CS" )
        i = 1
        for table in tables :
            rc._make_chemical_shifts( table, listid = i )
            i += 1

# there's always peaks
#  be nice, sort them by experiment
#
        tables = rc._get_tables( "spectral_peak_list" )
        if tables is None : raise Exception( "no peaks" )
        num = 1
        for table in sorted( tables, cmp = lambda x, y : cmp( int( x["experiment_id"] ), int( y["experiment_id"] ) ) ) :
            rc._make_peak_list( table, localid = num )
            num += 1

        if not rc._is_sane() : raise Exception( "Failed sanity check" )

        return rc

    ################################################################################################
    #
    def __init__( self, stardict, entrydb, entrydata, entryid, verbose = False ) :
        self._dict = stardict
        self._star = entrydb
        self._json = entrydata
        self._id = entryid
        self._verbose = bool( verbose )
        self._title = None
        self._get_title()
        if self._title is None : raise Exception( "no molecule name" )

# FIXME! this is also in entrydir.py, should be refactored
#  (however here we always want unix path separator)
#
#        self._setid = setid
#        self._setname = "set%02d"  % (self._setid,)
#        self._datadir = "nmr/%s" % (self._setname,)
#        self._peakdir = "%s/transitions" % (self._datadir,)
#        self._specdir = "%s/spectra" % (self._datadir,)

    ################################################################################################
    #
    def to_star( self, out, err = sys.stderr ) :
        errs = []
        starobj.StarWriter.pretty_print( entry = self._star, dictionary = self._dict, out = out,
                errlist = errs, verbose = self._verbose )
        if len( errs ) > 0 :
            err.write( "------------------- unparse errors ---------------------\n" )
            for e in errs :
                err.write( str( e ) )
                err.write( "\n" )

    ################################################################################################
    # title is mol_common_name of sample component of type solute
    #
    def _get_title( self ) :
        for i in self._json :
            for j in i :
                for k in j.keys() :
                    if (k == "sample_component") and (isinstance( j[k], list )) :
                        for l in j[k] :
                            if l["type"] == "solute" :
                                self._title = l["mol_common_name"]
                                return

    ################################################################################################
    # incoming data is the list of saveframe categories,
    #   which is a list of saverames,
    #     which is a dict of tables (value is a list) and free tags (value is a primitive type)
    #
    def _get_tables( self, category ) :
        rc = []
        for i in self._json :
            for j in i :
                for k in j.keys() :
                    if k == "sf_category" :
                        if j[k] == category :
                            rc.append( j )
                            break
        if len( rc ) < 1 : return None
        return rc

    ################################################################################################
    # not all tables we need are in the input, so just create all and let the unparser sort them out
    #
    def _create_tables( self ) :
        self._star.create_tables( dictionary = self._dict, db = self._star._db, use_types = True,
                verbose = self._verbose )

    ################################################################################################
    #
    #
    def _make_entry_info( self, table, files ) :

        assert isinstance( table, dict )
        if not "authors" in table.keys() : raise Exception( "no entry authors" )
        if files is not None: assert isinstance( files, dict )

        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "entry_information" )

        today = datetime.date.today().isoformat()
        cols = { "Sf_category"                 : "'entry_information'",
                 "Sf_framecode"                : ":name",
                 "Sf_ID"                       : sfid,
                 "ID"                          : "'%s'" % (self._id,),
                 "Title"                       : ":title",
                 "Type"                        : ":type",
                 "Version_type"                : ":version_type",
                 "Submission_date"             : "'%s'" % (today,),
                 "Accession_date"              : "'%s'" % (today,),
                 "Last_release_date"           : "'%s'" % (today,),
                 "Original_release_date"       : "'%s'" % (today,),
                 "Origination"                 : ":origination",
                 "NMR_STAR_version"            : "'%s'" % (self._dict.version,),
                 "Original_NMR_STAR_version"   : "'%s'" % (self._dict.version,),
                 "Experimental_method"         : ":experimental_method",
                 "Experimental_method_subtype" : ":experimental_method_subtype",
                 "DOI"                         : "'%s'" % (self.DOI % (self._id.upper(),),) }

        params = { "name"                        : table["name"],
                   "title"                       : self._title,
                   "type"                        : table["type"],
                   "version_type"                : table["version_type"],
                   "origination"                 : table["origination"],
                   "experimental_method"         : table["experimental_method"],
                   "experimental_method_subtype" : table["experimental_method_subtype"] }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entry" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# entry authors
#FIXME: there's never ORCID or Family_title in this list. As of 20180911
#
        cols = { "Ordinal"         : ":num",
                 "Given_name"      : ":name",
                 "Family_name"     : ":sur",
                 "First_initial"   : ":ini",
                 "Middle_initials" : ":mini",
                 "Entry_ID"        : "'%s'" % (self._id,),
                 "Sf_ID"           : sfid }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entry_author" ("%s") values (%s)' % (colstr,valstr,)
        params.clear()
        for author in table["authors"] :
            params["num"] = author["ordinal"]
            params["name"] = author["given_name"]
            params["sur"] = author["family_name"]
            params["ini"] = author["first_initial"]
            if "middle_initials" in author.keys() :
                params["mini"] = author["middle_initials"]
            else :
                params["mini"] = None
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# experimental methods
#  we can live without them. probbly
#
        if "experimental_methods" in table.keys() :
            cols = { "ID"          : ":id",
                     "Method"      : ":met",
                     "Subtype"     : ":sub",
                     "Entry_ID"    : "'%s'" % (self._id,),
                     "Sf_ID"       : sfid }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Entry_experimental_methods" ("%s") values (%s)' % (colstr,valstr,)
            params.clear()
            for meth in table["experimental_methods"] :
                params["id"] = meth["id"]
                params["met"] = meth["type"]
                params["sub"] = meth["subtype"]
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# entry source -- there's only one
#
        if "source" in table.keys() :
            cols = { "ID"                     : ":id",
                     "Project_name"           : ":prj",
                     "Organization_full_name" : ":org",
                     "Organization_initials"  : ":oa",
                     "Entry_ID"               : "'%s'" % (self._id,),
                     "Sf_ID"                  : sfid }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Entry_src" ("%s") values (%s)' % (colstr,valstr,)
            params.clear()
            cnt = 1
            for src in table["source"] :
                params["id"] = cnt
                cnt += 1
                params["prj"] = src["project"]
                params["org"] = src["orgname_full"]
                params["oa"] = src["orgname_abbrev"]
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# entry files
#
        if files is not None :
            cols = { "ID"       : ":id",
                     "Name"     : ":uri",
                     "URI"      : ":uri",
                     "Format"   : ":type",
                     "Details"  : ":kind",
                     "Entry_ID" : "'%s'" % (self._id,),
                     "Sf_ID"    : sfid }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Auxiliary_files" ("%s") values (%s)' % (colstr,valstr,)
            params.clear()
            cnt = 1
            for f in files.keys() :
                params["id"] = cnt
                cnt += 1
                params["uri"] = files[f]
                if f == "SVG_N" :
                    params["type"] = "SVG"
                    params["kind"] = "molecule image with all atom labels"
                else :
                    params["type"] = f
                    if f in ("SVG","PNG") :
                        params["kind"] = "molecule image"
                    elif f in ("MOL","SDF") :
                        params["kind"] = "molecule structure file"
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# Data_set, Datum, and Release are inserted by release script
#

    ################################################################################################
    #
    #
    def _make_citation( self, table, citid ) :

        assert isinstance( table, dict )
        if not "authors" in table.keys() : raise Exception( "no citation authors (%s)" % (params["title"],) )

        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "citations" )

        cols = { "Sf_category"        : "'citations'",
                 "Sf_framecode"       : ":name",
                 "Sf_ID"              : ":sfid",
                 "ID"                 : ":id",
                 "Class"              : ":class",
                 "Type"               : ":type",
                 "Title"              : ":title",
                 "DOI"                : ":doi",
                 "Status"             : ":status",
                 "PubMed_ID"          : ":pmid",
                 "Journal_abbrev"     : ":jcode",
                 "Journal_name_full"  : ":jname",
                 "Journal_volume"     : ":jvol",
                 "Journal_issue"      : ":jisue",
                 "Journal_ISSN"       : ":jissn",
                 "Page_first"         : ":fpage",
                 "Page_last"          : ":lpage",
                 "Year"               : ":year",
                 "Entry_ID"           : "'%s'" % (self._id,) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Citation" ("%s") values (%s)' % (colstr,valstr,)

        params = { "name"  : table["name"],
                   "sfid"  : sfid,
                   "class" : table["class"],
                   "id"    : citid,
                   "type"  : table["type"] }

        if table["type"] == "BMRB only" :
            params["title"] = self._title
            params["doi"] = self.DOI % (self._id.upper(),)
            params["status"] = "published"
            params["pmid"] = None
            params["jcode"] = None
            params["jname"] = None
            params["jvol"] = None
            params["jisue"] = None
            params["jissn"] = None
            params["fpage"] = None
            params["lpage"] = None
            params["year"] = datetime.date.today().year
        elif table["type"] == "journal" :
            params["title"] = table["title"]
            params["doi"] = table["doi"]
            params["status"] = table["status"]
            params["pmid"] = table["pmid"]
            params["jcode"] = table["journal_abbrev"]
            params["jname"] = table["journal_name_full"]
            if "jouurnal_volume" in table.keys() :
                params["jvol"] = table["jouurnal_volume"]
            else :
                params["jvol"] = None
            params["jisue"] = table["journal_issue"]
            if "journal_issn" in table.keys() :
                params["jissn"] = table["journal_issn"]
            else :
                params["jissn"] = None
            params["fpage"] = table["page_first"]
            params["lpage"] = table["page_last"]
            params["year"] = table["year"]
        else :
            raise Exception( "don't know how to make citation of type %s" % (table["type"],) )

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# authors: like entry authors except Larry has family title
#
        cols = { "Ordinal"         : ":num",
                 "Given_name"      : ":name",
                 "Family_name"     : ":sur",
                 "First_initial"   : ":ini",
                 "Middle_initials" : ":mini",
                 "Family_title"    : ":moe",
                 "Entry_ID"        : "'%s'" % (self._id,),
                 "Citation_ID"     : citid,
                 "Sf_ID"           : sfid }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Citation_author" ("%s") values (%s)' % (colstr,valstr,)
        params.clear()
        for author in table["authors"] :
            params["num"] = author["ordinal"]
            params["name"] = author["given_name"]
            params["sur"] = author["family_name"]
            params["ini"] = author["first_initial"]
            if "middle_initials" in author.keys() :
                params["mini"] = author["middle_initials"]
            else :
                params["mini"] = None
            if "family_title" in author.keys() :
                params["moe"] = author["family_title"]
            else :
                params["moe"] = None
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #
    #
    def _make_nat_src( self, table, entityid = 1 ) :

        assert isinstance( table, dict )
        if not "natural_source" in table.keys() : raise Exception( "no natural_source table" )

        lclid = 1
        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "natural_source" )

        cols = { "Sf_category"  : "'natural_source'",
                 "Sf_framecode" : ":name",
                 "Sf_ID"        : sfid,
                 "ID"           : lclid,
                 "Entry_ID"     : "'%s'" % (self._id,) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity_natural_src_list" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, params = { "name" : table["name"] }, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# TODO: post-cook: insert entity name for entity by id. it's always 1/entity_1 but why not DIR.
#
        cols = { "Sf_ID"                      : sfid,
                 "ID"                         : ":id",
                 "Entity_ID"                  : entityid,
                 "Entry_ID"                   : "'%s'" % (self._id,),
                 "NCBI_taxonomy_ID"           : "'na'",
                 "Type"                       : ":type",
                 "Common"                     : "'yes'",
                 "Organism_name_scientific"   : "'na'",
                 "Organism_name_common"       : "'na'",
                 "Superkingdom"               : "'na'",
                 "Genus"                      : "'na'",
                 "Species"                    : "'na'",
                 "Entity_natural_src_list_ID" : lclid }

        params = {}
        cnt = 0
        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity_natural_src" ("%s") values (%s)' % (colstr,valstr,)

# there's only one
#
        for row in table["natural_source"] :
            cnt += 1
            params.clear()
            params["type"] = row["type"]
            params["id"] = cnt

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    # as above only with diffierent names
    #
    def _make_expt_src( self, table, entityid = 1 ) :

        assert isinstance( table, dict )
        if not "experimental_source" in table.keys() : raise Exception( "no experimental_source table" )

        lclid = 1
        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "experimental_source" )

        cols = { "Sf_category"  : "'experimental_source'",
                 "Sf_framecode" : ":name",
                 "Sf_ID"        : sfid,
                 "ID"           : lclid,
                 "Entry_ID"     : "'%s'" % (self._id,) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity_experimental_src_list" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, params = { "name" : table["name"] }, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# TODO: post-cook: insert entity name for entity by id. it's always 1/entity_1 but why not DIR.
#
        cols = { "Sf_ID"                           : sfid,
                 "ID"                              : ":id",
                 "Entity_ID"                       : entityid,
                 "Entry_ID"                        : "'%s'" % (self._id,),
                 "Production_method"               : ":meth",
                 "Entity_experimental_src_list_ID" : lclid }

        params = {}
        cnt = 0
        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity_experimental_src" ("%s") values (%s)' % (colstr,valstr,)

# there's only one
#
        for row in table["experimental_source"] :
            cnt += 1
            params.clear()
            params["meth"] = row["production_method"]
            params["id"] = cnt

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    # software
    #
    def _make_software( self, table, softid = 1 ) :
        assert isinstance( table, dict )
        if not "software_name" in table.keys() : raise Exception( "no software name" )

        if not "name" in table.keys() :
            sfname = "software_%d" % (softid,)
        else :
            sfname = table["name"]

        sfid = self._star.insert_saveframe( name = sfname, entryid = self._id, category = "software" )

        cols = { "Sf_category"  : "'software'",
                 "Sf_framecode" : ":sfname",
                 "Sf_ID"        : sfid,
                 "ID"           : softid,
                 "Entry_ID"     : "'%s'" % (self._id,),
                 "Name"         : ":name",
                 "Version"      : ":vers" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Software" ("%s") values (%s)' % (colstr,valstr,)
        params = { "sfname" : sfname, "name" : table["software_name"] }
        if "software_version" in table.keys() :
            params["vers"] = table["software_version"]
        else :
            params["vers"] = None
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# vendor
#  there's usually only one
#
        if "vendor" in table.keys() :
            cols = { "Sf_ID"              : sfid,
                     "Software_ID"        : softid,
                     "Entry_ID"           : "'%s'" % (self._id,),
                     "Name"               : ":name",
                     "Electronic_address" : ":url" }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Vendor" ("%s") values (%s)' % (colstr,valstr,)

            for row in table["vendor"] :
                params.clear()
                params["name"] = row["name"]
                if "url" in row.keys() :
                    params["url"] = row["url"]
                else :
                    params["url"] = None

                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# task
#  there's usually many
#
        if "task" in table.keys() :
            cols = { "Sf_ID"       : sfid,
                     "Software_ID" : softid,
                     "Entry_ID"    : "'%s'" % (self._id,),
                     "Task"        : ":name" }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Task" ("%s") values (%s)' % (colstr,valstr,)

            for row in table["task"] :
                params.clear()
                params["name"] = row["task"]

                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# citation
#  need sfid & label
#
        if "citation" in table.keys() :
            cols = { "Sf_ID"          : sfid,
                     "Software_ID"    : softid,
                     "Entry_ID"       : "'%s'" % (self._id,),
                     "Citation_ID"    : ":citid",
                     "Citation_label" : ":citlab" }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Software_citation" ("%s") values (%s)' % (colstr,valstr,)

            qry = 'select "ID" from "Citation" where "Sf_framecode"=:name'
            for row in table["citation"] :
                citid = None
                rs = self._star.query( qry, params = { "name" : row["name"] } )

# there can be only one
#
                for r in rs : citid = r[0]

                if citid is not None :
                    params.clear()
                    params["citid"] = citid
                    params["citlab"] = row["name"]
                    if self._verbose :
                        sys.stdout.write( sql )
                        sys.stdout.write( "\n" )
                        pprint.pprint( params )
                    rc = self._star.execute( sql, params = params, commit = True )
                    if self._verbose :
                        sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #
    def _make_spectrometer( self, table, localid = 1 ) :
        assert isinstance( table, dict )
        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "NMR_spectrometer" )
        cols = { "Sf_category"    : "'NMR_spectrometer'",
                 "Sf_framecode"   : "'%s'" % (table["name"],),
                 "Sf_ID"          : sfid,
                 "ID"             : localid,
                 "Entry_ID"       : "'%s'" % (self._id,),
                 "Model"          : ":mdl",
                 "Manufacturer"   : ":mfg",
                 "Field_strength" : ":fld",
                 "Details"        : ":dtl" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "NMR_spectrometer" ("%s") values (%s)' % (colstr,valstr,)
        params = { "mdl" : table["model"], "mfg" : table["manufacturer"], "fld" : table["field_strength"] }
        if "details" in table.keys() :
            params["dtl"] = table["details"]
        else :
            params["dtl"] = None
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #
    #
    def _make_chem_comp( self, table, files ) :

        assert isinstance( table, dict )
        if not "atoms" in table.keys() : raise Exception( "no atoms" )
        if not "bonds" in table.keys() : raise Exception( "no bonds" )

        assert isinstance( files, dict )

        m = re.search( "^[^\d]+(\d+)$", self._id )
        if not m : raise Exception( "id: bad pattern" )

        compid = "BMET%s" % (m.group( 1 ),)
        saveframe_name = re.sub( r"_+", "_", re.sub( r"[^A-Za-z0-9_]", "_", table["name"] ) )
        sfid = self._star.insert_saveframe( name = saveframe_name, entryid = self._id,
                category = "chem_comp" )

        molformat = "MOL"
        molfile = files[molformat]
        imgformat = "SVG"
        imgfile = files[imgformat]

        cols = { "Sf_category"  : "'chem_comp'",
                 "Sf_framecode" : ":sfname",
                 "Sf_ID"        : sfid,
                 "Entry_ID"     : "'%s'" % (self._id,),
                 "ID"           : "'%s'" % (compid,),
                 "Provenance"   : "'BMRB'",
                 "Name"         : ":molname",
                 "Type"         : "'non-polymer'",
                 "BMRB_code"    : "'%s'" % (compid,),
                 "Initial_date" : "'%s'" % (datetime.date.today().isoformat(),),
                 "Number_atoms_all" : table["num_atoms_all"],
                 "Number_atoms_nh"  : table["num_atoms_nh"],
                 "InChI_code"       : "'%s'" % (table["inchi_code"],),
                 "Formal_charge"    : table["charge"],
                 "Paramagnetic"     : "'%s'" % (table["paramagnetic"],),
                 "Formula"          : ":formula",
                 "Formula_weight"   : table["weight"],
                 "Formula_mono_iso_wt_nat" : table["weight_monoisotopic_nat"],
                 "Formula_mono_iso_wt_13C" : table["weight_monoisotopic_c13"],
                 "Formula_mono_iso_wt_15N" : table["weight_monoisotopic_n15"],
                 "Formula_mono_iso_wt_13C_15N" : table["weight_monoisotopic_c13n15"],
                 "Image_file_name"             : "'%s'" % (imgfile,),
                 "Image_file_format"           : "'%s'" % (imgformat,),
                 "Struct_file_name"            : "'%s'" % (molfile,),
                 "Struct_file_format"          : "'%s'" % (molformat,) }

        vendor = None
        if "vendor" in table.keys() :
            vendor = table["vendor"]
            cols["Vendor"] = "'%s'" % (vendor,)
        vendor_code = None
        if "vendor_product_code" in table.keys() :
            vendor_code = table["vendor_product_code"]
            cols["Vendor_product_code"] = "'%s'" % (vendor_code,)

# TODO:
#  set _Chem_comp.Aromatic based on bonds & atoms

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Chem_comp" ("%s") values (%s)' % (colstr,valstr,)

        params = { "sfname" : saveframe_name, "molname" : self._title, "formula" : table["formula"] }
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# synonyms
#
        if "common_name" in table.keys() :
            cols = { "Sf_ID"    : sfid,
                     "Entry_ID" : "'%s'" % (self._id,),
                     "Comp_ID"  : "'%s'" % (compid,),
                     "Name"     : ":name",
                     "Type"     : ":type" }
            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Chem_comp_common_name" ("%s") values (%s)' % (colstr,valstr,)
            params = {}
            for row in table["common_name"] :
                params.clear()

#FIXME: assume it's either a dict w/ { "name" : x, "type" : y } or a list of names
#
                if isinstance( row, dict ) :
                    params["name"] = row["name"]
                    params["type"] = row["type"]
                else :
                    params["name"] = row
                    params["type"] = "name"
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# circa 2017 entries have alatis inchi string in the index but not in sdf
# from mid-2018 on alatis adds trhe <alatis_inchi{ tag to sdf,
#  that's pulled out into descriptrs table upstream
#
        if "descriptors" in table.keys() :
            found = False
            for row in table["descriptors"] :
                if row["program"] == "ALATIS" :
                    found = True
                    break
            if not found :
                table["descriptors"].append( { "type" : "InChI", "program" : "ALATIS",
                    "descriptor" : table["inchi_code"] } )
        else :
            table["descriptors"] = [ { "type" : "InChI", "program" : "ALATIS",
                    "descriptor" : table["inchi_code"] } ]

# descriptors -- these should always be there, actially, w/ at least the inchi string
#
        if "descriptors" in table.keys() :
            cols = { "Sf_ID"           : sfid,
                     "Entry_ID"        : "'%s'" % (self._id,),
                     "Comp_ID"         : "'%s'" % (compid,),
                     "Descriptor"      : ":name",
                     "Type"            : ":type",
                     "Program"         : ":prog",
                     "Program_version" : ":vers" }
            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Chem_comp_descriptor" ("%s") values (%s)' % (colstr,valstr,)
            params = {}
            for row in table["descriptors"] :
                params.clear()
                params["name"] = row["descriptor"]
                params["type"] = row["type"]
                params["prog"] = row["program"]
                if "version" in row.keys() :
                    params["vers"] = row["version"]
                else :
                    params["vers"] = "na"
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# identifiers: exactly like above
# except for the exceptions
# special case: NCI compunds have "NCI" for vendor and "NSCxxxxxx" for catalog number
# circa 2017 entries have alatis inchi string in the index but not in sdf
#
        if (vendor_code is not None) and (str( vendor ).upper() == "NCI") :
            if str( vendor_code ).strip().upper().startswith( "NSC" ) :
                nsc_num = str( vendor_code ).strip()[3:]
                if "identifiers" in table.keys() :
                    found = False
                    for row in table["identifiers"] :
                        if (row["type"] == "NSC NUMBER") and (row["identifier"] == nsc_num) :
                            found = True
                            break
                    if not found :
                        table["identifiers"].append( { "type" : "NSC NUMBER",
                                "identifier" : nsc_num, "program" : "na", "version" : "na" } )
                else :
                    table["identifiers"] = [ { "type" : "NSC NUMBER",
                                "identifier" : nsc_num, "program" : "na", "version" : "na" } ]
            else :
                raise Exception( "Vendor is NCI but can't get NSC code" )

        if "identifiers" in table.keys() :

            cols = { "Sf_ID"           : sfid,
                     "Entry_ID"        : "'%s'" % (self._id,),
                     "Comp_ID"         : "'%s'" % (compid,),
                     "Identifier"      : ":name",
                     "Type"            : ":type",
                     "Program"         : ":prog",
                     "Program_version" : ":vers" }
            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Chem_comp_identifier" ("%s") values (%s)' % (colstr,valstr,)
            params = {}
            for row in table["identifiers"] :
                params.clear()
                params["name"] = row["identifier"]
                params["type"] = row["type"]
                params["prog"] = row["program"]
                if "version" in row.keys() :
                    params["vers"] = row["version"]
                else :
                    params["vers"] = "na"
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )


# atoms
#
        cols = { "Sf_ID"         : sfid,
                 "Entry_ID"      : "'%s'" % (self._id,),
                 "Comp_ID"       : "'%s'" % (compid,),
                 "Atom_ID"       : ":name",
                 "Type_symbol"   : ":type",
                 "Stereo_config" : ":stereo",
                 "Charge"        : ":charge",
                 "Aromatic_flag" : ":arom",
                 "PDBX_ordinal"  : ":num" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Chem_comp_atom" ("%s") values (%s)' % (colstr,valstr,)
        aromatic = False
        params = {}
        for row in table["atoms"] :
            params.clear()
            params["name"] = row["id"]
            params["type"] = row["type"]
            params["stereo"] = row["stereo"]
            params["charge"] = row["charge"]
            if row["aromatic"] :
                params["arom"] = "yes"
                aromatic = True
            else :
                params["arom"] = "no"
            params["num"] = row["ordinal"]
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# bonds
#
        cols = { "Sf_ID"         : sfid,
                 "Entry_ID"      : "'%s'" % (self._id,),
                 "Comp_ID"       : "'%s'" % (compid,),
                 "Type"          : ":type",
                 "Value_order"   : ":ord",
                 "Atom_ID_1"     : ":atm1",
                 "Atom_ID_2"     : ":atm2",
                 "Aromatic_flag" : ":arom",
                 "Stereo_config" : ":stereo",
                 "Ordinal"       : ":num",
                 "ID"            : ":num" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Chem_comp_bond" ("%s") values (%s)' % (colstr,valstr,)

        params = {}
        for row in table["bonds"] :
            params.clear()
            params["num"] = row["id"]
            params["type"] = row["type"]
            params["stereo"] = row["stereo"]
            params["ord"] = row["ord"]
            params["atm1"] = row["atom1"]
            params["atm2"] = row["atom2"]
            params["arom"] = row["arom"]
            if row["arom"] == "yes" :
                aromatic = True

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# aromatic fix
#
        if aromatic :
            sql = """update "Chem_comp" set "Aromatic"='yes'"""
        else :
            sql = """update "Chem_comp" set "Aromatic"='no'"""
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# atom nomenclature -- only exists in 8 pre-alatis entries ATM.
#
        if "atom_nomenclature" in table.keys() :
            cols = { "Sf_ID"         : sfid,
                     "Entry_ID"      : "'%s'" % (self._id,),
                     "Comp_ID"       : "'%s'" % (compid,),
                     "Atom_ID"       : ":atm",
                     "Atom_name"     : ":alt",
                     "Naming_system" : ":sys" }
            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Atom_nomenclature" ("%s") values (%s)' % (colstr,valstr,)
            params = {}
            for row in table["atom_nomenclature"] :
                params.clear()
                params["atm"] = row["atom_id"]
                params["alt"] = row["atom_name"]
                params["sys"] = row["naming_system"]
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# DB links
#
        if "dblinks" in table.keys() :
            cols = { "Sf_ID"               : sfid,
                     "Entry_ID"            : "'%s'" % (self._id,),
                     "Comp_ID"             : "'%s'" % (compid,),
                     "Author_supplied"     : "'no'",
                     "Database_code"       : ":db",
                     "Accession_code"      : ":acc",
                     "Accession_code_type" : ":actype",
                     "Entry_relation_type" : "'same molecule'" }
            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Chem_comp_db_link" ("%s") values (%s)' % (colstr,valstr,)
            params = {}
            for row in table["dblinks"] :
                params.clear()
                params["db"] = row["db_code"]
                params["acc"] = row["acc_code"]
                params["actype"] = row["acc_type"]
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    # entity is copied from chem comp
    # this has to run after chem comp's in the database
    #
    def _make_entity( self, entityid = 1 ) :

        sfname = "entity_%d" % (entityid,)

        sfid = self._star.insert_saveframe( name = sfname, entryid = self._id, category = "entity" )

        params = {}

# there's only one chem comp
#
        sql = 'select "Name","ID","Sf_framecode","Paramagnetic","Formula_weight" from "Chem_comp"'
        rc = self._star.query( sql )
        for row in rc :
            params["compname"] = row[0]
            params["compid"] = row[1]
            params["complabel"] = row[2]
            params["paramag"] = row[3]
            params["mass"] = row[4]

        cols = { "Sf_category"                     : "'entity'",
                 "Sf_framecode"                    : "'%s'" % (sfname,),
                 "Sf_ID"                           : sfid,
                 "ID"                              : entityid,
                 "Entry_ID"                        : "'%s'" % (self._id,),
                 "Type"                            : "'non-polymer'",
                 "Ambiguous_conformational_states" : "'no'",
                 "Ambiguous_chem_comp_sites"       : "'no'",
                 "Nstd_linkage"                    : "'no'",
                 "Number_of_monomers"              : 1,
                 "Number_of_nonpolymer_components" : 1,
                 "Thiol_state"                     : "'not present'",
                 "Name"                            : ":compname",
                 "Nonpolymer_comp_ID"              : ":compid",
                 "Nonpolymer_comp_label"           : ":complabel",
                 "Paramagnetic"                    : ":paramag",
                 "Formula_weight"                  : ":mass" }

# ???
#
# _Entity.Nstd_monomer                     no
# _Entity.Nstd_chirality                   no


        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

#
#
        cols = { "Sf_ID"       : sfid,
                 "ID"          : 1,
                 "Auth_seq_ID" : 1,
                 "Entity_ID"   : entityid,
                 "Entry_ID"    : "'%s'" % (self._id,),
                 "Comp_ID"     : ":compid",
                 "Comp_label"  : ":complabel" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity_comp_index" ("%s") values (%s)' % (colstr,valstr,)

# there's only one
#

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    # assembly is copied from chem comp & entity
    # this has to run after entity is in the database
    # in theory nmr-star allows for more than one assembly...
    #
    def _make_assembly( self, assemblyid = 1 ) :

        sfname = "assembly_%d" % (assemblyid,)

        sfid = self._star.insert_saveframe( name = sfname, entryid = self._id, category = "assembly" )

        params = {}

# there's only one entity
#
        sql = 'select "Name","ID","Sf_framecode","Paramagnetic","Thiol_state" from "Entity"'
        rs = self._star.query( sql )
        for row in rs :
            params["compname"] = row[0]
            params["eid"] = row[1]
            params["elabel"] = row[2]
            params["paramag"] = row[3]
            params["thiol"] = row[4]

        cols = { "Sf_category"           : "'assembly'",
                 "Sf_framecode"          : "'%s'" % (sfname,),
                 "Sf_ID"                 : sfid,
                 "ID"                    : assemblyid,
                 "Entry_ID"              : "'%s'" % (self._id,),
                 "Number_of_components"  : 1,
                 "Thiol_state"           : ":thiol",
                 "Name"                  : ":compname",
                 "Paramagnetic"          : ":paramag" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Assembly" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

#
#
        cols = { "Sf_ID"                      : sfid,
                 "ID"                         : 1,
                 "Entry_ID"                   : "'%s'" % (self._id,),
                 "Assembly_ID"                : assemblyid,
                 "Experimental_data_reported" : "'yes'",
                 "Physical_state"             : "'native'",
                 "Conformational_isomer"      : "'no'",
                 "Chemical_exchange_state"    : "'no'",
                 "Entity_assembly_name"       : ":compname",
                 "Entity_ID"                  : ":eid",
                 "Entity_label"               : ":elabel" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Entity_assembly" ("%s") values (%s)' % (colstr,valstr,)

# there's only one
#

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
            pprint.pprint( params )
        rc = self._star.execute( sql, params = params, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# atoms are copied from chem comp
#  comp_index_id = seq_num is always 1, hard-coded instead of doing a join on comp_index table
#
        entityid = params["eid"]
        cols = { "Sf_ID"              : sfid,
                 "Entry_ID"           : "'%s'" % (self._id,),
                 "Assembly_ID"        : assemblyid,
                 "Entity_assembly_ID" : 1,
                 "Comp_index_ID"      : 1,
                 "Seq_ID"             : 1,
                 "Entity_ID"          : ":eid",
                 "Comp_ID"            : ":lab",
                 "Atom_ID"            : ":atm",
                 "Type_symbol"        : ":nuc",
                 "Assembly_atom_ID"   : ":num" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Atom" ("%s") values (%s)' % (colstr,valstr,)

# wtf. 1st 2 rows are returned twice. even with distinct
#
        atoms = set()
        qry = 'select "Comp_ID","Atom_ID","Type_symbol","PDBX_ordinal" from "Chem_comp_atom" order by "PDBX_ordinal"'
        rs = self._star.query( qry, newcursor = True )
        for row in rs :
            if row[1] in atoms : continue
            else : atoms.add( row[1] )
            params.clear()
            params["eid"] = entityid
            params["lab"] = row[0]
            params["atm"] = row[1]
            params["nuc"] = row[2]
            params["num"] = row[3]

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #
    def _make_experiment_list( self, table ) :
        assert isinstance( table, dict )
        if not "experiment" in table.keys() : raise Exception( "no experiments" )
        if not "experiment_file" in table.keys() : raise Exception( "no experiment files" )
        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "experiment_list" )
        localid = 1
        cols = { "Sf_category"    : "'experiment_list'",
                 "Sf_framecode"   : "'%s'" % (table["name"],),
                 "Sf_ID"          : sfid,
                 "ID"             : localid,
                 "Entry_ID"       : "'%s'" % (self._id,) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Experiment_list" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# experiments
#
        cols = { "Sf_ID"                       : sfid,
                 "Experiment_list_ID"          : localid,
                 "Entry_ID"                    : "'%s'" % (self._id,),
                 "Sample_state"                : "'isotropic'",
                 "NMR_spectrometer_label"      : ":spect_lbl",
                 "Sample_label"                : ":smpl_lbl",
                 "Sample_condition_list_label" : ":cond_lbl",
                 "ID"                          : ":id",
                 "Name"                        : ":name",
                 "Raw_data_flag"               : ":flg" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Experiment" ("%s") values (%s)' % (colstr,valstr,)
        params = {}
        for row in table["experiment"] :
            params.clear()
            params["name"] = row["name"]
            params["id"] = row["id"]
            params["flg"] = row["raw_data_flag"]
            params["spect_lbl"] = row["nmr_spectrometer_label"]
            params["smpl_lbl"] = row["sample_label"]
            params["cond_lbl"] = row["sample_condition_list_label"]

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        sql = 'update "Experiment" set "Sample_ID"=' \
            + '(select "ID" from "Sample" where "Sf_framecode"="Experiment"."Sample_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Experiment" set "Sample_condition_list_ID"=(select "ID" from ' \
            + '"Sample_condition_list" where "Sf_framecode"="Experiment"."Sample_condition_list_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Experiment" set "NMR_spectrometer_ID"=(select "ID" from ' \
            + '"NMR_spectrometer" where "Sf_framecode"="Experiment"."NMR_spectrometer_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

# files
#
        cols = { "Sf_ID"              : sfid,
                 "Experiment_list_ID" : localid,
                 "Entry_ID"           : "'%s'" % (self._id,),
                 "Experiment_ID"      : ":exptid",
                 "Name"               : ":name",
                 "Type"               : ":type",
                 "Content"            : ":cont",
                 "Directory_path"     : ":path",
                 "Details"            : ":dtl" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Experiment_file" ("%s") values (%s)' % (colstr,valstr,)
        params = {}
        pat = re.compile( r"\s+" )
        for row in table["experiment_file"] :
            params.clear()
            params["exptid"] = row["experiment_id"]
            params["type"] = row["type"]
            params["cont"] = row["content"]

            if "details" in row.keys() :
                params["dtl"] = row["details"]
            else :
                params["dtl"] = None

# text/directory is time-domain data, text/xml is peak list, image/png is spectrum pix
#  could go by content instead
# path is relative to entry directory
# FIXME:
#  peak lists and images are renamed to experiment names (with s/\s/_/) by entrydirs.py
#
            dataset_id = 1
            if "dataset_id" in row.keys() :
                dataset_id = int( row["dataset_id"] )

            if row["type"].lower() == "text/directory" :
                params["path"] = DATASET_DIR % (dataset_id,)
                params["name"] = row["name"]
            else :
                found = False
                for e in table["experiment"] :
                    if e["id"] == row["experiment_id"] :
                        found = True
                        expname = pat.sub( "_", e["name"] )
                        if row["type"].lower() == "text/xml" :
                            dstdir = PEAKFILE_DIR  % (dataset_id,)
                            params["path"] = "%s/%s" % (dstdir,expname,)
                            params["name"] = "%s.xml" % (expname,)
                        elif row["type"].lower() == "image/png" :
                            dstdir = SPECTRA_DIR  % (dataset_id,)
                            params["path"] = "%s/%s" % (dstdir,expname,)
                            params["name"] = "%s.png" % (expname,)
                        else :
                            raise Exception( "Can't figure path for %s" % (row["type"],) )
                        break
                if not found :
                            raise Exception( "Can't fiind experiment for %s" % (row["name"],) )

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #
    def _make_sample( self, table, sampleid = 1, assemblyid = 1, entityid = 1 ) :

        assert isinstance( table, dict )
        if not "sample_component" in table.keys() : raise Exception( "no sample components" )

        try :
            lclid = int( sampleid )
        except :
            lclid = 1

#FIXME: change to name in the other place
#
        if "name" in table.keys() : name = table["name"]
        else :  name = table["sf_framecode"]

        sfid = self._star.insert_saveframe( name = name, entryid = self._id,
                category = "sample" )

        cols = { "Sf_category"  : "'sample'",
                 "Sf_framecode" : ":name",
                 "Sf_ID"        : sfid,
                 "Entry_ID"     : "'%s'" % (self._id,),
                 "ID"           : lclid,
                 "Type"         : "'%s'" % (table["type"],) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Sample" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, params = { "name" : name }, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# components
# POST-COOK: resolve entity and assembly labels from IDs
#
        cols = { "ID"                      : ":id",
                 "Mol_common_name"         : ":molname",
                 "Isotopic_labeling"       : ":isolab",
                 "Assembly_ID"             : assemblyid,
                 "Entity_ID"               : entityid,
                 "Type"                    : ":type",
                 "Concentration_val"       : ":conc",
                 "Concentration_val_units" : ":cuni",
                 "Vendor"                  : ":vend",
                 "Vendor_product_name"     : ":vname",
                 "Vendor_product_code"     : ":vcode",
                 "Entry_ID"                : "'%s'" % (self._id,),
                 "Sample_ID"               : lclid,
                 "Sf_ID"                   : sfid }

        params = {}
        cnt = 0
        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Sample_component" ("%s") values (%s)' % (colstr,valstr,)
        for row in table["sample_component"] :
            cnt += 1
            params.clear()
            params["id"] = row["id"]
            params["molname"] = row["mol_common_name"]
            params["conc"] = row["concentration_val"]
            params["cuni"] = row["concentration_val_units"]
            params["type"] = row["type"]
            if row["type"] == "solute" :
                vendor = row["vendor"]
                params["vend"] = vendor
                vcode = row["vendor_product_code"]
                params["vcode"] = vcode
                params["vname"] = row["vendor_product_name"]
                params["isolab"] = row["isotopic_labeling"]
            else :
                params["vend"] = None
                params["vcode"] = None
                params["vname"] = None
                params["isolab"] = None

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        return (vendor, vcode)

    ################################################################################################
    #
    def _make_sample_conditions( self, table, condsid = 1 ) :

        assert isinstance( table, dict )
        if not "sample_condition_variable" in table.keys() : raise Exception( "no sample cocondition vars" )

        try :
            lclid = int( condsid )
        except :
            lclid = 1

#FIXME: change to name in the other place
#
        if "name" in table.keys() : name = table["name"]
        else :  name = table["sf_framecode"]

        sfid = self._star.insert_saveframe( name = name, entryid = self._id,
                category = "sample_conditions" )

        cols = { "Sf_category"  : "'sample_conditions'",
                 "Sf_framecode" : ":name",
                 "Sf_ID"        : sfid,
                 "Entry_ID"     : "'%s'" % (self._id,),
                 "ID"           : lclid }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Sample_condition_list" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, params = { "name" : name }, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# data
#
        cols = { "Type"                     : ":type",
                 "Val"                      : ":value",
                 "Val_err"                  : ":error",
                 "Val_units"                : ":units",
                 "Entry_ID"                 : "'%s'" % (self._id,),
                 "Sample_condition_list_ID" : lclid,
                 "Sf_ID"                    : sfid }

        params = {}
        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Sample_condition_variable" ("%s") values (%s)' % (colstr,valstr,)
        for row in table["sample_condition_variable"] :
            params.clear()
            params["type"] = row["type"]
            params["value"] = row["val"]
            if "val_err" in row.keys() :
                params["error"] = row["val_err"]
            else :
                params["error"] = None
            params["units"] = row["val_units"]
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )


    ################################################################################################
    #
    def _make_cs_reference( self, table, refid = 1 ) :
        assert isinstance( table, dict )

        try :
            localid = int( refid )
        except :
            localid = 1

        if not "chem_shift_ref" in table.keys() : raise Exception( "no CS reference loop" )
        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "chem_shift_reference" )
        cols = { "Sf_category"  : "'chem_shift_reference'",
                 "Sf_framecode" : "'%s'" % (table["name"],),
                 "Sf_ID"        : sfid,
                 "ID"           : localid,
                 "Entry_ID"     : "'%s'" % (self._id,) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Chem_shift_reference" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        cols = { "Sf_ID"                   : sfid,
                 "Chem_shift_reference_ID" : localid,
                 "Entry_ID"                : "'%s'" % (self._id,),
                 "Atom_type"               : ":nuc",
                 "Atom_isotope_number"     : ":iso",
                 "Mol_common_name"         : ":ref",
                 "Atom_group"              : ":grp",
                 "Chem_shift_units"        : ":unit",
                 "Chem_shift_val"          : ":val",
                 "Ref_method"              : ":meth",
                 "Ref_type"                : ":type",
                 "Indirect_shift_ratio"    : ":rat" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Chem_shift_ref" ("%s") values (%s)' % (colstr,valstr,)
        params = {}

        for row in table["chem_shift_ref"] :
            params["nuc"] = row["element"]
            params["iso"] = row["isotope"]
            params["ref"] = row["molecule"]
            params["grp"] = row["group"]
            params["unit"] = row["units"]
            params["val"] = row["val"]
            params["meth"] = row["method"]
            params["type"] = row["type"]
            params["rat"] = row["ratio"]

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #
    def _make_chemical_shifts( self, table, listid = 1, assemblyid = 1, entityid = 1, seqid = 1 ) :

        assert isinstance( table, dict )
        try :
            localid = int( listid )
        except :
            localid = 1

        if not "atom_chem_shift" in table.keys() : raise Exception( "no CS loop" )
        if not "chem_shift_experiment" in table.keys() : raise Exception( "no CS expt loop" )
        if not "chem_shift_reference_label" in table.keys() : raise Exception( "no CS ref label" )
        if not "sample_condition_list_label" in table.keys() : raise Exception( "no conds label" )

        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "assigned_chemical_shifts" )
        cols = { "Sf_category"              : "'assigned_chemical_shifts'",
                 "Sf_framecode"             : "'%s'" % (table["name"],),
                 "Sf_ID"                    : sfid,
                 "ID"                       : localid,
                 "Entry_ID"                 : "'%s'" % (self._id,),
                 "Chem_shift_reference_label"  : "'%s'" % (table["chem_shift_reference_label"],),
                 "Sample_condition_list_label" : "'%s'" % (table["sample_condition_list_label"],) }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Assigned_chem_shift_list" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        sql = 'update "Assigned_chem_shift_list" set "Sample_condition_list_ID"=' \
            + '(select "ID" from "Sample_condition_list" where "Sf_framecode"=' \
            + '"Assigned_chem_shift_list"."Sample_condition_list_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Assigned_chem_shift_list" set "Chem_shift_reference_ID"=' \
            + '(select "ID" from "Chem_shift_reference" where "Sf_framecode"=' \
            + '"Assigned_chem_shift_list"."Chem_shift_reference_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

# expts
#
        cols = { "Sf_ID"                       : sfid,
                 "Assigned_chem_shift_list_ID" : localid,
                 "Entry_ID"                    : "'%s'" % (self._id,),
                 "Sample_label"                : ":sampleid",
                 "Experiment_ID"               : ":expid",
                 "Experiment_name"             : ":expname" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Chem_shift_experiment" ("%s") values (%s)' % (colstr,valstr,)
        params = {}

        for row in table["chem_shift_experiment"] :
            params.clear()
            params["expid"] = row["experiment_id"]
            params["expname"] = row["experiment_name"]
            if "sample_label" in row.keys() :
                params["sampleid"] = row["sample_label"]

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        sql = 'update "Chem_shift_experiment" set "Sample_ID"=' \
            + '(select "ID" from "Sample" where "Sf_framecode"="Chem_shift_experiment"."Sample_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

# shifts
#  assume the input is sorted already.
#  ids are wrong but they match ids in ambiguity table, need to renumber and map
#
        cols = { "Sf_ID"                       : sfid,
                 "Assigned_chem_shift_list_ID" : localid,
                 "Entry_ID"                    : "'%s'" % (self._id,),
                 "Entity_assembly_ID"          : assemblyid,
                 "Entity_ID"                   : entityid,
                 "Comp_index_ID"               : seqid,
                 "Seq_ID"                      : seqid,
                 "ID"                          : ":num",
                 "Resonance_ID"                : ":rid",
                 "Atom_ID"                     : ":atm",
                 "Val"                         : ":val",
                 "Ambiguity_code"              : ":amb" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Atom_chem_shift" ("%s") values (%s)' % (colstr,valstr,)
        params = {}
        num = 1
        shift_numbers = {}
        for row in table["atom_chem_shift"] :
            params.clear()
            params["num"] = num
            shift_numbers[int( row["id"] )] = num
            params["atm"] = row["atom_id"]
            params["val"] = row["value"]
            params["amb"] = row["ambiguity_code"]
            if "resonance_id" in row.keys() :
                params["rid"] = row["resonance_id"]
            else :
                params["rid"] = None

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

            num += 1

        sql = 'update "Atom_chem_shift" set "Comp_ID"=(select "ID" from "Chem_comp")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

#        sql = 'update "Atom_chem_shift" set "Assembly_atom_ID"=' \
#            + '(select "Assembly_atom_ID" from "Atom" where ' \
#            + '"Entity_assembly_ID"="Atom_chem_shift"."Entity_assembly_ID" and ' \
#            + '"Entity_ID"="Atom_chem_shift"."Entity_ID" and ' \
#            + '"Comp_index_ID"="Atom_chem_shift"."Comp_index_ID" and ' \
#            + '"Atom_ID"="Atom_chem_shift"."Atom_ID")'

        sql = 'update "Atom_chem_shift" set "Assembly_atom_ID"=' \
            + '(select "Assembly_atom_ID" from "Atom" where "Entity_assembly_ID"=%s ' % (assemblyid,)
        sql += 'and "Entity_ID"=%s and "Comp_index_ID"=%s ' % (entityid, seqid)
        sql += 'and "Atom_ID"="Atom_chem_shift"."Atom_ID")'

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Atom_chem_shift" set "Atom_type"=' \
            + '(select "Type_symbol" from "Chem_comp_atom" where ' \
            + '"Comp_ID"="Atom_chem_shift"."Comp_ID" and ' \
            + '"Atom_ID"="Atom_chem_shift"."Atom_ID")'

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = """update "Atom_chem_shift" set "Atom_isotope_number"=1 where "Atom_type"='H'"""

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = """update "Atom_chem_shift" set "Atom_isotope_number"=13 where "Atom_type"='C'"""

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

# ambiguities
#
        if "ambiguous_atom_chem_shift" in table.keys() :
            cols = { "Sf_ID"                       : sfid,
                     "Assigned_chem_shift_list_ID" : localid,
                     "Entry_ID"                    : "'%s'" % (self._id,),
                     "Ambiguous_shift_set_ID"      : ":set",
                     "Atom_chem_shift_ID"          : ":shift" }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Ambiguous_atom_chem_shift" ("%s") values (%s)' % (colstr,valstr,)
            params = {}

            for row in table["ambiguous_atom_chem_shift"] :
                params.clear()
                params["set"] = row["set_id"]
                params["shift"] = shift_numbers[int( row["shift_id"] )]

                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

            sql = 'update "Atom_chem_shift" set "Ambiguity_set_ID"=' \
                + '(select "Ambiguous_shift_set_ID" from "Ambiguous_atom_chem_shift" ' \
                + 'where "Ambiguous_atom_chem_shift"."Entry_ID"="Atom_chem_shift"."Entry_ID" ' \
                + 'and "Ambiguous_atom_chem_shift"."Assigned_chem_shift_list_ID"="Atom_chem_shift"."Assigned_chem_shift_list_ID" ' \
                + 'and "Ambiguous_atom_chem_shift"."Atom_chem_shift_ID"="Atom_chem_shift"."ID")'

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
            rc = self._star.execute( sql, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

# all processing is done in topspin anyway, but be a little clever
#  (but not look for "chemical shift assignement" because that's not in the task)
#
        sql = 'insert into "Chem_shift_software" ("Sf_ID","Entry_ID","Assigned_chem_shift_list_ID",' \
            + '"Software_ID","Software_label") ' \
            + """select %d,'%s',%d,s."ID",s."Sf_framecode" """ % (sfid,self._id,localid,) \
            + """from "Software" s,"Task" t where t."Sf_ID"=s."Sf_ID" and t."Task"='data analysis'"""

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    ################################################################################################
    #

    def _make_peak_list( self, table, localid = 1 ) :

        assert isinstance( table, dict )
        assert "assigned_chem_shift_list_label" in table.keys()
        assert "sample_condition_list_label" in table.keys()
        assert "sample_label" in table.keys()

        if not "spectral_dim" in table.keys() : raise Exception( "spectral dim" )
        if not "spectral_transition" in table.keys() : raise Exception( "no transition list" )
        if not "spectral_transition_char" in table.keys() : raise Exception( "no transition shifts" )
        if not "spectral_transition_general_char" in table.keys() : raise Exception( "no transition intensities" )

        sfid = self._star.insert_saveframe( name = table["name"], entryid = self._id,
                category = "spectral_peak_list" )
        cols = { "Sf_category"                   : "'spectral_peak_list'",
                 "Sf_framecode"                  : "'%s'" % (table["name"],),
                 "Sf_ID"                         : sfid,
                 "ID"                            : localid,
                 "Entry_ID"                      : "'%s'" % (self._id,),
                 "Sample_label"                     : "'%s'" % (table["sample_label"],),
                 "Sample_condition_list_label"      : "'%s'" % (table["sample_condition_list_label"],),
                 "Assigned_chem_shift_list_label"   : "'%s'" % (table["assigned_chem_shift_list_label"],),
                 "Number_of_spectral_dimensions" : "'%s'" % (table["number_of_spectral_dimensions"],),
                 "Experiment_ID"                 : "'%s'" % (table["experiment_id"],),
                 "Experiment_name"               : "'%s'" % (table["experiment_name"],),
                 "Details"                       : ":dtl" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Spectral_peak_list" ("%s") values (%s)' % (colstr,valstr,)
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, params = { "dtl" : table["details"] }, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        sql = 'update "Spectral_peak_list" set "Sample_ID"=' \
            + '(select "ID" from "Sample" where "Sf_framecode"="Spectral_peak_list"."Sample_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Spectral_peak_list" set "Sample_condition_list_ID"=' \
            + '(select "ID" from "Sample_condition_list" where "Sf_framecode"=' \
            + '"Spectral_peak_list"."Sample_condition_list_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Spectral_peak_list" set "Assigned_chem_shift_list_ID"=' \
            + '(select "ID" from "Assigned_chem_shift_list" where "Sf_framecode"=' \
            + '"Spectral_peak_list"."Assigned_chem_shift_list_label")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Spectral_peak_list" set "Chem_shift_reference_ID"=' \
            + '(select "Chem_shift_reference_ID" from "Assigned_chem_shift_list"' \
            +  'where "ID"="Spectral_peak_list"."Assigned_chem_shift_list_ID")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

        sql = 'update "Spectral_peak_list" set "Chem_shift_reference_label"=' \
            + '(select "Chem_shift_reference_label" from "Assigned_chem_shift_list"' \
            +  'where "ID"="Spectral_peak_list"."Assigned_chem_shift_list_ID")'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

# as with chem shifts except here there is exact 'peak picking' task
#
        sql = 'insert into "Spectral_peak_software" ("Sf_ID","Entry_ID","Spectral_peak_list_ID",' \
            + '"Software_ID","Software_label") ' \
            + """select %d,'%s',%d,s."ID",s."Sf_framecode" """ % (sfid,self._id,localid,) \
            + """from "Software" s,"Task" t where t."Sf_ID"=s."Sf_ID" and t."Task"='peak picking'"""

        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rc = self._star.execute( sql, commit = True )
        if self._verbose :
            sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# dimensions
#
        cols = { "Sf_ID"                 : sfid,
                 "Spectral_peak_list_ID" : localid,
                 "Entry_ID"              : "'%s'" % (self._id,),
                 "ID"                    : ":num",
                 "Atom_type"             : ":atm",
                 "Atom_isotope_number"   : ":iso",
                 "Spectral_region"       : ":reg",
                 "Sweep_width"           : ":sw",
                 "Sweep_width_units"     : ":units" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Spectral_dim" ("%s") values (%s)' % (colstr,valstr,)
        params = {}
        num = 1
        shift_numbers = {}
        for row in table["spectral_dim"] :
            params.clear()
            params["num"] = row["id"]
            params["atm"] = row["atom_type"]
            params["iso"] = row["atom_isotope_number"]
            params["reg"] = row["spectral_region"]
            params["sw"] = row["sweep_width"]
            params["units"] = row["sweep_width_units"]

            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# and transitions
#
        cols = { "Sf_ID"                 : sfid,
                 "Spectral_peak_list_ID" : localid,
                 "Entry_ID"              : "'%s'" % (self._id,),
                 "ID"                    : ":num" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Spectral_transition" ("%s") values (%s)' % (colstr,valstr,)

        for row in table["spectral_transition"] :
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = { "num" : row["id"] }, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        cols = { "Sf_ID"                  : sfid,
                 "Spectral_peak_list_ID"  : localid,
                 "Entry_ID"               : "'%s'" % (self._id,),
                 "Spectral_transition_ID" : ":num",
                 "Intensity_val"          : ":val",
                 "Measurement_method"     : ":meth" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Spectral_transition_general_char" ("%s") values (%s)' % (colstr,valstr,)
        params = {}
        for row in table["spectral_transition_general_char"] :
            params.clear()
            params["num"] = row["spectral_transition_id"]
            params["val"] = row["intensity_val"]
            params["meth"] = row["measurement_method"]
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

        cols = { "Sf_ID"                  : sfid,
                 "Spectral_peak_list_ID"  : localid,
                 "Entry_ID"               : "'%s'" % (self._id,),
                 "Spectral_transition_ID" : ":num",
                 "Spectral_dim_ID"        : ":dim",
                 "Chem_shift_val"         : ":val" }

        colstr = '","'.join( str( k ) for k in cols.keys() )
        valstr = ",".join( str( v ) for v in cols.values() )
        sql = 'insert into "Spectral_transition_char" ("%s") values (%s)' % (colstr,valstr,)
        params = {}
        for row in table["spectral_transition_char"] :
            params.clear()
            params["num"] = row["spectral_transition_id"]
            params["val"] = row["chem_shift_val"]
            params["dim"] = row["spectral_dim_id"]
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
                pprint.pprint( params )
            rc = self._star.execute( sql, params = params, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# some atoms may be labeled
#
        if "assigned_spectral_transition" in table.keys() :
            cols = { "Sf_ID"                  : sfid,
                     "Spectral_peak_list_ID"  : localid,
                     "Entry_ID"               : "'%s'" % (self._id,),
                     "Spectral_transition_ID" : ":num",
                     "Spectral_dim_ID"        : ":dim",
                     "Val"                    : ":val",
                     "Atom_ID"                : ":atm" }

            colstr = '","'.join( str( k ) for k in cols.keys() )
            valstr = ",".join( str( v ) for v in cols.values() )
            sql = 'insert into "Assigned_spectral_transition" ("%s") values (%s)' % (colstr,valstr,)
            params = {}
            for row in table["assigned_spectral_transition"] :
                params.clear()
                params["num"] = row["spectral_transition_id"]
                params["val"] = row["val"]
                params["dim"] = row["spectral_dim_id"]
                params["atm"] = row["atom_id"]
                if self._verbose :
                    sys.stdout.write( sql )
                    sys.stdout.write( "\n" )
                    pprint.pprint( params )
                rc = self._star.execute( sql, params = params, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

# missing pieces
#
            sql = 'update "Assigned_spectral_transition" set "Assigned_chem_shift_list_ID"=1,' \
                + '"Entity_assembly_ID"=1,"Entity_ID"=1,"Comp_index_ID"=1'
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
            rc = self._star.execute( sql, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

            sql = 'update "Assigned_spectral_transition" set "Comp_ID"=(select "ID" from "Chem_comp")'
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
            rc = self._star.execute( sql, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows updated\n" % (rc.rowcount,) )

            sql = 'update "Assigned_spectral_transition" set "Assembly_atom_ID"=' \
                + '(select "Assembly_atom_ID" from "Atom" where "Entity_assembly_ID"=1 ' \
                + 'and "Entity_ID"=1 and "Comp_index_ID"=1 and "Atom_ID"="Assigned_spectral_transition"."Atom_ID")'
            if self._verbose :
                sys.stdout.write( sql )
                sys.stdout.write( "\n" )
            rc = self._star.execute( sql, commit = True )
            if self._verbose :
                sys.stdout.write( "%d rows inserted\n" % (rc.rowcount,) )

    # sanity checks
    #
    def _is_sane( self ) :

# get main inchi string and formula
#
        sql = 'select "InChI_code","Formula" from "Chem_comp"'
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rs = self._star.query( sql )
        inchi = None
        formula = None
        for row in rs :
            inchi = row[0]
            formula = row[1]
        if inchi is None :
            sys.stderr.write( "ERR: InChI_code is NULL\n" )
            return False
        inchi = str( inchi ).strip()
        if inchi == "" :
            sys.stderr.write( "ERR: empty InChI_code\n" )
            return False
        if formula is None :
            sys.stderr.write( "ERR: Formula is NULL\n" )
            return False
        formula = str( formula ).strip()
        if inchi == "" :
            sys.stderr.write( "ERR: empty Formula\n" )
            return False

        if self._verbose :
            sys.stdout.write( "Formula: %s\nInChI: %s\n" % (formula,inchi,) )

        pat = re.compile( r"^InChI=[^/]+/([^/]+)/" )
        m = pat.search( inchi )
        if not m :
            sys.stderr.write( "ERR: InChI_code %s does not match pattern\n" % (inchi,) )
            return False


# this fails for e.g. InChI=1S/C23H39NO19.Na vs C23H38NNaO19
#        if m.group( 1 ) != formula :
#            sys.stderr.write( "ERR: formula in InChI_code %s does not match Formula %s\n" % (m.group( 1 ),formula,) )
#            return False

# now check descriptors
#
        sql = """select "Descriptor","Program" from "Chem_comp_descriptor" where "Type"='InChI'"""
        if self._verbose :
            sys.stdout.write( sql )
            sys.stdout.write( "\n" )
        rs = self._star.query( sql )
        inchis = {}
        for row in rs :
            inchis[row[1]] = row[0]

        if self._verbose :
            pprint.pprint( inchis )

        if not "ALATIS" in inchis.keys() :
            sys.stderr.write( "ERR: no ALATIS InChI in decriptors\n" )
            return False

        if inchis["ALATIS"] != inchi :
            sys.stderr.write( "ERR: ALATIS InChI in decriptors does not match InChI_code\n" )
            return False

        sql = 'delete from "Chem_comp_descriptor" where "Program"=:prog and "Descriptor"=:inc'
        for prog in inchis.keys() :
            if prog == "ALATIS" : continue

            if inchis[prog] == inchis["ALATIS"] :

                if self._verbose :
                    sys.stdout.write( "%s %s %S\n" % (sql,prog,inchis[prog],) )
                rc = self._star.execute( sql, params = { "prog" : prog, "inc" : inchis[prog] }, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows deleted\n" % (rc.rowcount,) )

            m = pat.search( inchis[prog] )
            if not m :

# just get rid of it
#
                sys.stderr.write( "ERR: %s InChI does not match pattern, deleting\n" % (prog,) )
                if self._verbose :
                    sys.stdout.write( "%s %s %S\n" % (sql,prog,inchis[prog],) )
                rc = self._star.execute( sql, params = { "prog" : prog, "inc" : inchis[prog] }, commit = True )
                if self._verbose :
                    sys.stdout.write( "%d rows deleted\n" % (rc.rowcount,) )

            if m.group( 1 ) != formula :
                sys.stderr.write( "WARN: %s InChI does not match formula\n" % (prog,) )
                if self._verbose :
                    sys.stdout.write( "%s %s %S\n" % (sql,prog,inchis[prog],) )
#                rc = self._star.execute( sql, params = { "prog" : prog, "inc" : inchis[prog] }, commit = True )
#                if self._verbose :
#                    sys.stdout.write( "%d rows deleted\n" % (rc.rowcount,) )

        return True

#
#
#
if __name__ == '__main__':

    verbose = False
    cp = ConfigParser.SafeConfigParser()
    cp.read( sys.argv[1] )

    dat = None
    with open( sys.argv[2], "rU" ) as fin :
        dat = json.load( fin )

    e = StarMaker.from_nmrfam( config = cp, data = dat, id = "temp123", verbose = verbose )
    e.to_star( out = sys.stdout, err = sys.stderr )

#
#
