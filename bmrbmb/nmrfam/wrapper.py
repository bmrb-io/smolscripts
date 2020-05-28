#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  nmrfam_to_bmrb.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit

# convert files created by NMRFAM to NMR-STAR-ish data struct
# that is later merged with boilerplate stuff, inflated to fit BMRB completeness requirements,
# and converted to NMR-STAR
#
# IN:
#
# 1. catalog file: JSON dict, as of 20180801 contains
# "cid"          : PubChem CID (if available),
# "vendor"       : e.g. "Sigma-Aldrich",
# "catalog_num"  : vendor's ID,
# "catalog_name" : vendor's name,
# "solvent"      : one of the ones in sample/__init__.py, e.g. "methanol"
#    (NOTE: mixed solvents need special treatment)
# "paramagnetic" : "no",
# "sdf"          : ALATIS'ed SDF,
# "timedomain"   : name of the tarball file w/ TopSpin expt. directories,
# "botlog"       : name of nmrbot log file,
# "shifts"       : name of assignments file (CSV),
# "molpic"       : name of image file w/ atom labels,
# "inchi"        : ALATIS InChI string,
# "PEAKFILES"    : a ( filename : experiment name } dict of peak list filenames
# "PNGFILES"     : a { filename : experiment name } dict of spectra image filenames
#
# Experiment names are in topspin/__init__.py
#
# 2. all files referenced in the above catalog.
#
# OUT:
#
# a list of saveframe categories, each is
#  a list of saveframes, each is
#    a dict of "free tags" and "loop tables"
#      each loop table is a list of dicts (rows)
#
# get it as .data or as .json
#

from __future__ import absolute_import

import os
import sys
import json
import re
import pprint

# self
#
_HERE = os.path.split( __file__ )[0]
sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
import bmrbmb
import bmrbmb.boilerplate
import bmrbmb.chemcomp
import bmrbmb.www
import bmrbmb.sample
import bmrbmb.topspin
import bmrbmb.chemshifts

#
#
#
class FromNMRFAM( object ) :


    # main
    # indexfile is the catalog described above
    #
    @classmethod
    def make_entry( cls, indexfile, data_name = "set01", verbose = False ) :

        infile = os.path.realpath( indexfile )
        if not os.path.exists( infile ) :
            raise IOError( "File not found: %s" % (infile,) )
        indir = os.path.split( infile )[0]

        with open( infile, "rU" ) as f :
            idx = json.load( f )

        rc = cls( index = idx, workdir = indir, data_name = data_name, verbose = verbose )
        rc.precheck()

        boilerplate = True

        if boilerplate :
            rc._make_entry_information()
            rc._make_citations()

            rc._make_natural_source()
            rc._make_experimental_source()

        rc._make_chem_comp()
        rc._make_sample()

        if boilerplate :
            rc._make_software()

# FIXME (see comment in the method)
#
        rc._make_spectrometer()

        rc._parse_experiments()
        rc._parse_topspin_tarball()
        rc._make_experiments()
        rc._make_cs_reference()
        rc._make_cs_assignments()
        rc._make_peak_lists()

        return rc

    # index is a dictionary read from json index file
    # workdir is where all the files are
    # storage is a list of lists.
    # data_name is "set01", "set02", etc.
    #  nested lists are lists of saveframes, unique ones: with only one element.
    #  this shouldn't grow too big on small molecule entries.
    #
    def __init__( self, index, workdir, data_name = "set01", verbose = False ) :
        self._idx = index
        self._indir = os.path.realpath( workdir )
        if not os.path.isdir( workdir ) :
            raise IOError( "Not a directory: %s" % (self._indir,) )
        if data_name is None :
            self._set_num = ""
        else :
            self._set_num = str( data_name ).strip()
        self._verbose = bool( verbose )
        self._dat = []

# chem comp is created as part of a pre-check so we don't do this twice
#
        self._mol = None
        self._chem_comp = None
        self._inchi = None
        self._expts = {}
        self._exptfiles = {}

        self._star = None

    # entry data
    #
    @property
    def data( self ) :
        if len( self._dat ) < 1 : return None
        return self._dat

    def __json__( self ) :
        if len( self._dat ) < 1 : return json.dumps( None )
        return json.dumps( self._dat )
    json = property( __json__ )

    ################################################################################################
    # data processing part

    ################################################################################################
    # check that we have everything we need
    #
    def precheck( self ) :
        if self._verbose : sys.stdout.write( "%s.precheck()\n" % (self.__class__.__name__,) )

        if not "sdf" in self._idx.keys() :
            raise Exception( "No sdf in index file" )


        sdf = os.path.realpath( os.path.join( self._indir, self._idx["sdf"] ) )
        if not os.path.exists( sdf ) :
            raise IOError( "File not found: %s" % (sdf,) )

# sample
#
        if not "solvent" in self._idx.keys() :
            raise Exception( "No solvent in index file" )
        if not "catalog_name" in self._idx.keys() :
            raise Exception( "No catalog_name in index file" )
        if not all( ord( c ) < 128 for c in self._idx["catalog_name"] ) :
            raise Exception( "Non-ASCII byte in catalog_name in index file" )
        if not "catalog_num" in self._idx.keys() :
            raise Exception( "No catalog_num in index file" )
        if not "vendor" in self._idx.keys() :
            raise Exception( "No vendor in index file" )

        if not "solute_amount" in self._idx.keys() :
            sys.stderr.write( "WARN: no solute amount in index file, assumed 100\n" )
            self._idx["solute_amount"] = 100
        if not "solute_units" in self._idx.keys() :
            sys.stderr.write( "WARN: no solute_units in index file, assumed mM\n" )
            self._idx["solute_units"] = "mM"

# nmrbot log is needed for experiment lists etc.
#
        if not "botlog" in self._idx.keys() :
            raise Exception( "No botlog in index file" )

        botlog = os.path.realpath( os.path.join( self._indir, self._idx["botlog"] ) )
        if not os.path.exists( botlog ) :
            raise IOError( "File not found: %s" % (botlog,) )

# topspin tarball is needed for experiment files
#
        if not "timedomain" in self._idx.keys() :
            raise Exception( "No timedomain in index file" )

        tarball = os.path.realpath( os.path.join( self._indir, self._idx["timedomain"] ) )
        if not os.path.exists( tarball ) :
            raise IOError( "File not found: %s" % (tarball,) )

# CSV w/ assignments
#
        if not "shifts" in self._idx.keys() :
            raise Exception( "No shifts in index file" )

        table = os.path.realpath( os.path.join( self._indir, self._idx["shifts"] ) )
        if not os.path.exists( table) :
            raise IOError( "File not found: %s" % (table,) )

# SVG w/ atom labels
#
        if not "molpic" in self._idx.keys() :
            raise Exception( "No molpic in index file" )

        pic = os.path.realpath( os.path.join( self._indir, self._idx["molpic"] ) )
        if not os.path.exists( pic ) :
            raise IOError( "File not found: %s" % (molpic,) )

# peak files
# I suppose we can live w/o those
#
        if not "PEAKFILES" in self._idx.keys() :
            sys.stderr.write( "WARN: no PEAKFILES in index file\n" )
            self._idx["PEAKFILES"] = {}
        for f in self._idx["PEAKFILES"].keys() :
            fname = os.path.realpath( os.path.join( self._indir, f ) )
            if not os.path.exists( fname ) :
                raise IOError( "File not found: %s" % (fname,) )
            if not self._idx["PEAKFILES"][f] in bmrbmb.topspin.EXPERIMENTS.values() :
                raise Exception( "Unknown experiment: %s for file %s" % (expt,fname,) )

# spectra images
# ditto
#
        if not "PNGFILES" in self._idx.keys() :
            sys.stderr.write( "WARN: no PNGFILES in index file\n" )
            self._idx["PNGFILES"] = {}
        for f in self._idx["PNGFILES"].keys() :
            fname = os.path.realpath( os.path.join( self._indir, f ) )
            if not os.path.exists( fname ) :
                raise IOError( "File not found: %s" % (fname,) )
            if not self._idx["PNGFILES"][f] in bmrbmb.topspin.EXPERIMENTS.values() :
                raise Exception( "Unknown experiment: %s for file %s" % (expt,fname,) )

# chem comp from sdf
#

        self._mol = bmrbmb.chemcomp.Molecule.from_file( filename = sdf, verbose = self._verbose )
        self._chem_comp = self._mol.chem_comp

# inchi string is needed for web queries
#
        if "descriptors" in self._chem_comp.keys() :
            for i in self._chem_comp["descriptors"] :
                if (i["type"] == "InChI") and (i["program"] == "ALATIS") :
                    self._inchi = i["descriptor"]
        if self._inchi is None :
            if not "inchi" in self._idx.keys() :
                raise Exception( "No InChI string in SDF or index file" )
            self._inchi = self._idx["inchi"]

# top part is one-liners
#
#
    ################################################################################################
    #
    def _make_entry_information( self ) :
        if self._verbose : sys.stdout.write( "%s.make_entry_information()\n" % (self.__class__.__name__,) )
        self._dat.append( bmrbmb.boilerplate.get_saveframe( "entry_information" ) )

    ################################################################################################
    #
    def _make_citations( self ) :
        if self._verbose : sys.stdout.write( "%s.make_citations()\n" % (self.__class__.__name__,) )
        self._dat.append( bmrbmb.boilerplate.get_saveframe( "citation" ) )

# assembly and entity are made from chem comp when making star file
#
    ################################################################################################
    #
    def _make_natural_source( self ) :
        if self._verbose : sys.stdout.write( "%s.make_natural_source()\n" % (self.__class__.__name__,) )
        self._dat.append( bmrbmb.boilerplate.get_saveframe( "natural_source" ) )

    ################################################################################################
    #
    def _make_experimental_source( self ) :
        if self._verbose : sys.stdout.write( "%s.make_experimental_source()\n" % (self.__class__.__name__,) )
        self._dat.append( bmrbmb.boilerplate.get_saveframe( "experimental_source" ) )

# chem comp is missing a few bits
#
    ################################################################################################
    # num is always 1 here
    #
    def _make_chem_comp( self, num = 1 ) :
        if self._verbose : sys.stdout.write( "%s.make_chem_comp()\n" % (self.__class__.__name__,) )

        if self._chem_comp is None :
            raise Exception( "Run precheck() first!" )

        self._chem_comp["name"] = "chem_comp_%s" % (num,)
        self._chem_comp["inchi_code"] = self._inchi
        self._chem_comp["paramagnetic"] = self._idx["paramagnetic"]

# web stuff
#
        if not "dblinks" in self._chem_comp.keys() :
            self._chem_comp["dblinks"] = []
        bmrbid = bmrbmb.www.bmrbid( self._inchi, self._verbose )
        if not bmrbid is None :
            self._chem_comp["dblinks"].append( bmrbid )

        if not "identifiers" in self._chem_comp.keys() :
            self._chem_comp["identifiers"] = []
        iupac_name = bmrbmb.www.cactus_iupac_name( self._inchi, self._verbose )
        if not iupac_name is None :
            for name in iupac_name :
                self._chem_comp["identifiers"].append( name )

        if "cid" in self._idx.keys() :
            pug = bmrbmb.www.get_pubchem_data( self._idx["cid"], self._verbose )
            if pug is not None :
                if "identifiers" in pug.keys() :
                    for i in pug["identifiers"] :
                        self._chem_comp["identifiers"].append( i )
                if "common_name" in pug.keys() :
                    if not "common_name" in self._chem_comp.keys() :
                        self._chem_comp["common_name"] = []
                    for i in pug["common_name"] :
                        self._chem_comp["common_name"].append( i )

# make it a list of one
#
        self._dat.append( [self._chem_comp] )

    ################################################################################################
    # ALATIS map: for entries made before/without ALATIS
    # atom_ids throughout the entry are pre-alatis
    # atom_name is alatis and nomenclature is ALATIS
    #
    def add_alatis_map( self ) :
        if self._verbose : sys.stdout.write( "%s.add_alatis_map()\n" % (self.__class__.__name__,) )

        if self._chem_comp is None :
            raise Exception( "Run precheck() first!" )

        alt = self._mol.alatis_map()
        if alt is None :
            return

        if not "atom_nomenclature" in self._chem_comp.keys() :
            self._chem_comp["atom_nomenclature"] = []

        for row in alt :
            self._chem_comp["atom_nomenclature"].append( { "atom_id" : row["atom_id"],
                            "atom_name" : row["alatis_id"], "naming_system" : "ALATIS" } )

        return

# sample & conditions
#
    ################################################################################################
    #
    def _make_sample( self ) :
        if self._verbose : sys.stdout.write( "%s.make_sample()\n" % (self.__class__.__name__,) )

        sm = bmrbmb.sample.Sample.create( solute = self._idx["catalog_name"],
                                        solvent = self._idx["solvent"],
                                        vendor = self._idx["vendor"],
                                        catalog_name = self._idx["catalog_name"],
                                        catalog_num = self._idx["catalog_num"],
                                        concentration = self._idx["solute_amount"],
                                        units = self._idx["solute_units"] )

        self._dat.append( sm.data )

    ################################################################################################
    # software is boilerplate
    #
    def _make_software( self ) :
        if self._verbose : sys.stdout.write( "%s.make_software()\n" % (self.__class__.__name__,) )
        self._dat.append( bmrbmb.boilerplate.get_saveframe( "software" ) )

# FIXME
# there is only one spectrometer
# this should be extended to allow a choice and/or multiple instruments as well.
# but that would also require linking experiments to the instruments and dealing with multiple
# nmrbot log files and timedoamin data sets.
#
    ################################################################################################
    #
    def _make_spectrometer( self ) :
        if self._verbose : sys.stdout.write( "%s.make_spectrometer()\n" % (self.__class__.__name__,) )
        self._dat.append( bmrbmb.boilerplate.get_saveframe( "spectrometer" ) )

    ################################################################################################
    # read experiment list from NMRBot log file
    # this also reads stuff that goes into peak list saveframes.
    # experiment list, OTOH, needs filenames from another file
    #
    def _parse_experiments( self ) :
        if self._verbose : sys.stdout.write( "%s.parse_experiments()\n" % (self.__class__.__name__,) )

        botlog = os.path.realpath( os.path.join( self._indir, self._idx["botlog"] ) )
        bl = bmrbmb.topspin.NMRBotLog.parse( filename = botlog, verbose = self._verbose )
        self._expts = bl.data

    ################################################################################################
    # read experiment files from tarball
    #
    def _parse_topspin_tarball( self ) :
        if self._verbose : sys.stdout.write( "%s.parse_topspin_tarball()\n" % (self.__class__.__name__,) )
        tarball = os.path.realpath( os.path.join( self._indir, self._idx["timedomain"] ) )
        files = bmrbmb.topspin.ExperimentFiles.read( tarfile = tarball, verbose = self._verbose )
        self._exptfiles = files.data

    ################################################################################################
    # experiment_list saveframe
    #
    def _make_experiments( self ) :
        if self._verbose : sys.stdout.write( "%s._make_experiments()\n" % (self.__class__.__name__,) )

        if len( self._expts ) < 1 :
            self._parse_experiments()
        if len( self._exptfiles ) < 1 :
            self._parse_topspin_tarball()
        if (len( self._expts ) < 1) or (len( self._exptfiles ) < 1) :
            raise Exception( "No experiments/experiment directories!" )

# add missing bits
#
        for i in self._expts.keys() :
            for j in self._exptfiles.keys() :
                if (str( i ) == str( self._exptfiles[j]["experiment_id"] )) \
                and (self._exptfiles[j]["content"] == "NMR experiment directory") :
                    self._expts[i]["raw_data_flag"] = "yes"
                    self._exptfiles[j]["details"] = self._expts[i]["name"]
                    break

# expts
#
        sf = { "sf_category" : "experiment_list", "name" : "experiment_list" }
        sf["experiment"] = []
        for i in sorted( self._expts.keys(), cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :
            if "raw_data_flag" in self._expts[i].keys() :
                sf["experiment"].append( { "id" : i, "name" : self._expts[i]["name"], "raw_data_flag" : "yes" } )
            else :
                sf["experiment"].append( { "id" : i, "name" : self._expts[i]["name"], "raw_data_flag" : "no" } )

# files: topspin
#
        sf["experiment_file"] = []
        for i in sorted( self._exptfiles.keys(), cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :
            sf["experiment_file"].append( self._exptfiles[i] )

# files: peaklists
# FIXME: type text/xml is hardcoded here
#
        for (f, e) in self._idx["PEAKFILES"].items() :
            for i in self._expts.keys() :
                if self._expts[i]["name"] == e :
                    sf["experiment_file"].append( { "experiment_id" : i, "name" : f, "type" : "text/xml" ,
                            "details" : "TopSpin XML", "content" : "Peak list" } )
                    break

# files: pix
# FIXME: type image/png is hardcoded here
#
        for (f, e) in self._idx["PNGFILES"].items() :
            for i in self._expts.keys() :
                if self._expts[i]["name"] == e :
                    sf["experiment_file"].append( { "experiment_id" : i, "name" : f, "type" : "image/png" ,
                            "content" : "Spectral image" } )
                    break

# there's only one but we return the list of one for consistency
#
        self._dat.append( [sf] )

    ################################################################################################
    # pick the right one for reference compound
    #
    def _make_cs_reference( self ) :
        if self._verbose : sys.stdout.write( "%s.make_cs_reference()\n" % (self.__class__.__name__,) )

        ref = None
        for i in self._dat :
            for j in i :
                if "sample_component" in j.keys() :
                    for k in j["sample_component"] :
                        if k["type"] == "reference" :
                            ref = k["mol_common_name"]
                            break

        if ref is None :
            raise Exception( "Can't find reference compound in sample" )

        rc = None
        refs = bmrbmb.boilerplate.get_saveframe( "chem_shift_reference" )
        for i in refs :
            if i["id"] == ref :
                rc = i

        if rc is None :
            raise Exception( "No CS reference for %s" % (ref,) )

# override the default name
#
        rc["name"] = "chem_shift_ref_%s" % (self._set_num,)
        self._dat.append( [rc] )

    ################################################################################################
    #
    def _parse_assignments( self ) :
        if self._verbose : sys.stdout.write( "%s.parse_assignments()\n" % (self.__class__.__name__,) )

        csv = os.path.realpath( os.path.join( self._indir, self._idx["shifts"] ) )
        cs = bmrbmb.chemshifts.AtomChemShift.readfile( filename = csv, num = self._set_num, verbose = self._verbose )
        return cs.data

# experiments are all hardcoded to be "all"
#

    ################################################################################################
    #
    #
    def _add_cs_experiments( self, saveframe ) :
        if self._verbose : sys.stdout.write( "%s.all_cs_experiments()\n" % (self.__class__.__name__,) )

        if len( self._expts ) < 1 :
            self._parse_experiments()

        saveframe["chem_shift_experiment"] = []
        for i in sorted( self._expts.keys(), cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :
            saveframe["chem_shift_experiment"].append( { "experiment_id" : i, "experiment_name" : self._expts[i]["name"] } )

    ################################################################################################
    # read assignments from csv, add experiments
    #  note: parse csv returns a single dict, we add experiments withut worrying about adding to the
    #  right saveframe. then wrap in a list of one.
    #
    def _make_cs_assignments( self ) :
        if self._verbose : sys.stdout.write( "%s.make_cs_assignments()\n" % (self.__class__.__name__,) )

        dat = self._parse_assignments()
        self._add_cs_experiments( saveframe = dat )
        self._dat.append( [dat] )

    ################################################################################################
    # build list of peak saveframe dicts
    #
    def _make_peak_lists( self ) :
        if self._verbose : sys.stdout.write( "%s.make_peak_lists()\n" % (self.__class__.__name__,) )

        rc = []
        for f in self._idx["PEAKFILES"].keys() :
            fname = os.path.realpath( os.path.join( self._indir, f ) )
            if not os.path.exists( fname ) :
                raise IOError( "File not found: %s" % (fname,) )
            expt = self._idx["PEAKFILES"][f]
            if not expt in bmrbmb.topspin.EXPERIMENTS.values() :
                raise Exception( "Unknown experiment: %s" % (expt,) )

            rc.append( self._add_peak_list( exptname = expt, filename = fname ) )

        if len( rc ) > 0 :
            self._dat.append( rc )

    ################################################################################################
    # return a "saveframe" dict
    # could propably live without these and print warnings instead of throwing exceptions below
    #
    def _add_peak_list( self, exptname, filename ) :
        if self._verbose : sys.stdout.write( "%s._add_peak_lists(%s,%s)\n" % (self.__class__.__name__,exptname,filename,) )

        rc = { "sf_category" : "spectral_peak_list",
            "name" : "spectral_peaks_%s_%s" % (re.sub( r"\s+", "_", exptname ),self._set_num,),
            "experiment_name" : exptname }

        dim = []
        cnt = 0
        dtl = ""

# spectral_dim and top-level detail is frem experiments
#
        for i in self._expts.keys() :
            if self._expts[i]["name"] == exptname :
                rc["experiment_id"] = str( i ).strip()
                rc["experiment_name"] = exptname
                for j in sorted( self._expts[i]["dims"].keys() ) :
                    m = re.search( r"(\d+)([A-Za-z]+)", self._expts[i]["dims"][j]["nuc"] )
                    if not m :
                        raise Exception( "Bad nucleus %s in %s" % (self._expts[i]["dims"][j]["nuc"],filename,) )
                    cnt += 1

# there's nowhere to put topspin's dimension labels: F1, F2, etc., so we put them in details at top level
#
                    if dtl == "" : dtl = "%s: %s" % (j, self._expts[i]["dims"][j]["nuc"])
                    else : dtl = "%s, %s: %s" % (dtl, j, self._expts[i]["dims"][j]["nuc"])

                    sd = None
                    if m.group( 2 ) == "H" : sd = "Full H"
                    elif m.group( 2 ) == "C" : sd = "Full C"
                    else :
                        raise Exception( "Unsupported atom type %s in %s" % (m.group( 2 ),filename,) )

# FIXME ppm, 2 dimensions, 'relaitive height' are hardcoded
#
                    dim.append( { "id" : cnt, "atom_type" : m.group( 2 ), "atom_isotope_number" : m.group( 1 ),
                                "spectral_region" : sd, "sweep_width" : self._expts[i]["dims"][j]["sw"],
                                "sweep_width_units" : "ppm" } )

                break

        if len( dim ) < 1 :
            raise Exception( "No spectral dimensions in %s" % (filename,) )

        rc["details"] = dtl
        rc["number_of_spectral_dimensions"] = cnt
        rc["spectral_dim"] = dim

# the peaks
#
        trans = []
        trans_gen = []
        trans_char = []
        trans_ass = []

        peaks = bmrbmb.topspin.TopSpinPeakList.parse( filename = filename, verbose = self._verbose )
        if peaks.data is None :
            raise Exception( "No transitions in %s" % (filename,) )

        if str( peaks.data["experiment"] ).strip() != rc["experiment_id"] :

            sys.stderr.write( "**** peaks ****\n" )
            pprint.pprint( peaks.data, stream = sys.stderr )
            sys.stderr.write( "**** experiments ****\n" )
            pprint.pprint( self._expts, stream = sys.stderr )
            raise Exception( "Expt ID in peaks %s (%s) does not match id from index %s" \
                    % (filename,peaks.data["experiment"],rc["experiment_id"],) )

        for i in peaks.data["peaks"].keys() :
            trans.append( { "id" : str( i ).strip() } )
            trans_gen.append( { "spectral_transition_id" : str( i ).strip(),
                    "intensity_val" : peaks.data["peaks"][i]["intensity"],
                    "measurement_method" : "relative height" } )

            trans_char.append( { "spectral_transition_id" : str( i ).strip(),
                    "spectral_dim_id" : 1,
                    "chem_shift_val" : peaks.data["peaks"][i]["shift1"] } )
# 2D?
#
            if "shift2" in peaks.data["peaks"][i].keys() :
                trans_char.append( { "spectral_transition_id" : str( i ).strip(),
                        "spectral_dim_id" : 2,
                        "chem_shift_val" : peaks.data["peaks"][i]["shift2"] } )

# has atom label?
# there is never an atom label on 2D peaks, as of the time of this writing
#
            if "atom" in peaks.data["peaks"][i].keys() :
                trans_ass.append( { "spectral_transition_id" : str( i ).strip(),
                        "spectral_dim_id" : 1,
                        "val" : peaks.data["peaks"][i]["shift1"],
                        "atom_id" : peaks.data["peaks"][i]["atom"] } )

        if (len( trans ) < 1) or (len( trans_gen ) < 1) or (len( trans_char ) < 1) :
            raise Exception( "No transitions in %s" % (filename,) )

        rc["spectral_transition"] = trans
        rc["spectral_transition_general_char"] = trans_gen
        rc["spectral_transition_char"] = trans_char
        if len( trans_ass ) > 0 :
            rc["assigned_spectral_transition"] = trans_ass

# JIC
#
        peaks = None

        return rc

#
#
#
if __name__ == '__main__':

    dat = FromNMRFAM.make_entry( indexfile = sys.argv[1], verbose = False )
    dat.add_alatis_map()
    sys.stdout.write( json.dumps( dat.data, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )


#
#
