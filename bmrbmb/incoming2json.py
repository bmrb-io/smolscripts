#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  incoming2json.py
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
import re
import pprint

# fscking python can't import without "from ." in a package, and with "from ." outside of it.
# so don't run __main__ without editing this first:
#
if __name__ != '__main__':
    from . import boilerplate
    from . import chemcomp
    from . import chemshifts
    from . import nmrfam
    from . import nmrstar
    from . import sample
    from . import topspin
    from . import www

"""
Read files in incoming directory and make a partial BMRB entry.
This creates a data structure that can be dumped into JSON or passed
to NMR-STAR builder to create a complete BMRB entry NMR-STAR file.

This version assumes multiple samples but same conditions, spectrometer, etc.
"""

#
#
#
class FAMtoJSN( object ) :

# datasets are indexed 1, 2, ... in the input
# this is more human-readable.
# NOTE only 99 datasets/entry.
#
    DATASET_MASK = "set%02d"

    # main
    # indexfile is the catalog file, see nmrfam2incoming.py
    #
    @classmethod
    def make_entry( cls, indexfile, verbose = False ) :

        infile = os.path.realpath( indexfile )
        if not os.path.exists( infile ) :
            raise IOError( "File not found: %s" % (infile,) )
        indir = os.path.split( infile )[0]

        with open( infile, "rU" ) as f :
            idx = json.load( f )

        rc = cls( index = idx, workdir = indir, verbose = verbose )
        rc.make_dat( boilerplate = True )

        return rc

    # index is a dictionary read from json index file
    # workdir is where all the files are
    # storage is a list of lists.
    #  nested lists are lists of saveframes, unique ones: with only one element.
    #  this shouldn't grow too big on small molecule entries.
    #
    def __init__( self, index, workdir, verbose = False ) :
        self._idx = index
        self._indir = os.path.realpath( workdir )
        if not os.path.isdir( workdir ) :
            raise IOError( "Not a directory: %s" % (self._indir,) )
        self._verbose = bool( verbose )
        self._dat = []
        self._mol = None
        self._chem_comp = None
        self._inchi = None

# need to keep track of these
#
        self._samples = {}
        self._conditions = {}
        self._csrs = {}
        self._spectrometers = {}
        self._experiments = {}
        self._shifts = {}

    # entry data
    #
    @property
    def data( self ) :
        if len( self._dat ) < 1 : return None
        return self._dat

    #
    #
    @property
    def __json__( self ) :
        return json.dumps( self._dat, indent = 4, sort_keys = True )

    ################################################################################################
    # data processing part
    #
    # boilerplate is a flag the triggers generation of BMRB cruft -- if false, this only processes
    #  real data from NMRFAM
    #
    def make_dat( self, boilerplate = True ) :
        if self._verbose : sys.stdout.write( "%s.make_dat()\n" % (self.__class__.__name__,) )

        if boilerplate :
            self._make_entry_information()
            self._make_citations()

            self._make_natural_source()
            self._make_experimental_source()

            self._make_software()

        self._make_chem_comp()

# FIXME (see comment in the method)
#
        self._make_spectrometer()

# order is impotant here: make_samples() creates map for cs_ref
#
        self._make_samples()
        self._make_sample_conditions()
        self._make_cs_reference()

# make_experients() creates map for CS and peaks
#
        self._make_experiments()
        self._make_cs_assignments()
        self._make_peak_lists()

# one-liners
#
    ################################################################################################
    #
    def _make_entry_information( self ) :
        if self._verbose : sys.stdout.write( "%s.make_entry_information()\n" % (self.__class__.__name__,) )
        self._dat.append( boilerplate.get_saveframe( "entry_information" ) )
    #
    #
    def _make_citations( self ) :
        if self._verbose : sys.stdout.write( "%s.make_citations()\n" % (self.__class__.__name__,) )
        self._dat.append( boilerplate.get_saveframe( "citation" ) )

# assembly and entity are made from chem comp when making star file
#
    #
    def _make_natural_source( self ) :
        if self._verbose : sys.stdout.write( "%s.make_natural_source()\n" % (self.__class__.__name__,) )
        self._dat.append( boilerplate.get_saveframe( "natural_source" ) )
    #
    #
    def _make_experimental_source( self ) :
        if self._verbose : sys.stdout.write( "%s.make_experimental_source()\n" % (self.__class__.__name__,) )
        self._dat.append( boilerplate.get_saveframe( "experimental_source" ) )
    #
    #
    def _make_software( self ) :
        if self._verbose : sys.stdout.write( "%s.make_software()\n" % (self.__class__.__name__,) )
        self._dat.append( boilerplate.get_saveframe( "software" ) )

# FIXME
# there is only one spectrometer
# this should be extended to allow a choice and/or multiple instruments as well.
#
    #
    def _make_spectrometer( self ) :
        if self._verbose : sys.stdout.write( "%s.make_spectrometer()\n" % (self.__class__.__name__,) )
        instrument = boilerplate.get_saveframe( "spectrometer" )
        self._dat.append( instrument )
        for i in self._idx["data"].keys() :
            self._spectrometers[i] = instrument[0]["name"]

#
# chem comp is missing a few bits
#
    ################################################################################################
    # num is always 1 here
    #
    def _make_chem_comp( self, num = 1 ) :
        if self._verbose : sys.stdout.write( "%s.make_chem_comp()\n" % (self.__class__.__name__,) )

        sdf = os.path.join( self._indir, self._idx["sdf"] )
        self._mol = chemcomp.Molecule.from_file( filename = sdf, verbose = self._verbose )
        self._chem_comp = self._mol.chem_comp

        if "descriptors" in self._chem_comp.keys() :
            for i in self._chem_comp["descriptors"] :
                if (i["type"] == "InChI") and (i["program"] == "ALATIS") :
                    self._inchi = i["descriptor"]
        if self._inchi is None :
            if not "inchi" in self._idx.keys() :
                raise Exception( "No InChI string in SDF or index file" )
            self._inchi = self._idx["inchi"]


        self._chem_comp["name"] = "chem_comp_%s" % (num,)
        self._chem_comp["inchi_code"] = self._inchi
        self._chem_comp["paramagnetic"] = self._idx["paramagnetic"]
        if "vendor" in self._idx.keys() :
            self._chem_comp["vendor"] = self._idx["vendor"]
            self._chem_comp["vendor_product_code"] = self._idx["catalog_num"]
        if "catalog_name" in self._idx.keys() :
            self._chem_comp["name"] = self._idx["catalog_name"]
            if not "common_name" in self._chem_comp.keys() :
                self._chem_comp["common_name"] = []
            self._chem_comp["common_name"].append( self._idx["catalog_name"] )

# web stuff
#
        if not "dblinks" in self._chem_comp.keys() :
            self._chem_comp["dblinks"] = []
        bmrbid = www.bmrbid( self._inchi, self._verbose )
        if not bmrbid is None :
            self._chem_comp["dblinks"].append( bmrbid )

        if not "identifiers" in self._chem_comp.keys() :
            self._chem_comp["identifiers"] = []
        iupac_name = www.cactus_iupac_name( self._inchi, self._verbose )
        if not iupac_name is None :
            for name in iupac_name :
                self._chem_comp["identifiers"].append( name )

        if "cid" in self._idx.keys() :
            pug = www.get_pubchem_data( self._idx["cid"], self._verbose )
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
    #
    # the flip side of the above. Original version is for when ALATIS atom names are
    # different from what's in the entry. This is for when the entry file has ALATIS
    # names but other data files may have other atom names.
    #
    def add_author_atomname_map( self ) :
        if self._verbose : sys.stdout.write( "%s.add_author_atomname_map()\n" % (self.__class__.__name__,) )

        if self._chem_comp is None :
            raise Exception( "Run precheck() first!" )

        alt = self._mol.alatis_map()
        if alt is None :
            return

        if not "atom_nomenclature" in self._chem_comp.keys() :
            self._chem_comp["atom_nomenclature"] = []

        for row in alt :
            self._chem_comp["atom_nomenclature"].append( { "atom_id" : row["alatis_id"],
                            "atom_name" : row["atom_id"], "naming_system" : "author" } )

        return

    ################################################################################################
    # non-unique data.
    # Currently based on samples only.
    #
    # map dataset to
    # - sample conditions
    # - sample
    # - cs reference
    # for experiments and assigned CS.
    # In addition, spectral peaks link to
    # - assigned chemical shifts
    #
    ################################################################################################
    # Need to know if reference compound is TMS or DSS
    #
    def _make_samples( self ) :
        if self._verbose : sys.stdout.write( "%s.make_samples()\n" % (self.__class__.__name__,) )
        samples = []
        for i in self._idx["data"].keys() :
            sm = sample.Sample.create( num = i,
                    solute = self._idx["catalog_name"],
                    vendor = self._idx["vendor"],
                    catalog_name = self._idx["catalog_name"],
                    catalog_num = self._idx["catalog_num"],
                    solvent = self._idx["data"][i]["solvent"],
                    concentration = self._idx["data"][i]["solute_amount"],
                    units = self._idx["data"][i]["solute_units"] )

            samples.append( sm.data )
            self._samples[i] = sm.data["name"]
            for j in sm.data["sample_component"] :
                if j["type"] == "reference" :
                    self._csrs[i] = j["mol_common_name"]
        self._dat.append( samples )

    ################################################################################################
    # D2O and water include pH, others don't.
    # Moreover if we have D2O and water, they may have different pH.
    # cheat: since pH units are different for water and D2O, just make separate conditions
    # and pretend that makes sense.
    #
    def _make_sample_conditions( self ) :
        if self._verbose : sys.stdout.write( "%s.make_sample_conditions()\n" % (self.__class__.__name__,) )

        conditions = {}
        for i in self._idx["data"].keys() :
            if self._idx["data"][i]["solvent"].lower() == "water" :
                if "pH" in self._idx["data"][i] :
                    conditions[i] = "water_ph_%s" % self._idx["data"][i]["pH"]
                else :
                    conditions[i] = "water"
            elif self._idx["data"][i]["solvent"].lower() == "d2o" :
                if "pH" in self._idx["data"][i] :
                    conditions[i] = "d2o_ph_%s" % self._idx["data"][i]["pH"]
                else :
                    conditions[i] = "d2o"
        else :
            conditions[i] = "no_ph"

# saveframes we need to make
#
        conds = set( conditions.values() )
        pat = re.compile( r"^((?:water)|(?:d2o))_ph_([0-9.]+)$" )
        lst = []
        i = 0
        for j in conds :
            i += 1
            if j in ( "d2o", "water", "no_ph" ) :
                sc = sample.Conditions.create( num = i, solvent = j, verbose = self._verbose )
            else :
                m = pat.search( j )
                if not m :
                    raise Exception( "Bad solvent %s" % (j,) )
                    sc = sample.Conditions.create( num = i, solvent = m.group( 1 ), pH = m.group( 2 ), verbose = self._verbose )

            lst.append( sc.data )

            for k in conditions.keys() :
                if conditions[k] == j :
                    self._conditions[k] = sc.data["name"]
                    break

        self._dat.append( lst )

    ################################################################################################
    # pick boilerplate SF for reference compound
    # update the dataset-to-reference mapping
    # *The mapping must already exist*
    #
    def _make_cs_reference( self ) :
        if self._verbose : sys.stdout.write( "%s.make_cs_reference()\n" % (self.__class__.__name__,) )

        assert len( self._csrs ) > 0

        refs = []
        for sf in boilerplate.get_saveframe( "chem_shift_reference" ) :
            for i in self._csrs.keys() :
                if sf["id"] == self._csrs[i] :
                    sfname = "chem_shift_reference_%s" % (self._csrs[i],)
                    sf["name"] = sfname
                    self._csrs[i] = sfname
                    refs.append( sf )
                    break

        self._dat.append( refs )

    ################################################################################################
    # experiments in each set are numbered 1..10 but in this list they're numbered sequentially
    # They link to sample, sample conditions, and spectrometers, those must exist
    #
    def _make_experiments( self ) :
        if self._verbose : sys.stdout.write( "%s.make_experiments()\n" % (self.__class__.__name__,) )

        assert len( self._samples ) > 0
        assert len( self._conditions ) > 0
        assert len( self._spectrometers ) > 0

        exptnum = 0

        sf = { "sf_category" : "experiment_list", "name" : "experiment_list" }
        sf["experiment"] = []
        sf["experiment_file"] = []

        for i in self._idx["data"].keys() :

            self._experiments[i] = []

            botlog = os.path.realpath( os.path.join( self._indir, self._idx["data"][i]["botlog"] ) )
            bl = topspin.NMRBotLog.parse( filename = botlog, verbose = self._verbose )

            tarball = os.path.realpath( os.path.join( self._indir, self._idx["data"][i]["timedomain"] ) )
            files = topspin.ExperimentFiles.read( tarfile = tarball, verbose = self._verbose )

            if (len( bl.data ) < 1) or (len( files.data ) < 1) :
                raise Exception( "No experiments/experiment directories in set %s" % (i,) )

# need to add this first: 'details' maps back to expt id later on
#
            for j in sorted( bl.data.keys() ) :
                exptnum += 1
                self._experiments[i].append( { "id" : exptnum, "name" : bl.data[j]["name"],
                        "orig_id" : j, "dims" : bl.data[j]["dims"] } )

                for k in files.data.keys() :
                    if str( j ) == str( files.data[k]["experiment_id"] ) :
                        files.data[k]["details"] = bl.data[j]["name"]
                        if files.data[k]["content"] == "NMR experiment directory" :
                            bl.data[j]["raw_data_flag"] = "yes"
                        break

# make saveframe
#  experiments
#
            for j in sorted( bl.data.keys(), cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :

                expid = None
                for k in self._experiments[i] :
                    if k["name"] == bl.data[j]["name"] :
                        expid = k["id"]
                        break
                if expid is None :
                    raise Exception( "Can't find experiment id for %s:%s" % (i,bl.data[j]["name"]) )

                expt = { "id" : expid, "name" : bl.data[j]["name"], "sample_label" : self._samples[i],
                        "sample_condition_list_label" : self._conditions[i],
                        "nmr_spectrometer_label" : self._spectrometers[i] }
                if "raw_data_flag" in bl.data[j].keys() :
                    expt["raw_data_flag"] = "yes"
                else :
                    expt["raw_data_flag"] = "no"
                sf["experiment"].append( expt )

#  experiement files
#

            for j in sorted( files.data.keys(), cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :
                expid = None
                for k in self._experiments[i] :
                    if k["name"] == files.data[j]["details"] :
                        expid = k["id"]
                        break
                if expid is None :
                    raise Exception( "Can't find experiment id for file %s:%s" % (i,j) )

                exptfile = files.data[j]
                exptfile["experiment_id"] = expid
                exptfile["dataset_id"] = i
                sf["experiment_file"].append( exptfile )

# peak list
# NOTE type text/xml is hardcoded
#
            for (j,v) in sorted( self._idx["data"][i]["PEAKFILES"].items(),
                        cmp = lambda a, b : cmp( a[1], b[1] ) ) :
                expid = None
                for k in self._experiments[i] :
                    if k["name"] == v :
                        expid = k["id"]
                        break
                if expid is None :
                    raise Exception( "Can't find experiment id for peak file %s:%s" % (i,j) )

                exptfile = { "experiment_id" : expid, "name" : j, "type" : "text/xml",
                            "details" : "TopSpin XML", "content" : "Peak list",
                            "dataset_id" : i }
                sf["experiment_file"].append( exptfile )

# pictures
# NOTE type image/png is hardcoded
#
            for (j,v) in sorted( self._idx["data"][i]["PNGFILES"].items(),
                        cmp = lambda a, b : cmp( a[1], b[1] ) ) :
                expid = None
                for k in self._experiments[i] :
                    if k["name"] == v :
                        expid = k["id"]
                        break
                if expid is None :
                    raise Exception( "Can't find experiment id for peak file %s:%s" % (i,j) )

                exptfile = { "experiment_id" : expid, "name" : j, "type" : "image/png",
                            "content" : "Spectral image", "dataset_id" : i }
                sf["experiment_file"].append( exptfile )


# there's only one but we return the list of one for consistency
#
        self._dat.append( [sf] )

    ################################################################################################
    # read assignments from csv, add experiments
    #  parse csv returns a single dict
    #
    def _make_cs_assignments( self ) :
        if self._verbose : sys.stdout.write( "%s.make_cs_assignments()\n" % (self.__class__.__name__,) )

# assignments need to be sorted by chem_comp_atom.ordinal
# (so this has to come after make_chem_comp)
# NOTE: there can be only one chem comp
#
        chem_comp = None
        for i in self._dat :
            for j in i :
                if j["sf_category"] == "chem_comp" :
                    chem_comp = j
        if chem_comp is None :
            raise Exception( "No chem comp" )
        if not "atoms" in chem_comp.keys() :
            raise Exception( "No atoms in chem comp" )

        data = []

        for i in self._idx["data"].keys() :
            csv = os.path.realpath( os.path.join( self._indir, self._idx["data"][i]["shifts"] ) )
            try :
                cs = chemshifts.Resonance.readfile( filename = csv, verbose = self._verbose )
            except :
                sys.stderr.write( "Parsing %s:\n" % (csv,) )
                raise

# make saveframe
#
            sfname = "chem_shift_%s" % ((nmrstar.DATASET_MASK % (int( i ),)),)
            self._shifts[i] = sfname
            sf = { "name" : sfname }
            sf["sf_category"] = "assigned_chemical_shifts"
            sf["chem_shift_reference_label"] = self._csrs[i]
            sf["sample_condition_list_label"] = self._conditions[i]

# assignemnts
# check
#
            for s in cs.shifts_table :
                found = False
                for a in chem_comp["atoms"] :
                    if a["id"] == s["atom_id"] :
                        found = True
                        break
                if not found :
                    raise Exception( "atom %s is not found in chem. comp" % (s["atom_id"],) )
# sort
#
            sf["atom_chem_shift"] = []
            num = 0
            ambi = False
            for a in sorted( chem_comp["atoms"], cmp = lambda x, y : cmp( int( x["ordinal"] ), int( y["ordinal"] ) ) ) :
                for s in cs.shifts_table :
                    if s["atom_id"] == a["id"] :
                        num += 1
                        sf["atom_chem_shift"].append( { "id" : num, 
                                                "atom_id" : s["atom_id"],
                                                "ambiguity_code" : s["ambiguity_code"],
                                                "resonance_id" : s["resonance_id"],
                                                "value" : s["value"] } )
                        if s["ambiguity_code"] > 3 :
                            ambi = True

# this table is redundant becasue we have rsonance ids, but whatever
#
            if ambi :
                ambilist = {}
                for s in sf["atom_chem_shift"] :
                    if s["ambiguity_code"] > 3 :
                        if s["resonance_id"] in ambilist.keys() :
                            ambilist[s["resonance_id"]].append( s["id"] )
                        else :
                            ambilist[s["resonance_id"]] = [s["id"]]

                num = 0
                sf["ambiguous_atom_chem_shift"] = []
                for j in sorted( ambilist.keys(), cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :
                    num += 1
                    for k in sorted( ambilist[j], cmp = lambda a, b : cmp( int( a ), int( b ) ) ) :
                        sf["ambiguous_atom_chem_shift"].append( { "set_id" : num, "shift_id" : k } )

# experiments
#
            sf["chem_shift_experiment"] = []

            for j in sorted( self._experiments[i], cmp = lambda a, b : cmp( int( a["id"] ), int( b["id"] ) ) ) :
                sf["chem_shift_experiment"].append( { "experiment_id" : j["id"], "experiment_name" : j["name"],
                                                    "sample_label" : self._samples[i] } )

            data.append( sf )

        self._dat.append( data )

    ################################################################################################
    # build list of peak saveframe dicts
    #
    def _make_peak_lists( self ) :
        if self._verbose : sys.stdout.write( "%s.make_peak_lists()\n" % (self.__class__.__name__,) )

        rc = []
        for dataset_num in self._idx["data"].keys() :
            for peakfile in self._idx["data"][dataset_num]["PEAKFILES"].keys() :
                fname = os.path.realpath( os.path.join( self._indir, peakfile ) )
                exptname  = self._idx["data"][dataset_num]["PEAKFILES"][peakfile]
                sf = { "sf_category" : "spectral_peak_list" }
                sf["name"] = "spectral_peaks_%s_set%02d" % (re.sub( r"\s+", "_", exptname ),int( dataset_num ),)
                sf["experiment_name"] = exptname
                sf["sample_label"] = self._samples[dataset_num]
                sf["sample_condition_list_label"] = self._conditions[dataset_num]
                sf["assigned_chem_shift_list_label"] = self._shifts[dataset_num]

                dim = []
                cnt = 0
                dtl = ""

# spectral_dim and top-level detail is frem experiments
#
                for experiment in self._experiments[dataset_num] :
                    if experiment["name"] == exptname :
                        sf["experiment_id"] = experiment["id"]
                        sf["experiment_name"] = exptname
                        for dim_num in sorted( experiment["dims"].keys() ) :
                            m = re.search( r"(\d+)([A-Za-z]+)", experiment["dims"][dim_num]["nuc"] )
                            if not m :
                                raise Exception( "Bad nucleus %s in %s" % \
                                    (experiment["dims"][dim_num]["nuc"],fname,) )
                            cnt += 1

# there's nowhere to put topspin's dimension labels: F1, F2, etc., so we put them in details at top level
#
                            if dtl == "" : dtl = "%s: %s" % (dim_num, experiment["dims"][dim_num]["nuc"])
                            else : dtl = "%s, %s: %s" % (dtl, dim_num, experiment["dims"][dim_num]["nuc"])

                            sd = None
                            if m.group( 2 ) == "H" : sd = "Full H"
                            elif m.group( 2 ) == "C" : sd = "Full C"
                            else :
                                raise Exception( "Unsupported atom type %s in %s" % (m.group( 2 ),fname,) )

# FIXME ppm, 2 dimensions, 'relaitive height' are hardcoded
#
                            dim.append( { "id" : cnt, "atom_type" : m.group( 2 ), "atom_isotope_number" : m.group( 1 ),
                                "spectral_region" : sd, "sweep_width" : experiment["dims"][dim_num]["sw"],
                                "sweep_width_units" : "ppm" } )

                        break

                if len( dim ) < 1 :
                    raise Exception( "No spectral dimensions in %s" % (fname,) )

                sf["details"] = dtl
                sf["number_of_spectral_dimensions"] = cnt
                sf["spectral_dim"] = dim

# the peaks
#
                trans = []
                trans_gen = []
                trans_char = []
                trans_ass = []

                peaks = topspin.TopSpinPeakList.parse( filename = fname, verbose = self._verbose )
                if peaks.data is None :
                    sys.stderr.write( "ERR: no transitions in %s\n" % (fname,) )
                    sys.stderr.write( "NOTE that this may be caused by empty extra peaklist element(s): topspin does that sometimes\n" )
                    raise Exception( "No transitions in %s" % (fname,) )

#                if str( peaks.data["experiment"] ).strip() != str( sf["experiment_id"] ).strip() :
                if str( peaks.data["experiment"] ).strip() != str( experiment["orig_id"] ).strip() :

                    sys.stderr.write( "ERR: in file %s\n" % (fname,) )
                    sys.stderr.write( "**** for expt name %s\n" % (exptname,) )
                    sys.stderr.write( "**** peaks ****\n" )
                    pprint.pprint( peaks.data, stream = sys.stderr )
                    sys.stderr.write( "**** experiments ****\n" )
                    pprint.pprint( self._experiments[dataset_num], stream = sys.stderr )
                    raise Exception( "Expt ID in peaks %s (%s) does not match id from index %s" \
                            % (fname,peaks.data["experiment"],experiment["orig_id"],) )

                for peak_num in peaks.data["peaks"].keys() :
                    trans.append( { "id" : str( peak_num ).strip() } )
                    trans_gen.append( { "spectral_transition_id" : str( peak_num ).strip(),
                    "intensity_val" : peaks.data["peaks"][peak_num]["intensity"],
                    "measurement_method" : "relative height" } )

                    trans_char.append( { "spectral_transition_id" : str( peak_num ).strip(),
                        "spectral_dim_id" : 1,
                        "chem_shift_val" : peaks.data["peaks"][peak_num]["shift1"] } )
# 2D?
#
                    if "shift2" in peaks.data["peaks"][peak_num].keys() :
                        trans_char.append( { "spectral_transition_id" : str( peak_num ).strip(),
                            "spectral_dim_id" : 2,
                            "chem_shift_val" : peaks.data["peaks"][peak_num]["shift2"] } )

# has atom label?
# there is never an atom label on 2D peaks, as of the time of this writing
#
                    if "atom" in peaks.data["peaks"][peak_num].keys() :
                        trans_ass.append( { "spectral_transition_id" : str( peak_num ).strip(),
                            "spectral_dim_id" : 1,
                            "val" : peaks.data["peaks"][peak_num]["shift1"],
                            "atom_id" : peaks.data["peaks"][peak_num]["atom"] } )

                if (len( trans ) < 1) or (len( trans_gen ) < 1) or (len( trans_char ) < 1) :
                    raise Exception( "No transitions in %s" % (filename,) )

                sf["spectral_transition"] = trans
                sf["spectral_transition_general_char"] = trans_gen
                sf["spectral_transition_char"] = trans_char
                if len( trans_ass ) > 0 :
                    sf["assigned_spectral_transition"] = trans_ass

                rc.append( sf )

        self._dat.append( rc )

    ################################################################################################

#
#
#
if __name__ == '__main__':

    import boilerplate
    import chemcomp
    import chemshifts
    import nmrfam
    import nmrstar
    import sample
    import topspin
    import www

    dat = FAMtoJSN.make_entry( indexfile = sys.argv[1], verbose = False )
#    dat.add_alatis_map()
    sys.stdout.write( json.dumps( dat.data, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )


#
#
