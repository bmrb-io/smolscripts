#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  schwalbe2json.py
#
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

if __name__ != "__main__" :

#    from . import boilerplate
    from . import chemcomp
    from . import chemshifts
    from . import nmrfam
#    from . import nmrstar
#    from . import sample
#    from . import topspin
    from . import www

"""
Read SDF file and make a BMRB chem comp in JSON format.
"""

#
#
#
class SDFtoJSN( object ) :

    @classmethod
    def run( cls, filename, verbose = False ) :
        cnv = cls( filename, verbose )
        cnv._make_chem_comp()
#        sys.stdout.write( cnv.__json__ )
        return cnv

    #
    def __init__( self, filename, verbose = False ) :
        self._fname = os.path.realpath( filename )
        self._verbose = bool( verbose )
        self._dat = []
        self._mol = None
        self._chem_comp = None
        self._inchi = None

        self._concentration = None
        self._volume = None

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
    # num is always 1 here
    #
    def _make_chem_comp( self, num = 1 ) :
        if self._verbose : sys.stdout.write( "%s.make_chem_comp()\n" % (self.__class__.__name__,) )

        self._mol = chemcomp.Molecule.from_file( filename = self._fname, verbose = self._verbose )
        self._chem_comp = self._mol.chem_comp

# the above will parse out alatis inchi
#
        if "descriptors" in self._chem_comp.keys() :
            for i in self._chem_comp["descriptors"] :
                if (i["type"] == "InChI") and (i["program"] == "ALATIS") :
                    self._inchi = i["descriptor"]
        self._chem_comp["inchi_code"] = self._inchi

        self._chem_comp["paramagnetic"] = "FIXME"

# SDF meta
#
        for (key, val) in self._mol.iter_sdfdata() :
            if str( key ).strip().upper() == "MOLECULE_NAME" :
                if (val is None) or (str( val ).strip() == "" ) :
                    self._chem_comp["name"] = "chem_comp_%s" % (num,)
                else : self._chem_comp["name"] = str( val ).strip()
            if str( key ).strip().upper() == "VENDOR" :
                if (val is None) or (str( val ).strip() == "" ) :
                    pass
                else : self._chem_comp["vendor"] = str( val ).strip()
            if str( key ).strip().upper() == "COMPOUND_ID" :
                if (val is None) or (str( val ).strip() == "" ) :
                    pass
                else : self._chem_comp["vendor_product_code"] = str( val ).strip()
            m = re.search( r"^concentration_(.+)$", str( key ).strip(), re.I )
            if m :
                if (val is None) or (str( val ).strip() == "" ) :
                    pass
                else :
                    self._concentration = { "units" : m.group( 1 ), "value" : str( val ).strip() }
            m = re.search( r"^volume_(.+)$", str( key ).strip(), re.I )
            if m :
                if (val is None) or (str( val ).strip() == "" ) :
                    pass
                else :
                    self._volume = { "units" : m.group( 1 ), "value" : str( val ).strip() }

# web stuff
#
        if not "dblinks" in self._chem_comp.keys() :
            self._chem_comp["dblinks"] = []
        bmrbid = www.bmrbid( self._inchi, self._verbose )
        if not bmrbid is None :
            self._chem_comp["dblinks"].append( bmrbid )

        if not "identifiers" in self._chem_comp.keys() :
            self._chem_comp["identifiers"] = []

#FIXME: need inchi from chem_comp
        iupac_name = www.cactus_iupac_name( self._inchi, self._verbose )
        if not iupac_name is None :
            for name in iupac_name :
                self._chem_comp["identifiers"].append( name )

#        if "cid" in self._idx.keys() :
#            pug = www.get_pubchem_data( self._idx["cid"], self._verbose )
#            if pug is not None :
#                if "identifiers" in pug.keys() :
#                    for i in pug["identifiers"] :
#                        self._chem_comp["identifiers"].append( i )
#                if "common_name" in pug.keys() :
#                    if not "common_name" in self._chem_comp.keys() :
#                        self._chem_comp["common_name"] = []
#                    for i in pug["common_name"] :
#                        self._chem_comp["common_name"].append( i )

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

###############################################################################################
#
# assembly
#
def add_assembly( dat ) :
    comp = dat[0][0]

# fix up if not set:
#
    if not "aromatic" in comp.keys() :
        arom = False 
        for a in comp["atoms"] :
            if a["aromatic"] :
                arom = True
                break
# JIC
#
        if not arom :
            for b in comp["bonds"] :
                if b["arom"] == "yes" :
                    arom = True
                    break

        comp["aromatic"] = arom

# entity: boilerplate
#
    enty = {}
    enty["formula_weight"] = comp["weight"]
    enty["nonpolymer_comp_label"] = comp["name"]
    enty["paramagnetic"] = comp["paramagnetic"]
    enty["id"] = 1
    enty["name"] = "entity_1"
    enty["sf_category"] = "entity" 
    enty["thiol_state"] = "not present"
    enty["type"] = "non-polymer"
    enty["ambiguous_conformational_states"] = "no"
    enty["ambiguous_chem_comp_sites"] = "no"
    enty["nstd_linkage"]  = "no"
    enty["number_of_monomers"] = 1
    enty["number_of_nonpolymer_components"] = 1
    enty["entity_comp_index"] = [ {
        "auth_seq_id": 1,
        "comp_label": comp["name"],
        "entity_id": 1,
        "id": 1
    } ]

    dat.insert( 0, [enty] )

# assembly: copy of atoms
#

    assy = {}
    assy["id"] = 1
    assy["name"] = "assembly_1"
    assy["sf_category"] = "assembly"
    assy["number_of_components"] = 1
    assy["paramagnetic"] = enty["paramagnetic"]
    assy["thiol_state"] = enty["thiol_state"]
    assy["entity_assembly"] = [ {
        "id" : 1, "assembly_id" : 1, "experimental_data_reported" : "yes",
        "physical_state" : "native", "conformational_isomer" : "no",
        "chemical_exchange_state" : "no", "entity_assembly_name" : comp["name"],
        "entity_id" : 1
    } ]

    assy["atoms"] = []
    for a in comp["atoms"] :
        assy["atoms"].append( {
            "assembly_atom_id": a["ordinal"],
            "id": a["id"], 
            "type_symbol": a["type"],
            "assembly_id": 1, 
            "comp_index_id": 1, 
            "entity_assembly_id": 1, 
            "entity_id": 1, 
            "seq_id": 1, 
        } )

    dat.insert( 0, [assy] )


################################################################################################
#
#
#
if __name__ == '__main__':

    import chemcomp
    import www

    cnv = SDFtoJSN.run( sys.argv[1] )
    cnv.add_alatis_map()

    add_assembly( cnv.data )

    if (cnv._volume is not None) or (cnv._concentration is not None) :
        sample = { "name" : "sample", "sf_category" :  "sample" }
        if cnv._volume is not None :
            sample |= { "volume" : cnv._volume }
        if cnv._concentration is not None :
            sample |= { "concentration" : cnv._concentration }

    cnv.data.append( [sample] )

    sys.stdout.write( json.dumps( cnv.data, indent = 4, sort_keys = True ) )


#
#
