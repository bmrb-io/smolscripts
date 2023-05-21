#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  molecule.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit


from __future__ import absolute_import

import sys
import os
import json
import pprint

_HERE = os.path.split( __file__ )[0]
sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
import bmrbmb.chemcomp

# combine stuff extracted by OpenBabel and RDKit into a (partial) chem comp
#
class Molecule( object ) :

    # main
    #
    @classmethod
    def from_file( cls, filename, format = "sdf", verbose = False ) :
        if verbose : sys.stdout.write( "%s.from_file(%s)\n" % (cls.__name__,filename,) )
        infile = os.path.realpath( filename )

# if fucking python decides to convert the filename string to unicode somewhere along the way,
# c++ bindings in the toolkits below will fail
#
# since unicode handling is different between 2 & 3, code below probably makes it v.2 only.
# until someone fixes this shit.
#

#        if isinstance( infile, unicode ) :
#            infile = infile.encode( "ascii" )
        obm = bmrbmb.chemcomp.OBmolecule.from_file( filename = infile, format = format, verbose = verbose )
        rdm = bmrbmb.chemcomp.RDmolecule.from_file( filename = infile, format = format, verbose = verbose )

        mol = cls( obmol = obm, rdmol = rdm, verbose = verbose )

        return mol

    #
    #
    @classmethod
    def from_string( cls, string, format = "sdf", verbose = False ) :
        if verbose : sys.stdout.write( "%s.from_string()\n" % (cls.__name__,) )
        if not format == "sdf" :
            raise NotImplementedError( "from_string(): %s" % (format,) )
        obm = bmrbmb.chemcomp.OBmolecule.from_string( string = string, format = format, verbose = verbose )
        rdm = bmrbmb.chemcomp.RDmolecule.from_string( string = string, format = format, verbose = verbose )

        mol = cls( obmol = obm, rdmol = rdm, verbose = verbose )

        return mol


    #
    #
    def __init__( self, obmol, rdmol, verbose = False ) :

        if verbose : sys.stdout.write( "%s.__init__()\n" % (self.__class__.__name__,) )

# this requires importing a whole lot of stuff, and resp. obj constructors do this anyway
#
#        assert isinstance( obmol, pybel.Molecule )
#        assert isinstance( rdmol, rdkit.Chem.rdchem.Mol )

        self._obmol = obmol
        self._rdmol = rdmol
        self._verbose = bool( verbose )

    # wrapper for RDKit MOL writer
    #
    def __mol__( self ) :
        if self._verbose : sys.stdout.write( "%s.__mol__()\n" % (self.__class__.__name__,) )
        return self._rdmol.mol
    mol = property( __mol__ )

    # merge RDKit and OpenBabel output
    #
    def __smiles__( self ) :
        if self._verbose : sys.stdout.write( "%s.__smiles__()\n" % (self.__class__.__name__,) )
        rc = self._rdmol.smiles
        rc.extend( self._obmol.smiles )
        return rc
    smiles = property( __smiles__ )

    # wrapper for OB
    #
    def __inchi__( self ) :
        return self._obmol.inchi
    inchi = property( __inchi__ )

    # chem comp stuff
    # can get big-ish, don't use on large molecules
    #
    @property
    def chem_comp( self ) :

        rc = { "sf_category" : "chem_comp" }

# atoms
#
        rc["atoms"] = list( self.iter_atoms() )

# bonds
#
        rc["bonds"] = list( self.iter_bonds() )

# smiles & inchi go into descriptors
#  merge with upstream ones as well
#
# names go into identifiers, we only have the upstream ones at this point
#
# squash duplicates while we're at it
#
        rc["dblinks"] = []
        identifiers = []
        descriptors = self.smiles
        descriptors.extend( self.inchi )
        systematic_name = None

        for (key, val) in self.iter_sdfdata() :

            k = str( key ).strip().upper()
            v = str( val ).strip()

# if same name has several "types", prefer "systematic" kind
#
            if k.upper() == "PUBCHEM_IUPAC_SYSTEMATIC_NAME" :
                systematic_name = v

            found = False
            for i in descriptors :
                if i["descriptor"] == v :
                    found = True
                    break
            if not found :
                if k == "PUBCHEM_IUPAC_INCHI" :
                    descriptors.append( { "descriptor" : v, "type" : "InChI",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_IUPAC_INCHIKEY" :
                    descriptors.append( { "descriptor" : v, "type" : "InChI_KEY",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_OPENEYE_ISO_SMILES" :
                    descriptors.append( { "descriptor" : v, "type" : "SMILES_ISOMERIC",
                            "program" : "PUBCHEM_OPENEYE", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_OPENEYE_CAN_SMILES" :
                    descriptors.append( { "descriptor" : v, "type" : "SMILES_CANONICAL",
                            "program" : "PUBCHEM_OPENEYE", "version" : "na" } )
                    continue

# still here?
#
            found = False
            for i in identifiers :
                if i["identifier"] == v :
                    found = True
                    break
            if not found :
                if k == "PUBCHEM_IUPAC_SYSTEMATIC_NAME" :
                    identifiers.append( { "identifier" : v, "type" : "SYSTEMATIC NAME",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_IUPAC_TRADITIONAL_NAME" :
                    identifiers.append( { "identifier" : v, "type" : "TRADITIONAL NAME",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_IUPAC_NAME" :
                    identifiers.append( { "identifier" : v, "type" : "NAME",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_IUPAC_OPENEYE_NAME" :
                    identifiers.append( { "identifier" : v, "type" : "OPENEYE NAME",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue
                elif k == "PUBCHEM_IUPAC_CAS_NAME" :
                    identifiers.append( { "identifier" : v, "type" : "CAS NAME",
                            "program" : "PUBCHEM_IUPAC", "version" : "na" } )
                    continue

# still here?
#
            if k == "PUBCHEM_COMPOUND_CID" :
                rc["dblinks"].append( { "db_code" : "PubChem", "acc_type" : "cid", "acc_code" : v } )

# alatis
#
            if k == "ALATIS_STANDARD_INCHI" :
                descriptors.append( { "descriptor" : v, "type" : "InChI", "program" : "ALATIS", "version" : "na" } )

# endfor

        rc["descriptors"] = descriptors

        if systematic_name is not None :
            for i in range( len( identifiers ) ) :
                if identifiers[i]["identifier"] == systematic_name :
                    identifiers[i]["type"] = "SYSTEMATIC NAME"

        rc["identifiers"] = identifiers

# misc. stuff
#
        rc["formula"] = self.formula
        rc["charge"] = self.charge

# don't recalculate for every call
#
        cnts = self.atom_counts
        rc["num_atoms_all"] = cnts["all"]
        rc["num_atoms_nh"] = cnts["non-H"]

        wts = self.masses
        rc["weight"] = wts["mass"]
        rc["weight_monoisotopic_nat"] = wts["monoisotopic"]
        rc["weight_monoisotopic_c13"] = wts["c13"]
        rc["weight_monoisotopic_n15"] = wts["n15"]
        rc["weight_monoisotopic_c13n15"] = wts["c13n15"]

        return rc

    #  formatted as json
    #
    def __json__( self ) :
        return json.dumps( self.chem_comp )

    json = property( __json__ )

    # formula
    #
    def __formula__( self ) :
        return self._obmol._mol.formula
    formula = property( __formula__ )

    # charge
    #
    def __charge__( self ) :
        return self._obmol._mol.charge
    charge = property( __charge__ )

    # masses
    # we could fetch the table from NIST and calculate weights ourselves...
    #
    def __masses__( self ) :
        rc = { "mass" : self._obmol._mol.molwt, "monoisotopic" : self._obmol._mol.exactmass }
        rc["c13"] = self._obmol.monoisotopic_mass( carbon = 13 )
        rc["n15"] = self._obmol.monoisotopic_mass( nitrogen = 15 )
        rc["c13n15"] = self._obmol.monoisotopic_mass( carbon = 13, nitrogen = 15 )
        return rc
    masses = property( __masses__ )

    # atoms
    # NMR-STAR chem comp has fields for num atoms and num non-hydrogens
    #
    def __atom_counts__( self ) :
        rc = { "all" : 0, "non-H" : 0 }
        for atom in self.iter_atoms() :
            rc["all"] += 1
            if not atom["type"] in ("H","D","T") :
                rc["non-H"] += 1
        return rc
    atom_counts = property( __atom_counts__ )

    #
    #
    def iter_atoms( self ) :
        """atoms iterator, wrapper for RDmolecule.iter_atoms"""
        return self._rdmol.iter_atoms()

    # this one combines bond types from OpenBabel with the rest of bond info from RDKit.
    # bonds are matched on atom names which should be the same as long as nobody messed with the
    # OB/RD mols before calling this.
    #
    # return { "id" : N, "ord" : X, "arom" : yes if "ord" is "AROM", "stereo" : N|E|Z,
    #        "type" : covalent|amide|ester|carbonyl, "atom1" : I, "atom2" : J }
    #
    def iter_bonds( self ) :
        """bonds iterator"""

# fewer columns in OB
#
        obbonds = {}
        for bond in self._obmol.iter_bonds() :
            obbonds["%s<->%s" % (bond["atom1"],bond["atom2"],)] = bond["type"]

        for bond in self._rdmol.iter_bonds() :
            key = "%s<->%s" % (bond["atom1"],bond["atom2"],)
            if not key in obbonds.keys() :
                raise Exception( "No such bond in OpenBabel molecule: %s" % (key,) )

            bond["type"] = obbonds[key]
            yield bond

    # extra stuff in SDF input file
    #
    def iter_sdfdata( self ) :
        """wrapper for RDmolecule.iter_sdfdata"""
        return self._rdmol.iter_sdfdata()

    # alatis atom map: only needed if atom names in the assignments don't match alatis ordering
    # "atom_id" is from the input to alatis, "alatis_id" is the new name
    #
    def alatis_map( self ) :

        rc = []
        for (key, val) in self.iter_sdfdata() :

            if str( key ).strip().upper() == "ALATIS_MAP" :

# Hesam's output is fucked up. hope he does not change it
#
                first = True
                atmidx = 1
                newidx = 0
                oldidx = 2
                for line in val.splitlines() :
                    if first :
#                        vals = line.split( "\t" )
#                        for i in range( len (vals ) ) :
#                            if str( vals[i] ).strip().upper().startswith( "NEW" ) :
#                                newidx = i
#                            elif str( vals[i] ).strip().upper().startswith( "ATOM" ) :
#                                atmidx = i
#                            elif str( vals[i] ).strip().upper().startswith( "ORIG" ) :
#                                oldidx = i
#                            else :
#                                raise Exception( "Unknown column header in ALATIS_MAP" )
                        first = False
                        continue

                    vals = line.split()

                    rc.append( { "atom_id" : "%s%d" % (str( vals[atmidx] ).strip(), int( vals[oldidx] ) ),
                            "alatis_id" : "%s%d" % (str( vals[atmidx] ).strip(), int( vals[newidx] ) ) } )

        if len( rc ) < 1 : return None
        return rc

########################################################################################
#
#
#
if __name__ == '__main__':

    m = Molecule.from_file( sys.argv[1], verbose = True )
    m.add_alatis_map()
    pretty = json.loads( m.json )
    sys.stdout.write( json.dumps( pretty, indent = 4, sort_keys = True ) )
    sys.stdout.write( "\n" )

#
#
