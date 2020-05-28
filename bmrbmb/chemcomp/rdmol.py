#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  rdmol.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit


from __future__ import absolute_import

import os
import sys
import json
import re

import pprint

import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.Draw
import rdkit.Chem.inchi
import rdkit.Chem.rdchem
import rdkit.Chem.rdMolDescriptors
#bugbubug
from rdkit import rdBase

# RDKit 2018-03-1+
#import rdkit.Chem.rdCoordGen


# RDKit, unlike OpenBabel, keeps proper stereochemistry flags. So we use it to build atom and bond
# tables.
#
class RDmolecule( object ) :

    BOND_TYPES = {
        rdkit.Chem.rdchem.BondType.AROMATIC      : "AROM",
        rdkit.Chem.rdchem.BondType.DATIVE        : "DATIVE",
        rdkit.Chem.rdchem.BondType.DATIVEL       : "DATIVEL",
        rdkit.Chem.rdchem.BondType.DATIVEONE     : "DATIVEONE",
        rdkit.Chem.rdchem.BondType.DATIVER       : "DATIVER",
        rdkit.Chem.rdchem.BondType.DOUBLE        : "DOUB",
        rdkit.Chem.rdchem.BondType.FIVEANDAHALF  : "FIVEANDAHALF",
        rdkit.Chem.rdchem.BondType.FOURANDAHALF  : "FOURANDAHALF",
        rdkit.Chem.rdchem.BondType.HEXTUPLE      : "HEXTUPLE",
        rdkit.Chem.rdchem.BondType.HYDROGEN      : "HYDROGEN",
        rdkit.Chem.rdchem.BondType.IONIC         : "IONIC",
        rdkit.Chem.rdchem.BondType.ONEANDAHALF   : "ONEANDAHALF",
        rdkit.Chem.rdchem.BondType.OTHER         : "OTHER",
        rdkit.Chem.rdchem.BondType.QUADRUPLE     : "QUAD",
        rdkit.Chem.rdchem.BondType.QUINTUPLE     : "QUINTUPLE",
        rdkit.Chem.rdchem.BondType.SINGLE        : "SING",
        rdkit.Chem.rdchem.BondType.THREEANDAHALF : "THREEANDAHALF",
        rdkit.Chem.rdchem.BondType.THREECENTER   : "THREECENTER",
        rdkit.Chem.rdchem.BondType.TRIPLE        : "TRIP",
        rdkit.Chem.rdchem.BondType.TWOANDAHALF   : "TWOANDAHALF",
        rdkit.Chem.rdchem.BondType.UNSPECIFIED   : "UNSPECIFIED",
        rdkit.Chem.rdchem.BondType.ZERO          : "ZERO"
    }

    BOND_STEREO = {
        rdkit.Chem.rdchem.BondStereo.STEREOE : "E",
        rdkit.Chem.rdchem.BondStereo.STEREOZ : "Z",
        rdkit.Chem.rdchem.BondStereo.STEREONONE : "N",
        rdkit.Chem.rdchem.BondStereo.STEREOANY  : "UNK"
    }

# rdkit returns 'R', 'S', and 'N' for these in python so they're not needed anymore
#
#    ATOM_STEREO = (
#        rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
#        rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
#        rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
#    )



    # RDKit removed hydrogens by default. We're supposed to have them and in proper order already.
    #
    @classmethod
    def from_string( cls, string, format = "sdf", verbose = False ) :
        if verbose : sys.stdout.write( "%s.from_string()\n" % (cls.__name__,) )
        if format != "sdf" :
            raise NotImplementedError( "from_string(): %s" % (format,) )
        mol = rdkit.Chem.MolFromMolBlock( string, removeHs = False )
        return cls( mol, verbose )

    #
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

        if isinstance( infile, unicode ) :
            infile = infile.encode( "ascii" )

        if sys.version_info[0] == 2 :
            mol = rdkit.Chem.SupplierFromFilename( infile, removeHs = False ).next()
        elif sys.version_info[0] == 3 :
            mol = rdkit.Chem.SupplierFromFilename( infile, removeHs = False ).__next__()
        else : raise Exception( "golang is looking better every day" )
        return cls( mol, verbose )

    #
    #
    def __init__( self, molecule, verbose = False ) :
        assert isinstance( molecule, rdkit.Chem.rdchem.Mol )
        self._mol = molecule
        self._verbose = verbose

    #
    #
    def __mol__( self ) :
        if self._verbose : sys.stdout.write( "%s.__mol__()\n" % (self.__class__.__name__,) )
        return rdkit.Chem.MolToMolBlock( self._mol )
    mol = property( __mol__ )

    #
    #
    @property
    def version( self ) :
        return str( rdBase.rdkitVersion ).strip()

    #
    #
    @property
    def program( self ) :
        return "RDKit"

    # RDKit will output all H'es (since the molecule has them), that's not typically done in smiles.
    # so we remove them here
    #
    def to_smiles( self, kind = None ) :

        if self._verbose : sys.stdout.write( "%s.to_smiles(%s)\n" % (self.__class__.__name__,kind,) )
        rc = None
        m = rdkit.Chem.RemoveHs( self._mol )
        if kind == "isomeric" :
            rc = rdkit.Chem.MolToSmiles( m, isomericSmiles = True )
        elif kind == "canonical" :
            rc = rdkit.Chem.MolToSmiles( m, canonical = True )
        else :
            rc = rdkit.Chem.MolToSmiles( m, isomericSmiles = False, canonical = False )

        if rc is not None :
            rc = str( rc ).strip()
            if len( rc ) < 1 :
                rc = None

        return rc

    def __smiles__( self ) :
        rc = []
        can = self.to_smiles( kind = "canonical" )
        if can is not None :
            rc.append( { "descriptor" : can, "type" : "SMILES_CANONICAL", "program" : self.program,
                        "version" : self.version } )
        iso = self.to_smiles( kind = "isomeric" )
        if iso is not None :
            rc.append( { "descriptor" : can, "type" : "SMILES_ISOMERIC", "program" : self.program,
                        "version" : self.version } )
        smi = self.to_smiles()
        if smi is not None :
            rc.append( { "descriptor" : can, "type" : "SMILES", "program" : self.program,
                        "version" : self.version } )
        if len( rc ) < 1 :
            rc = None
        return rc

    smiles = property( __smiles__ )

    # chem comp stuff formatted as json, or none
    # can be big-ish, don't use on large molecules
    #
    def __json__( self ) :

        rc = {}

        atoms = list( self.iter_atoms() )
        if len( atoms ) > 0 :
            rc["chem_comp_atom"] = atoms

        bonds = list( self.iter_bonds() )
        if len( bonds ) > 0 :
            rc["chem_comp_bond"] = bonds

        sm = self.smiles
        if sm is not None :
            rc["chem_comp_descriptor"] = sm

        extras = list( self.iter_sdfdata() )
        if len( extras ) > 0 :
            rc["extra_properties"] = extras

        return json.dumps( rc )

    json = property( __json__ )

    #
    # pictures
    #
    def to_svg( self, title = None, full = True ) :
        if self._verbose : sys.stdout.write( "%s.to_svg(%s)\n" % (self.__class__.__name__,(full and "full" or "small")) )

        if title is None : t = ""
        if str( title ).strip()  == "" : t = ""
        if isinstance( title, unicode ) :
            t = title.encode( "ascii" )

        if full :
            rc = self._to_svg_full( title = t )
        else :
            rc = self._to_svg_small( title = t )

        rectpat = re.compile( r'<svg:rect\s.+</svg:rect>' )
        svg = re.sub( rectpat, '', rc )

        return svg

#
# RDKit 2018 has a new 2D geometry optimizer but it won't build on centos 7.
#
    #
    #
    def _to_svg_small( self, title = None ) :
        """return 'presentation' SVG"""
        if self._verbose : sys.stdout.write( "%s._to_svg_small()\n" % (self.__class__.__name__,) )

        dr = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG( 1000, 1000 )
        dr.SetFontSize( 0.7 ) # 0.27 )
#        op = dr.drawOptions()
#        for i in range( self._mol.GetNumAtoms() ) :
#            op.atomLabels[i] = self._mol.GetAtomWithIdx( i ).GetSymbol() + str( (i + 1) )

        m = rdkit.Chem.RemoveHs( self._mol )
        rdkit.Chem.AllChem.Compute2DCoords( m )

        dr.DrawMolecule( m )
        dr.FinishDrawing()
        return dr.GetDrawingText()

    # 'full' pix has all atoms with labels.
    #
    def _to_svg_full( self, title = None, newstyle = False, minimize = None ) :
        """return 'full' SVG"""
        if self._verbose : sys.stdout.write( "%s._to_svg_full()\n" % (self.__class__.__name__,) )


        dr = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG( 1000, 1000 )
        dr.SetFontSize( 0.7 ) # 0.27 )
        op = dr.drawOptions()
        op.clearBackground = False
        for i in range( mol.GetNumAtoms() ) :
            op.atomLabels[i] = mol.GetAtomWithIdx( i ).GetSymbol() + str( (i + 1) )

        rdkit.Chem.AllChem.Compute2DCoords( mol )

# rdkit 2018-03-1+
#
        if newstyle :
            rdkit.Chem.rdCoordGen.AddCoords( mol )
# self._mol.GetConformer( 0 ).GetPositions()[0] has generated coords

# sometimes you get better results if you do this:
#
        if minimize is not None :
            self.minimize( field = minimize )

        dr.DrawMolecule( mol )
        dr.FinishDrawing()
        return dr.GetDrawingText()

    # RDKit's implementation is ugly: there's separate call for every force field...
    #
    def minimize( self, field = "UFF" ) :
        """minimize energy using specified force field"""
        if self._verbose : sys.stdout.write( "%s.minimize(%s)\n" % (self.__class__.__name__,field,) )

        if not str( field ).upper() in ("UFF","MMFF94","MMFF94S") :
            return

        rc = rdkit.Chem.AllChem.EmbedMolecule( self._mol )
        if rc < 0 :
            rc = rdkit.Chem.AllChem.EmbedMolecule( self._mol, useRandomCoords = True )
        try :

# if at first it doesn't converge, try try try again
#
            if str( field ).upper() == "UFF" :
                if rdkit.Chem.AllChem.UFFOptimizeMolecule( self._mol ) == 1 :
                    rdkit.Chem.AllChem.UFFOptimizeMolecule( self._mol, maxIters = 1000 )

            if str( field ).upper() in ("MMFF94", "MMFF94S") :
                if rdkit.Chem.AllChem.MMFFOptimizeMolecule( self._mol, mmffVariant = field ) == 1 :
                    rdkit.Chem.AllChem.MMFFOptimizeMolecule( self._mol, mmffVariant = field, maxIters = 1000 )

        except ValueError :
            pass

###############################
    #
    # chem comp atom
    # return { "id" : XN, "type" : X, "ordinal" : N, "stereo" : N|R|S }
    # "aromatic flag" is set from bonds, so that's done somewhere else
    #
    def iter_atoms( self ) :
        """atoms iterator"""
        if self._verbose : sys.stdout.write( "%s.iter_atoms()\n" % (self.__class__.__name__,) )

# coming from ALATIS, we should have 3D MOL with all stereochemistry information
#
        stereos = rdkit.Chem.FindMolChiralCenters( self._mol, includeUnassigned = False )
        for atom in self._mol.GetAtoms() :
#            print "****** atom ", atom.GetIdx()
            rc = { "ordinal" : (atom.GetIdx() + 1), "type" : atom.GetSymbol() }
            rc["id"] = "%s%s" % (rc["type"], rc["ordinal"],)
            rc["stereo"] = 'N'
            for (idx, flag,) in stereos :
                if idx == (rc["ordinal"] - 1) :
                    if (flag is not None) and (str( flag ).strip() != "?") :
                        rc["stereo"] = str( flag ).strip()
                        break
            rc["charge"] = atom.GetFormalCharge()
            rc["aromatic"] = atom.GetIsAromatic()

            yield rc

    # chem comp bond
    # return { "id" : N, "ord" : X, "arom" : yes if "ord" is "AROM", "stereo" : N|E|Z,
    #        "atom1" : I, "atom2" : J }
    #
    # RDKit does not provide bond type: carbonyl, ester, etc. that I know of, that comes from OpenBabel
    #
    def iter_bonds( self ) :
        """bonds iterator"""
        if self._verbose : sys.stdout.write( "%s.iter_bonds()\n" % (self.__class__.__name__,) )

        for bond in self._mol.GetBonds() :
            rc = { "id" : (bond.GetIdx() + 1) }
            rc["arom"] = "no"
            rc["ord"] = self.BOND_TYPES[bond.GetBondType()]
            if rc["ord"] == "AROM" : rc["arom"] = "yes"
            rc["stereo"] = self.BOND_STEREO[bond.GetStereo()]
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            rc["atom1"] = "%s%d" % (atom1.GetSymbol(), (atom1.GetIdx() + 1))
            rc["atom2"] = "%s%d" % (atom2.GetSymbol(), (atom2.GetIdx() + 1))

            yield rc

    # extra stuff in SDF
    #
    def iter_sdfdata( self ) :
        """iterate over extra properties extracted from input (SDF)"""

        return self._mol.GetPropsAsDict().iteritems()

#
#
#
if __name__ == '__main__':
    m = RDmolecule.from_file( sys.argv[1], verbose = False )
#    print m._to_svg_small()
#    print m._to_svg_full()
#    print m.smiles

#    m.minimize( field = "MMFF94" )
#    print m.mol
    sys.stdout.write( m.json )


#
#
