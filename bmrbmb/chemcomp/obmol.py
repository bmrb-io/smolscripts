#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  obmol.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit


from __future__ import absolute_import

import os
import sys
import re
import json
from openbabel import openbabel
from openbabel import pybel

# wrapper for OpenBabel
# require:
# - pybel
# - openbabel
#
# NOTE that the input is assumed to come from ALATIS and contain all hydrogens as well as proper
# stereochemistry information. Don't add/remove hydrogens except when making a 'presentation' SVG
#
# OpenBabel tends to generate better depictions with all atoms and atom labels drawn, however for
# some of them RDKit works better. Or even Marvin.
#
class OBmolecule( object ) :
    """OpenBabel molecule"""

    # use at own risk: creating OBMols from formats other than MOL/SDF can result in anything
    #
    @classmethod
    def from_string( cls, string, format = "sdf", verbose = False ) :
        if verbose : sys.stdout.write( "%s.from_string()\n" % (cls.__name__,) )
        mol = pybel.readstring( format, string )
        return cls( mol, verbose )

    # retruns only the 1st molecule in the file
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

        if sys.version_info[0] == 2 :
            mol = pybel.readfile( format = format, filename = infile ).next()
        elif sys.version_info[0] == 3 :
            mol = pybel.readfile( format = format, filename = infile ).__next__()
        else : raise Exception( "golang is looking better every day" )
        return cls( mol, verbose )

    #
    #
    def __init__( self, molecule, verbose = False ) :
        assert isinstance( molecule, pybel.Molecule )
        self._mol = molecule
        self._verbose = verbose

    #
    #
    def __sdf__( self ) :
        if self._verbose : sys.stdout.write( "%s.__sdf__()\n" % (self.__class__.__name__,) )
        return self._mol.write( "sdf" )
    sdf = property( __sdf__ )

    #
    #
    @property
    def version( self ) :
        return str( openbabel.OBReleaseVersion() ).strip()

    #
    #
    @property
    def program( self ) :
        return "OpenBabel"

    # different toolkits can generate different smiles strings. OTOH they tend to generate the same
    # strings for "canonical", "isomeric", and "just smiles". Go figure.
    #
    # openbabel smiles string has trailing garbage
    #
    def to_smiles( self, kind = None ) :
        if self._verbose : sys.stdout.write( "%s.to_smiles(%s)\n" % (self.__class__.__name__,kind,) )
        rc = None
        if kind == "canonical" :
            rc = self._mol.write( format = "can" )
        else :
            rc = self._mol.write( format = "smi" )
        if rc is not None :
            return rc.strip().split()[0]
        return None

    def __smiles__( self ) :
        rc = []
        can = self.to_smiles( kind = "canonical" )
        if can is not None :
            rc.append( { "descriptor" : can, "type" : "SMILES_CANONICAL", "program" : self.program,
                        "version" : self.version } )
        smi = self.to_smiles()
        if smi is not None :
            rc.append( { "descriptor" : can, "type" : "SMILES", "program" : self.program,
                        "version" : self.version } )
        if len( rc ) < 1 :
            rc = None
        return rc

    smiles = property( __smiles__ )

    # InChI should be the same between all toolkits as long as nobody messes with the molecule
    #
    def __inchi__( self ) :
        inchi = self._mol.write( format = "inchi" )
        key = self._mol.write(  format = "inchikey" )
        return [ { "type" : "InChI", "descriptor" : inchi.strip(), "program" : self.program, "version" : self.version },
            { "type" : "InChI_KEY", "descriptor" : key.strip(), "program" : self.program, "version" : self.version } ]
    inchi = property( __inchi__ )

    # chem comp stuff formatted as json
    # can get big-ish, don't use on large molecules
    #
    def __json__( self ) :

        rc = {}
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

    # masses
    #
    def monoisotopic_mass( self, carbon = None, nitrogen = None ) :
        if (carbon is None) and (nitrogen is None) :
            return self._mol.exactmass

# openbabel's copy-constructor does not preserve atom and bond indices, but for this we don't care
# it's not clear that modifying the molecule "in-place" won't break something else
#
        mol = openbabel.OBMol( self._mol.OBMol )
        mol.BeginModify()
        for atom in openbabel.OBMolAtomIter( mol ) :
#            if atom.IsCarbon() :
            if atom.GetAtomicNum() == openbabel.Carbon :
                if carbon is not None :
                    atom.SetIsotope( int( carbon ) )
            elif atom.GetAtomicNum() == openbabel.Nitrogen :
                if nitrogen is not None :
                    atom.SetIsotope( int( nitrogen ) )
        mol.EndModify()
        return mol.GetExactMass()

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

        return re.sub( r"<title>.+</title>", ("<title>%s</title>" % (t,)), rc )

    # 'presentation' pix has only backbone atoms and no labels.
    # this returns an SVG with transparent background and dimenstions massaged to 20em high, no width
    #
    # Note that this changes the molecule: it removes protons. Adding them back is not guaranteed to
    # result in the same molecule (stereo chemistry bites)
    #
    def _to_svg_small( self, title = None ) :
        """return 'presentation' SVG"""
        if self._verbose : sys.stdout.write( "%s._to_svg_small()\n" % (self.__class__.__name__,) )

        if (title is not None) and (str( title ).strip() != "" ) :
            t = str( title ).strip()
        else :
            t = ""
        self._mol.OBMol.SetTitle( t )

        self._mol.removeh()

        svg = self._mol.write( format = "svg", opt = { "b" : "none" } )
        svg = svg.replace( 'width="200px"', '' ).replace( 'height="200px"', 'height="20em"' )

        return svg

    # 'full' pix has all atoms with labels.
    # this returns an SVG with transparent background and dimenstions massaged to 1000x1000px,
    # and more readable colour and size of atom labels
    #
    def _to_svg_full( self, title = None ) :
        """return 'full' SVG"""
        if self._verbose : sys.stdout.write( "%s._to_svg_full()\n" % (self.__class__.__name__,) )

        if (title is not None) and (str( title ).strip() != "" ) :
            t = str( title ).strip()
        else :
            t = ""
        self._mol.OBMol.SetTitle( t )

        svg = self._mol.write( format = "svg", opt = { "a" : None, "b" : "none", "i" : None } )

        svg = svg.replace( 'width="200px"', 'width="1000px"' ).replace( 'height="200px"', 'height="1000px"' )

# atom indices: change colour
#
        clrpat = re.compile( r'(<text\s.+fill=")rgb\(255,0,0\)("\s+stroke=")rgb\(255,0,0\)(".+>\d+</text>)' )
#        svg = re.sub( clrpat, r"\1rgb(0,135,0)\2rgb(0,135,0)\3", svg )
        svg = re.sub( clrpat, r"\1rgb(0,135,0)\2rgb(255,255,0)\3", svg )

# and font
#
        fntpat = re.compile( r'(<text\s.+font-size=)"12"(\s*>\d+</text>)' )
        svg = re.sub( fntpat, r'\1"8"\2', svg )

        strpat = re.compile( r'(<text\s.+stroke-width=)"1"(.+>\d+</text>)' )
        svg = re.sub( strpat, r'\1"0.3"\2', svg )

        cpat = re.compile( r'(<text\s.+fill=")rgb\(102,102,102\)("\s+stroke=")rgb\(102,102,102\)(".+>C</text>)' )
        svg = re.sub( cpat, r"\1rgb(0,0,0)\2rgb(0,0,0)\3", svg )

        hpat = re.compile( r'(<text\s.+fill=")rgb\(191,191,191\)("\s+stroke=")rgb\(191,191,191\)(".+>H</text>)' )
        svg = re.sub( hpat, r"\1rgb(100,100,100)\2rgb(100,100,100)\3", svg )

        return svg

    # making molecule images, 3D in particular, is a semi-manual process: you have to try a few
    # tricks and a few different programs to make it look good.
    #
    # energy minimization is one of them. Note that it changes the molecule.
    #
    def minimize( self, field = "MMFF94" ) :
        """minimize energy using specified force field"""
        if self._verbose : sys.stdout.write( "%s.minimize()\n" % (self.__class__.__name__,) )

        ff = openbabel.OBForceField.FindForceField( field )
        if ff == 0 :
            if self._verbose : sys.stdout.write( "No force field %s!\n" % (field,) )
            return

        ff.SetLogLevel( openbabel.OBFF_LOGLVL_NONE ) # LOW )
        ff.SetLogToStdErr()

        if ff.Setup( self._mol.OBMol ) == 0 :
            if self._verbose : sys.stdout.write( "Can't setup force field %s!\n" % (field,) )
            return

        ff.SteepestDescent( 2000 )
        ff.GetCoordinates( self._mol.OBMol )

    # chem comp bond
    # return { "id" : N, "type" : covalent|amide|ester|carbonyl, "atom1" : I, "atom2" : J }
    #
    # RDKit does not provide bond type: carbonyl, ester, etc. that I know of
    # OpenBabel does not do E|Z stereo that I know of
    #
    # all bonds are covalent by default
    #
    # atom1, atom2 should match between the 2 provided nobody messed with the molecule before calling
    # these iterators
    #
    def iter_bonds( self ) :
        """bonds iterator"""
        if self._verbose : sys.stdout.write( "%s.iter_bonds()\n" % (self.__class__.__name__,) )

#        et = openbabel.OBElementTable()
        for i in range( self._mol.OBMol.NumBonds() ) :
            bond = self._mol.OBMol.GetBond( i )
            rc = { "id" : (bond.GetIdx() + 1) }
            rc["type"] = "covalent"
            if bond.IsAmide() : rc["type"] = "amide"
            elif bond.IsEster() : rc["type"] = "ester"
            elif bond.IsCarbonyl() : rc["type"] = "carbonyl"

#            elt1 = et.GetSymbol( bond.GetBeginAtom().GetAtomicNum() )
            elt1 = openbabel.GetSymbol( bond.GetBeginAtom().GetAtomicNum() )
            if str( elt1 ).lower() == "xx" : elt1 = "x"
            rc["atom1"] = "%s%d" % (elt1,bond.GetBeginAtom().GetIdx(),)
#            elt2 = et.GetSymbol( bond.GetEndAtom().GetAtomicNum() )
            elt2 = openbabel.GetSymbol( bond.GetEndAtom().GetAtomicNum() )
            if str( elt2 ).lower() == "xx" : elt2 = "x"
            rc["atom2"] = "%s%d" % (elt2,bond.GetEndAtom().GetIdx(),)

            yield rc

    # extra stuff in SDF
    #
    def iter_sdfdata( self ) :
        """iterate over extra properties extracted from input (SDF)"""
        return self._mol.data.items()

#
#
#
if __name__ == '__main__':
    m = OBmolecule.from_file( sys.argv[1], verbose = True )
#    m.minimize()
#    print m.sdf
#    print
#    print m._to_svg_full()
    sys.stdout.write( m.json )

#
#
