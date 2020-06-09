#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  extras.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#
# extra stuff:
#  a "3D" PNG (pymol)
#  a MOL file
#  a "presenttion" PNG


from __future__ import absolute_import

import os
import sys
import tempfile

# self
#
#_HERE = os.path.split( __file__ )[0]
#sys.path.append( os.path.realpath( os.path.join( os.path.join( _HERE, ".." ), ".." ) ) )
#import bmrbmb.chemcomp

from .. import chemcomp

# class because I want to keep the molecule
#
#
class EntryExtras( object ) :

    #
    #
    @classmethod
    def makefiles( cls, sdf, verbose = False ) :
        rc = cls( sdf = sdf, verbose = verbose )
        rc.make_mol()
        rc.make_svg()
        rc.make_png()

    #
    #
    def __init__( self, sdf, verbose ) :
        infile = os.path.realpath( sdf )
        if not os.path.exists( infile ) :
            raise IOError( "Not found: %s" % (infile,) )
        self._infile = infile
        self._mol = chemcomp.Molecule.from_file( filename = infile, verbose = verbose )
        self._verbose = bool( verbose )

    # SDF to MOL: this only really strips off SDF "properties" from the end of MOL block.
    # SDF can contain more than one MOL (ours don't), apparently not all software out there can handle that.
    #
    def make_mol( self, outfile = None ) :
        if outfile is None :
            outfile = "%s.mol" % (os.path.splitext( self._infile )[0],)

        if self._infile == outfile :
            raise Exception( "%s: input and output are the same file" )

        with open( outfile, "wb" ) as out :
            out.write( self._mol.mol )

    # SVG is a "presentation" image to be used on webpages etc.
    # openbabel usually makes better pictures, though for "presentation" ones it rarely matters
    #
    def make_svg( self ) :

        outfile = "%s.svg" % (os.path.splitext( self._infile )[0],)
        with open( outfile, "wb" ) as out :
            out.write( self._mol._obmol.to_svg( full = False ) )

    # pymol: MOL to 3D PNG. The image is intended to be used instead of jsmol widget in lists etc.
    # where you'd end up creating lots of jsmol instances.
    # energy-minimizing first tends to help, RDKit does a better job of it
    #
    def make_png( self, minimize = True ) :

        if minimize :
            (fd, fname) = tempfile.mkstemp( dir = os.path.split( self._infile )[0], suffix = ".mol")
            molfile = fname
            self._mol._rdmol.minimize()
            os.write( fd, self._mol._rdmol.mol )
            os.close( fd )
        else :
            molfile = "%s.mol" % (os.path.splitext( self._infile )[0],)

        outfile = "%s.png" % (os.path.splitext( self._infile )[0],)
        chemcomp.make_image( molfile, outfile, verbose = self._verbose )

        if minimize :
            os.unlink( molfile )

#
#
#
if __name__ == '__main__':
    EntryExtras.makefiles( sdf = sys.argv[1], verbose = True )

#
#
