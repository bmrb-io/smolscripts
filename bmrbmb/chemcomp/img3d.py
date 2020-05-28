#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  img3d.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit

# this is a wrapper for pymol
# you need pymol installed and in PYMOL_PATH in __init__.py
#
from __future__ import absolute_import

import os
import sys
import time
import threading

PYMOL_PATH = "/opt/pymol/lib/python2.7/site-packages"
sys.path.append( PYMOL_PATH )
import pymol

# Sometimes some of the commented out lines work, and sometimes: better than not-coommented-out ones.
# Sometimes it generates a file and sometimes the file is zero bytes.
# Go figure.
#
def make_image( infile, outfile, verbose = False ) :
    molfile = os.path.realpath( infile )
    if not os.path.exists( infile ) :
        raise IOError( "File not found: %s" % (molfile,) )

    pngfile = os.path.realpath( outfile )

    pymol.pymol_argv = ["pymol", "-qc"]
    pymol.finish_launching()

    cmd = pymol.cmd
    cmd.load( molfile )
    cmd.hide( "everything" )

# sticks or lines: choose wisely
#
#    cmd.show( "sticks" )
#    cmd.set_bond( "stick_radius", 0.25, "all" )
#    cmd.set( "stick_radius", 0.25 )

    cmd.show( "lines" )
    cmd.set( "line_as_cylinders", 1 )
    cmd.set( "line_width", 8.0 )

#    cmd.show( "cartoon" )
#    cmd.set( "cartoon_discrete_colors", 1 )
#    cmd.set( "cartoon_fancy_helices", 1 )
#    cmd.set( "cartoon_side_chain_helper", "on" )

#    cmd.color( "red", "ss h" )
#    cmd.color( "yellow", "ss s" )
#    cmd.color( "green", "ss l+''" )

    cmd.util.cbaw()

    cmd.set( "ray_opaque_background", "off" )
    cmd.set( "ray_trace_mode",  1 )
    cmd.set( "antialias", 2 )
    cmd.set( "ray_trace_color", "grey" )

#    cmd.set( "nonbonded_as_cylinders", 1 )
#    cmd.set( "ribbon_as_cylinders", 1 )
#    cmd.set( "mesh_as_cylinders", 1 )
#    cmd.set( "dash_as_cylinders", 1 )
#    cmd.set( "render_as_cylinders", 1 )

    cmd.set( "orthoscopic", "off" )
    cmd.set( "valence", "off" )

    cmd.png( pngfile, width = 800, dpi = 300, ray = 1 )

# pymol fires off rendering in a background thread that doesn't finish when the above returns
#  and then you still don't know if it's done or not
#
    while threading.active_count() > 2 :
        time.sleep( 2 )
    cmd.quit()

    time.sleep( 5 )

    if not os.path.exists( pngfile ) :
        sys.stderr.write( "Failed: no file (%s)\n" % (pngfile,) )
        return False

    si = os.stat( outfile )
    if si.st_size < 1 :
        sys.stderr.write( "Failed: 0-byte file (%s)\n" % (pngfile,) )
        os.unlink( outfile )
        return False


#
#
#
if __name__ == '__main__':
    make_image( sys.argv[1], "out.png", verbose = True )

#
#
