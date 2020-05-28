#!/usr/bin/python -u
#

import sys
import os

_UP = os.path.realpath( os.path.join( os.path.split( __file__ )[0], ".." ) )
sys.path.append( _UP )
import bmrbmb

if __name__ == "__main__" :

    srcdir = os.path.realpath( sys.argv[1] )
    inchi = None
    with open( os.path.join( srcdir, "inchi.inchi" ), "rU" ) as fin :
        inchi = fin.read()
    if inchi is None : raise Exception( "no inchi" )

    idx = { "sdf" : os.path.join( srcdir, "alatis_output_compact.sdf" ),
            "paramagnetic" : "no",
            "inchi" : inchi }

    cnv = bmrbmb.incoming2json.FAMtoJSN( index = idx, workdir = srcdir )
    cnv._make_chem_comp()
    sys.stdout.write( cnv.__json__ )

#
#
