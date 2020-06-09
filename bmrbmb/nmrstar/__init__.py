#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#
#  __init__.py
#
#  Copyright 2018 Board of Regents university of Wisconsin - Madison
#    Dimitri Maziuk <dmaziuk@bmrb.wisc.edu>
#
#  This code is free: reuse what you like but give credit
#
# Add bits needed for a BMRB entry: assembly, entity
# Create entry directory structure and copy/rename files
#

from __future__ import absolute_import

import os
import sys
import sqlite3
import re

README = """
Directory contents:
-------------------
<BMRB ID>.sdf         -- SDF: source for the molecule info, usually downloaded from PubChem.
<BMRB ID>.mol         -- MOL version of the above.
<BMRB ID>.str         -- NMR-STAR file for the entry. NMR_Experiment_File and Experiment tables
                         link data directories to experiment details such as sample, conditions,
                         spectrometer information, etc.
<BMRB ID>_nom.svg     -- 2D molecule image with atom labels (matching atom names in the NMR-STAR file
                         and the order in SDF/MOL files).
<BMRB ID>.png         -- 3D molecule image
<BMRB ID>.svg         -- 2D molecule image

nmr                   -- NMR experiment data.
nmr/setNN             -- "set01" should always be there. E.g. concentration or pH study may contain
                         more that one set of NMR data.
                         Subdirectories are in TopSpin format: "1", "2", etc.
nmr/setNN/NMRbot*.txt -- NMRbot log file.
nmr/setNN/spectra     -- spectral images. If more than one, image named *_0 is the full spectrum.
nmr/setNN/transitions -- Peak lists.
"""

# datasets are indexed 1, 2, ... in the input
# this is more human-readable.
# NOTE only 99 datasets/entry.
#
DATASET_MASK = "set%02d"
DATASET_DIR = os.path.join( "nmr", DATASET_MASK )
PEAKFILE_DIR = os.path.join( DATASET_DIR, "transitions" )
SPECTRA_DIR = os.path.join( DATASET_DIR, "spectra" )

# sqlite3 database with one table: ids (id text primary key, dirname text, inchi text)
#
IDLIST = "/share/dmaziuk/projects/metabolomics_pipe/bmrbids.sqlt3"

# return bmrb id for dirname. inchi string gets added when creating a new record
#
def bmrbid( dirname, inchi, verbose = False ) :
    global IDLIST
    if not os.path.exists( os.path.realpath( IDLIST ) ) :
        raise IOError( "Not found: %s" % (os.path.realpath( IDLIST ),) )

    conn = sqlite3.connect( os.path.realpath( IDLIST ) )
    curs = conn.cursor()
    rc = None
    sql = "select id from ids where dirname=?"
    curs.execute( sql, (dirname,) )
    for row in curs :
        rc = row[0]
        break

    if rc is None :
        sql = "select id from ids where inchi=?"
        curs.execute( sql, (inchi,) )
        for row in curs :
            sys.stderr.write( "WARN: no %s in ID list, but %s has ID %s\n" % (dirname,inchi,row[0],) )

# make new one
#
        pat = re.compile( r"^bmse(\d+)$" )
        sql = "select id from ids"
        bid = -1
        curs.execute( sql )
        for row in curs :
            m = pat.search( row[0] )
            if not m :
                raise Exception( "ID %s does not match pattern" % (row[0],) )
            if int( m.group( 1 ) ) > bid :
                bid = int( m.group( 1 ) )

        if bid < 0 :
            raise Exception( "Can't add BMRB ID" )

        bid += 1
        rc = "bmse%06d" % (bid,)

        sql = "insert into ids (id,dirname,inchi) values (?,?,?)"
        curs.execute( sql, (rc,dirname,inchi,) )
        conn.commit()

    curs.close()
    conn.close()
    return rc

#
#
#
from .entrydir import EntryDir
from .extras import EntryExtras
from .star import StarMaker

#
#
#
__all__ = [ "README", "DATASET_MASK", "DATASET_DIR", "PEAKFILE_DIR", "SPECTRA_DIR",
        "bmrbid",
        "EntryDir",
        "EntryExtras",
        "StarMaker"
    ]
#
#
