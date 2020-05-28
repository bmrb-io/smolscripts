# smolscripts

Internal: processing scripts for creating metabolomics/small molecule BMRB entries from NMR facility data.

Some of this code is of more general use:
  * create BMRB entry from (carefully crafted) JSON,
  * create NMR-STAR molecule sections from SDF/MOL,
  * create molecule images from SFD/MOL. incl. WSGI web insterface for that.

## require:

 - pymol
 - rdkit
 - obabel

## files

`assign_bmseid.py and bmrbids[.*].sql|sqlt3` - table (sqlite3) of BMSE IDs and associated InChI strings, and a script to look up/append new one

`Dockerfile` builds a container for use @ NMRFAM
  - conda install of uwsgi, python 3.6, werkzeug (flask), rdkit, and openbabel
  - runs `wsgi.py` in uwsgi on port 9090


