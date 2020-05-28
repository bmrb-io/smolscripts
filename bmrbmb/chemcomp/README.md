## makes "chemical component" and related stuff

Uses `openbabel` and `rdkit` to make the chem. comp. Uses `pymol` to render the "3D" image.

`molecule.py` is the wrapper around `obmol.py` (OpenBabel) and `rdmol.py` (RDKit) wrappers.
Originally I needed both: OpenBabel for SVG depictions and RDKit for CIP stereochemistry.
With RDKit's recent drawing code I could deprecate OpenBabel altogether but OpenBabel is
easy. It's RDKit that's a pain to install.

`img3d.py` is a PyMol wrapper that makes a pretty "3D" picture of the molecule. Worked with 
open-source pymol last I tried (but we mainly used University-licensed commercial version).
