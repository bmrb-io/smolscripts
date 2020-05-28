## deal with assignment table

This is very much BMRB-NMRFAM-specific: the spectroscopist creates assigned chemical shifts table
in this format (and no other):
```
atom_ID,chem_shift,ambiguity,set_ID
C1,24.9588,1,
C2,"63.1833, 63.6985, 63.7526",4,1
C3,"63.1833, 63.6985, 63.7526",4,1
C4,"63.1833, 63.6985, 63.7526",4,1
C5,62.7662,1,
```

Code here turns it into NMR-STAR assignments and ambiguity tables. The tricky part is NMR-STAR requires one
chemical shift assigned to one IUPAC atom name. IRL there's _n_ peaks for _m_ atoms. 

  * If _n_ > _m_, e.g. for mixtures of isomers, like sugars, NMR-STAR requires the molecule to be defined 
completely differently. This code will just error out.

  * If _n_ <= _m_ (peaks overlapping) the code will simply repeat values from _n_ for each next atom.
It's up to the user to figure out that any value from _n_ applies to any atom in _m_ and the one-to-one
correspondence is an artifact of the data model.

