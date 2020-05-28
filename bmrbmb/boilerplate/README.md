## These are boilerplate values for BMRB metabolomics entries.

Each file is a JSON list of dicts corresp. to NMR-STAR saveframes.
There may be several related saveframes in one file, e.g. software
and citation for that software.

Each dict has name and "sf_category". It may have an "id" -- the 
"id" is special: when there's a choice between several saveframes
(e.g. choose DSS or TMS referencing), that is the key for picking
one.

  * `key - single value` pairs correspond to NMR-STAR tags at the 
saveframe level.

  * `Key - list` values correspond to NMR-STAR loops (tables). Tables are
lists of rows, rows are dicts of tag (column) - value pairs. 
(This is going into an SQL backend, row-wise is easier for inserts.)
