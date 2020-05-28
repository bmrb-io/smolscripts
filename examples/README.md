
`alatis_output_compact.*` and `to_json.py` shows how to use the code to make an NMR-STAR chem comp 
(as JSON) from SDF.

`lovastatin.str` is the output NMR-STAR file for BMRB entry BMSE001335 (just the STAR file, complete entry
directory can be found on the BMRB website), `lovastatin` subdir is the partial input: there's little point in 
keeping a 100MB timedomain data tarball in the repo, so that's been trimmed down to a single experiment.
