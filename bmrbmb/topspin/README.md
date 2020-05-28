## process files out of Bruker TopSpin

`botlog.py` - this one is actually for processig `NMRBot` log file. The log contains NMR experiment list
and some parameters: sweep widths etc. (that could be pulled out of TopSpin files but it was easier to
write them out from NMRBot).

`exptfiles.py` - a wrapper for `tar`/`gzip` that extracts timedomain data maing sure the paths are OK etc.

`topspinpeaks.py` - reads TopSpin `peaklist.xml` into NMR-STAR-like peak list structure(s). Peaks go into 
NMR-STAR "transition" list. If "annotation" field is present, it's assumend to contain atom name (convention
that NMRFAM spectroscopist follows), and "assigned transition" list is also created. 

*BUG*: TopSpin sometimes generates `peaklist.xml` with multiple peak list sections. So far I've only seen
them with one section containing data, the others are empty, but I've no idea if that is always true.
This code will error out on files like that, edit the file and delete the extra peak list sections.
