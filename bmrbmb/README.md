At this level, processing is a 2-step shuffle:

  1. Read/check input files and move them to "sessiondir"/incoming, renaming for uniformity.
  2. Make data structure from files in incoming. It follows NMR-STAR data model and can be dumped to JSON.

## files

`nmrfam2incoming.py` - step 1

`incoming2json.py` - step 2

`nmrfam2json.py.old` - initial version of the code (I'm too lazy to preserve local version history when migrating to github, I'll just comit the old file and remove it later)
