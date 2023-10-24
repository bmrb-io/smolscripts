#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate rdkit-env

/opt/scripts/bmrbmb/schwalbe2json.py $1 $2 $3
