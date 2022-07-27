#!/bin/sh

./clean.sh

echo "** Creating input files for small model **"
# Creating the input rss files
python3 small.py

echo "** Start cpu modelling for small model **"
# Running single node acoustic gpu modeling
rsAcousticmod2dgpu mod-small.cfg

# Plotting seismograms
python3 Ploting_seismo.py

# Plotting snapshots
python3 Ploting_snaps-single.py
