#!/bin/sh

./clean.sh

echo "** Creating input files for small model **"
# Creating the input rss files
python3 small.py

echo "** Start gpu modelling for small model **"
# Running single node acoustic gpu modeling
/usr/bin/time --format='Wall time: %e sec(s)' rsAcousticmod2dgpu mod-small.cfg

# Plotting seismograms
python3 Ploting_seismo.py

