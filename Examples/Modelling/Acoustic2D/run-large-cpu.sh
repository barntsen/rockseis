#!/bin/sh

./clean.sh

echo "** Creating input files for large model **"
# Creating the input rss files
python3 large.py


echo "** Start cpu modelling for large model **"
# Running single node acoustic cpu modeling
/usr/bin/time --format='Wall time: %e sec(s)' rsAcousticmod2d mod-large.cfg


# Plotting seismograms
python3 Ploting_seismo.py

