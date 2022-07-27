#!/bin/sh

./clean.sh

echo "** Creating input files for large model **"
# Creating the input rss files
python3 large.py


echo "** Start cpu modelling for large model **"
# Running single node acoustic gpu modeling
rsAcousticmod2dgpu mod-large.cfg

# Plotting seismograms
python3 Ploting_seismo.py

